#include "group.h"
#include "global.h"
#include <assert.h>
#include <iostream>
#include <cfloat>
#include "random.hpp"

extern randomG RANDOM;

Group::Group(list<Node> new_Nodes)
{
    Nodes.clear();
    Nodes = new_Nodes;
    AllNodeNum = Nodes.size();
}

Group::Group()
{
    Nodes.clear();
    AllNodeNum = 0 ;
}

Group::~Group()
{
    Nodes.clear();
}

double Group::getMin()
{
    double tmpMin = DBL_MAX;
    list<Node>::iterator iter;
    for(iter = Nodes.begin() ; iter != Nodes.end() ; ++iter)
       if(iter->getFitness() < tmpMin)
           tmpMin = iter->getFitness();
    return tmpMin;
}

double Group::getMax()
{
    double tmpMax = -1* DBL_MAX;
    list<Node>::iterator iter;
    for(iter = Nodes.begin() ; iter != Nodes.end() ; ++iter)
    {
       if(iter->getFitness() > tmpMax)
           tmpMax = iter->getFitness();
    }
    return tmpMax;
}

double Group::getMean()
{
    double tmpSum = 0;
    list<Node>::iterator iter;
    for(iter = Nodes.begin() ; iter != Nodes.end() ; ++iter)
           tmpSum += iter->getFitness();
    return (tmpSum / Nodes.size());
}

double Group::getVariance()
{
    double mean = getMean();
    double sum2 = 0;
    list<Node>::iterator iter;
    for(iter = Nodes.begin() ; iter != Nodes.end() ; ++iter)
           sum2 += iter->getFitness() * iter->getFitness();
    return (sum2 / Nodes.size() - mean * mean);
}

double Group::getNormVariance(double global_min, double global_max)
{
    double norm_mean = (getMean() - global_max) / (global_min - global_max);
    double sum2 = 0;
    list<Node>::iterator iter;
    for(iter = Nodes.begin() ; iter != Nodes.end() ; ++iter)
    {
           double tmp = (iter->getFitness() - global_max) / (global_min - global_max);
           sum2 +=  tmp * tmp;
    }
    return (sum2 / Nodes.size() - norm_mean * norm_mean);
}

double Group::getUCBVal(int total_played, double global_min, double global_max)
{
    //ucb_tuned1
    double v = getNormVariance(global_min, global_max) + sqrt(2*log(total_played) / Nodes.size());
    double V = (v>0.25) ? 0.25 : v;
    double norm_mean = (getMean() - global_max) / (global_min - global_max);
    // double norm_min = (getMin() - global_max) / (global_min - global_max);

    // double tmp = newMean - sqrt(log(total) / AllUsedNumber *V);
    return norm_mean + sqrt(V * log(total_played) / Nodes.size());
}

int Group::getSize()
{
    return Nodes.size();
}

Node Group::get_mean_node()
{
    Node tmp;
    list<Node>::iterator iter;
    for(iter = Nodes.begin() ; iter != Nodes.end() ; ++iter)
	   tmp.allele = tmp.allele + (*iter).allele;
    Node mean;
    mean.allele = tmp.allele / Nodes.size();
    assert(!mean.outofBound());
    return mean;
}

Node Group::get_mean_node(Eigen::VectorXd weight)
{
    // make sure weight and Nodes have same size 
    assert(weight.size() == Nodes.size());
    double weightSum = 0;
    int index = 0;
    sort_node();

    Node *mean = new Node();
    mean->allele.setZero();
    list<Node>::iterator iter;
    for(iter = Nodes.begin() ; iter != Nodes.end() ; ++iter)
    {
       mean->allele = mean->allele + (*iter).allele * weight[index];
       weightSum += weight[index];
       index ++;
    }
    if(mean->outofBound())
    {
        mean->intoBound(); 
        //這真的是很必不得已的作法，mean有得時候會算錯，然後跑到boundary外面
        //其實算錯的值都在~1e16，但是因為那時候sigma*D更小，所以就會卡住sample always outof bound
    }
    assert(fabs(weightSum - 1) == 0);
    assert(!mean->outofBound());
    return (*mean);
}

Node Group::get_best_node()
{
    Node best;
    list<Node>::iterator iter = Nodes.begin();
    double min = iter->getFitness();
    best = *iter;
    for( ; iter != Nodes.end() ; ++iter)
    {
        if(iter->getFitness() < min)
        {
            min = iter->getFitness();
            best = *iter;
        }
    }
    return best;
}

Node Group::get_worst_node()
{
    Node worst;
    list<Node>::iterator iter = Nodes.begin();
    double max = iter->getFitness();
    worst = *iter;
    for( ; iter != Nodes.end() ; ++iter)
    {
        if(iter->getFitness() > max)
        {
            max = iter->getFitness();
            worst = *iter;
        }
    }
    return worst;
}

Eigen::MatrixXd Group::node_matrix()
{
    sort_node();
    Eigen::MatrixXd tmp(dimension, Nodes.size());
    int index = 0;
    list<Node>::iterator iter;
    for(iter = Nodes.begin() ; iter != Nodes.end() ; ++iter)
    {
        tmp.col(index) = iter->allele;
        index ++;
    }
    return tmp;
}

void Group::add_node(Node a)
{
    Nodes.push_back(a);
    AllNodeNum ++;
    return ;
}

void Group::drop_node(Node a)
{
    list<Node>::iterator iter;
    for(iter = Nodes.begin() ; iter != Nodes.end() ; ++iter)
        if(iter->allele == a.allele)
            Nodes.erase(iter);
    return ;
}


void Group::replace_node(Node a, Node b)
{
    add_node(a);
    drop_node(b);
}

void Group::truncate_size(int size)
{
    list<Node>::iterator iter_st = Nodes.begin();
    list<Node>::iterator iter_end = Nodes.end();
    advance(iter_st, size);   
    Nodes.erase(iter_st, iter_end);
}

bool compare_node (Node& first, Node& second)
{
  return ( first.getFitness() < second.getFitness() );
}

void Group::sort_node()
{
    Nodes.sort(compare_node);
}

list<Node> Group::random_pick(int pick_num)
{
    if(pick_num > Nodes.size())
        return Nodes;
    else
    {
        list<Node> pick_nodes;
        int *intArry = new int[Nodes.size()];
        RANDOM.uniformArray(intArry , Nodes.size(), 0 , Nodes.size()-1);
        for(int i=0; i< pick_num; i++)
        {
            list<Node>::iterator iter = Nodes.begin();
            std::advance(iter, intArry[i]);
            pick_nodes.push_back(*iter);
        }
        return pick_nodes;
    }
}

void Group::print()
{
    list<Node>::iterator iter;
    cout << "========== group nodes: ===========" << endl;
    for(iter = Nodes.begin() ; iter != Nodes.end() ; ++iter)
    {
        iter->print();
    }
    cout << "===================================" << endl;
}

Group& Group::operator=(const Group rhs)
{
    if(this == &rhs)
        return *this;

    AllNodeNum = rhs.AllNodeNum;
    Nodes.clear();
    Nodes = rhs.Nodes;
    return *this;
}
