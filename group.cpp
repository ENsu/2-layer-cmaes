#include "group.h"
#include "global.h"
#include <assert.h>
#include <cfloat>

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
       if(iter->getFitness() > tmpMax)
           tmpMax = iter->getFitness();
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

double Group::getMean(double *weight, int weight_length)
{
    assert(weight_length == Nodes.size());
    // make sure weight and Nodes have same size
    double Sum = 0;
    double weightSum = 0;

    int index = 0;
    list<Node>::iterator iter;
    for(iter = Nodes.begin() ; iter != Nodes.end() ; ++iter)
    {
        Sum += weight[index] * iter->getFitness();
        weightSum += weight[index];
        index++;
    }
    assert(weightSum == 1);
    return Sum;

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

Eigen::MatrixXd Group::getCov()
{
    Eigen::MatrixXd Cmu;
    Cmu.setZero(dimension, dimension);
    Eigen::VectorXd tmp;

    list<Node>::iterator iter;
    for(iter = Nodes.begin() ; iter != Nodes.end() ; ++iter)
    {
        tmp = iter->allele;
        Cmu += tmp*tmp.transpose();
    }
    return Cmu;
}

Eigen::MatrixXd Group::getCov(double *weight, int weight_length)
{
    assert(weight_length == Nodes.size());
    // make sure weight and Nodes have same size
    Eigen::MatrixXd Cmu;
    Cmu.setZero(dimension, dimension);
    Eigen::VectorXd tmp;

    list<Node>::iterator iter;
    int index = 0;
    for(iter = Nodes.begin() ; iter != Nodes.end() ; ++iter)
    {
        tmp = iter->allele;
        Cmu += (tmp*tmp.transpose() * weight[index]);
        index ++;
    }
    return Cmu;
}

double Group::getUCBVal(int total_played, double global_min, double global_max)
{
    //ucb_tuned1
    double v = getNormVariance(global_min, global_max) + sqrt(2*log(total_played) / Nodes.size());
    cout << getNormVariance(global_min, global_max) << endl;
    cout << sqrt(2*log(total_played) / Nodes.size()) << endl;
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
    mean.fitness = mean.getFitness();
    return mean;
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

void Group::print()
{
    list<Node>::iterator iter;
    cout << "========== group nodes: ===========" << endl;
    for(iter = Nodes.begin() ; iter != Nodes.end() ; ++iter)
        iter->print();
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
