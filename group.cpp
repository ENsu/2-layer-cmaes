#include "group.h"
#include <assert.h>
#include <cfloat>


GROUP::GROUP(int new_gid, int dimension, list<Node> new_Nodes)
{
    gID = new_gid;
    dim = dimension;
    Nodes = new_Nodes;
    AllNodeNum = Nodes.size();
}

GROUP::GROUP(int new_gid, int dimension)
{
    gID = new_gid;
    dim = dimension;
    Nodes.clear();
    AllNodeNum = 0 ;
}

GROUP::~GROUP()
{
    Nodes.clear();
}

double GROUP::getMin()
{
    double tmpMin = DBL_MAX;
    list<Node>::iterator iter;
    for(iter = Nodes.begin() ; iter != Nodes.end() ; ++iter)
       if(iter->getFitness() < tmpMin)
           tmpMin = iter->getFitness();
    return tmpMin;
}

double GROUP::getMax()
{
    double tmpMax = -1* DBL_MAX;
    list<Node>::iterator iter;
    for(iter = Nodes.begin() ; iter != Nodes.end() ; ++iter)
       if(iter->getFitness() > tmpMax)
           tmpMax = iter->getFitness();
    return tmpMax;
}

double GROUP::getMean()
{
    double tmpSum = 0;
    list<Node>::iterator iter;
    for(iter = Nodes.begin() ; iter != Nodes.end() ; ++iter)
           tmpSum += iter->getFitness();
    return (tmpSum / Nodes.size());
}

double GROUP::getMean(double *weight)
{
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
    return (Sum / weightSum);

}

double GROUP::getVariance()
{
    double mean = getMean();
    double sum2 = 0;
    list<Node>::iterator iter;
    for(iter = Nodes.begin() ; iter != Nodes.end() ; ++iter)
           sum2 += iter->getFitness() * iter->getFitness();
    return (sum2 / Nodes.size() - mean * mean);
}

double GROUP::getNormVariance(double global_min, double global_max)
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

double GROUP::getUCBVal(int total_played, double global_min, double global_max)
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

int GROUP::getSize()
{
    return Nodes.size();
}

Node GROUP::get_mean_node()
{
    Node tmp(dim);
    list<Node>::iterator iter;
    for(iter = Nodes.begin() ; iter != Nodes.end() ; ++iter)
	   tmp.allele = tmp.allele + (*iter).allele;
    Node mean(dim);
    mean.allele = tmp.allele / Nodes.size();
    mean.fitness = mean.getFitness();
    return mean;
}

Node GROUP::get_best_node()
{
    Node best(dim) ;
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

Node GROUP::get_worst_node()
{
    Node worst(dim) ;
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

void GROUP::add_node(Node a)
{
    Nodes.push_back(a);
    AllNodeNum ++;
    return ;
}

void GROUP::drop_node(Node a)
{
    list<Node>::iterator iter;
    for(iter = Nodes.begin() ; iter != Nodes.end() ; ++iter)
        if(iter->allele == a.allele)
            Nodes.erase(iter);
    return ;
}


void GROUP::replace_node(Node a, Node b)
{
    add_node(a);
    drop_node(b);
}

bool compare_node (Node& first, Node& second)
{
  return ( first.getFitness() > second.getFitness() );
}

void GROUP::sort_node()
{
    Nodes.sort(compare_node);
}

void GROUP::print()
{
    list<Node>::iterator iter;
    cout << "========== group nodes: ===========" << endl;
    for(iter = Nodes.begin() ; iter != Nodes.end() ; ++iter)
        iter->print();
    cout << "===================================" << endl;
}
