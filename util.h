#include <iostream>
#include <list>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

list<Eigen::VectorXd> *k_means(Eigen::VectorXd *vectors, int vector_num, int group_num)
{
	int dimension = vectors[0].size();
	bool SwapFlag = true;
	int *group_id = new int[vector_num];
	Eigen::VectorXd *group_center = new Eigen::VectorXd[group_num];
	list<Eigen::VectorXd> *vectorGroup = new list<Eigen::VectorXd>[group_num];
	for(int i=0; i<group_num; i++) // initialize group center
	{
		group_center[i] = vectors[i];
	}
	double min_distance;
	int group_index;
	while(SwapFlag)
	{
		SwapFlag = false;
		// assign each vector to its nearest group
		for(int i=0; i<vector_num; i++)
		{
			min_distance = (vectors[i] - group_center[0]).norm();
			group_index = 0;
			for(int j=1; j<group_num; j++)
			{
				double distance = (vectors[i] - group_center[j]).norm();
				if(distance < min_distance)
				{
					min_distance = distance;
					group_index = j;
				}
			}
			if(group_index != group_id[i])
			{
				SwapFlag = true;
				group_id[i] = group_index;
			}
			vectorGroup[group_index].push_back(vectors[i]);
		}

		// update center value
		for(int i=0; i<group_num; i++)
		{
			Eigen::VectorXd mean;
			mean.setZero(dimension);
			list<Eigen::VectorXd>::iterator iter;
    		for(iter = vectorGroup[i].begin() ; iter != vectorGroup[i].end() ; ++iter)
				mean += *iter;
			mean /= vectorGroup[i].size();
			group_center[i] = mean;
			if(SwapFlag == true)
				vectorGroup[i].clear();
		}
	}
	delete group_id;
	//delete group_center; // this would cause error: pointer being freed was not allocated
	return vectorGroup;
}

void check_kmean_func(list<Eigen::VectorXd> *vectorGroup, int group_num)
{	
	int dimension = vectorGroup[0].begin()->size();
	Eigen::VectorXd *group_mean = new Eigen::VectorXd[group_num];
	for(int i=0; i< group_num; i++)
	{
		group_mean[i].setZero(dimension);
		cout << "-----group " << i << "------" << endl;
		list<Eigen::VectorXd>::iterator iter;
		for(iter = vectorGroup[i].begin() ; iter != vectorGroup[i].end() ; ++iter)
		{
			group_mean[i] += *iter;
			cout << "node:" << *iter << endl;
		}
		group_mean[i] /= vectorGroup[i].size();
		cout << "group mean = " << group_mean[i] << endl;
	}
	for(int i=0; i<group_num; i++)
	{
		list<Eigen::VectorXd>::iterator iter;
		for(iter = vectorGroup[i].begin() ; iter != vectorGroup[i].end() ; ++iter)
		{
			double min_distance = (*iter - group_mean[i]).norm();
			for(int j=0; j<group_num; j++)
			{
				assert(min_distance <= (*iter-group_mean[j]).norm());	
			}
		}
	}
	delete group_mean;
	return ;
}