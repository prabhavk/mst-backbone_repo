#ifndef MODELSELECTOR_H
#define MODELSELECTOR_H

#include <random>
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <boost/algorithm/string.hpp>
#include "SEM_MS.h"
#include "rootedPhylogeny_MS.h"
using namespace std;
class ModelSelector {
public:
	string sequenceFileName;
	string undirectedEdgeListFileName;
	string u_name;
	string v_name;
	string outputFilePrefix;
	std::chrono::system_clock::time_point start_time;
	std::chrono::system_clock::time_point current_time;
	ofstream logFile;
	vector <pair<int, int>> edgesForRooting;
	void PerformModelSelection();
	// void PerformModelSelectionForRootedTree();
	void PerformModelSelectionTest();
	void WriteTotalComputeTime();
	rootedPhylogeny_tree * RT;
	ModelSelector(string sequenceFileNameToSet, string undirectedEdgeListFileNameToSet, string u_nameToSet, string v_nameToSet) {
		this->sequenceFileName = sequenceFileNameToSet;
		this->undirectedEdgeListFileName = undirectedEdgeListFileNameToSet;
		this->u_name = u_nameToSet;
		this->v_name = v_nameToSet;
		outputFilePrefix = this->sequenceFileName + ".rootedAt_" + u_name + "_" + v_name;
		start_time = std::chrono::high_resolution_clock::now();		
		this->RT = new rootedPhylogeny_tree();
		this->logFile.open(outputFilePrefix + ".log");		

	}
	~ ModelSelector() {
		delete RT;
	}
};

// void ModelSelector::PerformModelSelection() {
// 	this->RT->ReadUndirectedEdgeList(this->undirectedEdgeListFileName);
// 	this->RT->SetSequenceFileName(this->sequenceFileName);
// 	// iterate over edges, root at each edge, and perform model selection for rooted tree
// }

void ModelSelector::PerformModelSelection() {
	int u_id; int v_id;	
	RT->ReadUndirectedEdgeList(this->undirectedEdgeListFileName);		
//	RT->ReadDirectedEdgeListForBifurcatingTree(this->directedEdgeListFileName);
	RT->SetSequenceFileName(this->sequenceFileName);
	
	u_id = RT->GetVertexId(this->u_name);
	v_id = RT->GetVertexId(this->v_name);
	pair <int, int> edge;
	if (u_id < v_id) {
		edge = pair<int,int>(u_id,v_id);
	} else {
		edge = pair<int,int>(v_id,u_id);
	}
	
	RT->RootTreeAtEdge(edge.first, edge.second);	
	RT->ReadSequenceFile(this->sequenceFileName);	
	RT->PerformModelSelectionUsingNelderMead();
	
	this->logFile << RT->optimizationLogString;
	cout << "BIC: "<< "\t" << setprecision(10) << RT->minimum_BIC_for_rooted_tree << endl;	
	this->logFile << "BIC: "<< "\t" << setprecision(10) << RT->minimum_BIC_for_rooted_tree << endl;	
//	RT->StoreGloballyOptimalEdgeLengths();
//	RT->StoreGloballyOptimalRateMatrices();
//	RT->StoreGloballyOptimalRateCategories();
	RT->WriteNewickFile(outputFilePrefix);
	RT->WriteEdgeList(outputFilePrefix);
	RT->WriteRateCategoryPerVertexAndModelParameters(outputFilePrefix);
}


void ModelSelector::PerformModelSelectionTest() {
	int u_id; int v_id;
	RT->ReadUndirectedEdgeList(this->undirectedEdgeListFileName);	
//	RT->ReadDirectedEdgeListForBifurcatingTree(this->directedEdgeListFileName);
	RT->SetSequenceFileName(this->sequenceFileName);
	u_id = RT->GetVertexId(this->u_name);
	v_id = RT->GetVertexId(this->v_name);
	pair <int, int> edge;
	if (u_id < v_id) {
		edge = pair<int,int>(u_id,v_id);
	} else {
		edge = pair<int,int>(v_id,u_id);
	}
	RT->RootTreeAtEdge(edge.first, edge.second);		
	RT->ReadSequenceFile(this->sequenceFileName);
	RT->PerformModelSelectionUsingNelderMead();
	cout << "BIC: "<< "\t" << RT->minimum_BIC_for_rooted_tree << endl;	
	this->logFile << "BIC: "<< "\t" << RT->minimum_BIC_for_rooted_tree << endl;	
//	RT->StoreGloballyOptimalEdgeLengths();
//	RT->StoreGloballyOptimalRateMatrices();
//	RT->StoreGloballyOptimalRateCategories();
	RT->WriteNewickFile(outputFilePrefix);
	RT->WriteEdgeList(outputFilePrefix);
	RT->WriteRateCategoryPerVertexAndModelParameters(outputFilePrefix);
}

void ModelSelector::WriteTotalComputeTime() {	
	this->current_time = std::chrono::high_resolution_clock::now();
	cout << "Total CPU time used is " << chrono::duration<double>(this->current_time - this->start_time).count() << " seconds\n";
	this->logFile << "Total CPU time used is " << chrono::duration<double>(this->current_time - this->start_time).count() << " seconds\n";
	this->logFile.close();
}

#endif