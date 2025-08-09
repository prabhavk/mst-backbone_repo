#ifndef mstBackbone_H
#define mstBackbone_H

#include <string>
#include "MST.h"
#include <tuple>
#include <iostream>
#include <stdio.h>
#include <pthread.h>
#include "SEM.h"
#include "utilities.h"
#include <cstring> 
// #include <boost/algorithm/string.hpp>
#include <chrono>
//#include "globalPhylogeny.h"
//#include "EM.h"
//#include "rootedPhylogeny.h"
//#include "optimizer.h"

// using namespace Eigen;
using namespace std;
// using namespace std::chrono;
class MSTBackbone
{
private:
	default_random_engine generator;
	vector <string> sequenceNames;
	map <string,unsigned char> mapDNAtoInteger;		
	ofstream mstBackboneLogFile;
	int numberOfLargeEdgesThreshold; // equal to num of hairs for the phc algorithm
	int numberOfHiddenVertices = 0;
	int edgeWeightThreshold;	
	chrono::system_clock::time_point start_time;
	chrono::system_clock::time_point current_time;
	chrono::system_clock::time_point t_start_time;
	chrono::system_clock::time_point t_end_time;
	chrono::system_clock::time_point m_start_time;
	chrono::system_clock::time_point m_end_time;
	// chrono::seconds timeTakenToComputeEdgeAndVertexLogLikelihoods;
	chrono::duration<double> timeTakenToComputeEdgeAndVertexLogLikelihoods;
	// chrono::seconds timeTakenToComputeGlobalUnrootedPhylogeneticTree;
	chrono::duration<double> timeTakenToComputeGlobalUnrootedPhylogeneticTree;
	// chrono::seconds timeTakenToComputeSubtree;
	chrono::duration<double> timeTakenToComputeSubtree;
	// chrono::seconds timeTakenToComputeSupertree;
	chrono::duration<double> timeTakenToComputeSupertree;
	// chrono::seconds timeTakenToRootViaEdgeLoglikelihoods;
	chrono::duration<double> timeTakenToRootViaEdgeLoglikelihoods;
	// chrono::seconds timeTakenToRootViaRestrictedSEM;
	chrono::duration<double> timeTakenToRootViaRestrictedSEM;
	string sequenceFileName;	
	string prefix_for_output_files;
	string ancestralSequencesString;
	string MSTFileName;
	string GMMparametersFileName;
	string distance_measure_for_NJ = "Hamming";
	bool apply_patch = false;
	bool grow_tree_incrementally = false;
	// Remove this
	int ComputeHammingDistance(string seq1, string seq2);
//	float ComputeVertexOrderPerturbedDistance(vector<unsigned char> recodedSeq1, vector<unsigned char> recodedSeq2);	
	int ComputeHammingDistance(vector<unsigned char> recodedSeq1, vector<unsigned char> recodedSeq2);
	int GetEdgeIndex (int vertexIndex1, int vertexIndex2, int numberOfVertices);
	MST_tree * M;
//	globalPhylogeny * T;
	SEM * T;
	SEM * t;
//	phylogeny_tree * P_ptr;
//	rootedPhylogeny_tree * RT_ptr;
//	void ComputeMST(string sequenceFileName);
//	void ComputeVMST(string sequenceFileName);
	void WriteOutputFiles();
	bool debug;
	bool verbose;
	bool localPhyloOnly;	
	bool modelSelection;
	string modelForRooting = "UNREST";
	string supertree_method;
	int numberOfVerticesInSubtree;
	string GetSequenceListToWriteToFile(map <string, vector <unsigned char>> compressedSeqMap, vector <vector <int> > sitePatternRepetitions);
	vector <string> must_have;
	vector <string> may_have;
public:
	void SetDNAMap();
	void SetThresholds();
	void MSTBackboneWithOneExternalVertex();
	void MSTBackboneWithFullSEMAndMultipleExternalVertices();	
	void RootSuperTree();
	void MSTBackboneWithRootSEMAndMultipleExternalVertices();
	void MSTBackboneOverlappingSets();
	void MSTBackboneOnlyLocalPhylo();	
	MSTBackbone(string sequenceFileNameToAdd, int subtreeSizeThresholdToset, string prefix_for_output_files_to_set, string distance_measure_for_NJ_to_set, bool verbose_flag_to_set, bool flag_root_supertree, string supertree_method_to_set) {
		// MSTBackbone(string sequenceFileNameToAdd, int subtreeSizeThresholdToset, string prefix_for_output_files_to_set, bool localPhyloOnly_to_set, bool modelSelection_to_set, string modelForRooting_to_set, bool useChowLiu_toset) {
		// bool localPhyloOnly = TRUE;
		// this->useChowLiu = useChowLiu_toset;
		// this->localPhyloOnly = localPhyloOnly_to_set;		
		// this->modelForRooting = modelForRooting_to_set;		
		start_time = chrono::high_resolution_clock::now();				
		this->prefix_for_output_files = prefix_for_output_files_to_set;
		// output files		
		this->mstBackboneLogFile.open(this->prefix_for_output_files + ".mstbackbone_log");
		this->sequenceFileName = sequenceFileNameToAdd;		
		this->supertree_method = supertree_method_to_set;		
		this->verbose = verbose_flag_to_set;
		this->distance_measure_for_NJ = distance_measure_for_NJ_to_set;
		cout << "Distance measure used for NJ is " << this->distance_measure_for_NJ << endl;
		this->mstBackboneLogFile << "Distance measure used for NJ is " << this->distance_measure_for_NJ << endl;
		if (flag_root_supertree) {
			cout << "Supertree will be rooted" << endl;
			this->mstBackboneLogFile << "Supertree will be rooted " << endl;
		} else {
			cout << "Supertree will not be rooted" << endl;
			this->mstBackboneLogFile << "Supertree will not be rooted " << endl;
		}				
		this->numberOfLargeEdgesThreshold = subtreeSizeThresholdToset;		
		mstBackboneLogFile << "Constraint size is set at\t" << this->numberOfLargeEdgesThreshold << endl;
		cout << "Constraint size is set at\t" << this->numberOfLargeEdgesThreshold << endl;
		mstBackboneLogFile << "Prefix for output files is \t" << this->prefix_for_output_files << endl;
		cout << "Prefix for output files is \t" << this->prefix_for_output_files << endl;
		MSTFileName = prefix_for_output_files + ".initial_MST";
		this->SetDNAMap();
		this->ancestralSequencesString = "";								
		this->m_start_time = std::chrono::high_resolution_clock::now();
		this->M = new MST_tree();
		this->M->ReadSequences(this->sequenceFileName);
		this->M->ComputeMST();
		cout << this->M->num_duplicated_sequences << " duplicate sequences found; duplicate sequences will be not be used by mst-backbone; instead they will be added to the supertree constructed by mst-backbone" << endl;
		this->mstBackboneLogFile << this->M->num_duplicated_sequences << " duplicate sequences found; duplicate sequences will be not be used by mst-backbone; instead they will be added to the supertree constructed by mst-backbone" << endl;
		if (true) { // replace with user input
			this->M->WriteToFile(MSTFileName);
		}		
		this->current_time = std::chrono::high_resolution_clock::now();
		cout << "Time taken to compute MST is " << chrono::duration<double>(this->current_time-this->m_start_time).count() << " seconds\n";
		this->mstBackboneLogFile << "Time taken to compute MST is " << chrono::duration<double>(this->current_time-this->m_start_time).count() << " seconds\n";
	    // Compute Chow-Liu tree using UNREST and get probability distribution for root position
		this->M->SetNumberOfLargeEdgesThreshold(this->numberOfLargeEdgesThreshold);
		this->T = new SEM(1,this->distance_measure_for_NJ,this->verbose);
		this->T->SetStream(this->mstBackboneLogFile);
		int numberOfInputSequences = (int) this->M->vertexMap->size();
		this->T->numberOfObservedVertices = numberOfInputSequences;
		this->m_start_time = std::chrono::high_resolution_clock::now();
		// timeTakenToComputeGlobalUnrootedPhylogeneticTree -= timeTakenToComputeEdgeAndVertexLogLikelihoods;		
		if (this->supertree_method == "mstbackbone") { // MST_BACKBONE BioRxiv v1 and v2 [rejected], Bioinformatics submission
			cout << "Supertree method is " << this->supertree_method << endl;
			this->MSTBackboneWithFullSEMAndMultipleExternalVertices();
		}
		
		this->current_time = std::chrono::high_resolution_clock::now();
		cout << "Time taken for computing unrooted supertree is " << chrono::duration<double>(this->current_time-this->m_start_time).count() << " seconds\n";
		this->mstBackboneLogFile << "Time taken for computing unrooted supertree is " << chrono::duration<double>(this->current_time-this->m_start_time).count() << " seconds\n";
		if (flag_root_supertree){
			cout << "Rooting supertree by fitting the GMM via EM" << endl;
			this->RootSuperTree();
		}
		this->current_time = std::chrono::high_resolution_clock::now();
		cout << "Total CPU time used is " << chrono::duration<double>(this->current_time-this->start_time).count() << " seconds\n";
		this->mstBackboneLogFile << "Total CPU time used is " << chrono::duration<double>(this->current_time-this->start_time).count() << " seconds\n";
		this->mstBackboneLogFile.close();
			}
	~MSTBackbone(){
		delete this->T;
		delete this->M;	
	}
};

void MSTBackbone::SetDNAMap() {
	this->mapDNAtoInteger["A"] = 0;
	this->mapDNAtoInteger["C"] = 1;
	this->mapDNAtoInteger["G"] = 2;
	this->mapDNAtoInteger["T"] = 3;
}

void MSTBackbone::MSTBackboneOnlyLocalPhylo() {
	vector <string> names;
	vector <vector <unsigned char> > sequences;
	vector <int> sitePatternWeights;
	vector <vector <int> > sitePatternRepetitions;	
	vector <int> idsOfVerticesToRemove;
	vector <int> idsOfVerticesToKeep;
	vector <int> idsOfExternalVertices;
	vector <int> idsOfVerticesForSEM;
	vector <tuple <int, string, vector <unsigned char>>> idAndNameAndSeqTupleForVerticesToAdd;	
	//----##############################################################---//
	//	1.	Initialize the global phylogenetic tree T as the empty graph   //
	//----##############################################################---//
//	cout << "Starting MST-backbone" << endl;
//	cout << "1.	Initialize the global phylogenetic tree T as the empty graph" << endl;
	int numberOfInputSequences = (int) this->M->vertexMap->size();		
	current_time = chrono::high_resolution_clock::now();
	timeTakenToComputeEdgeAndVertexLogLikelihoods = chrono::duration<double>(current_time-current_time);
	
	// Initialize global phylogeny
	// idsOfVerticesForSEM.clear();
	// for (pair <int, MST_vertex *> vIdAndPtr : * this->M->vertexMap) {
	// 	idsOfVerticesForSEM.push_back(vIdAndPtr.first);
	// }
	// tie (names, sequences, sitePatternWeights, sitePatternRepetitions) = this->M->GetCompressedSequencesSiteWeightsAndSiteRepeats(idsOfVerticesForSEM);	
	// this->T->sequenceFileName = this->sequenceFileName;
	// this->T->AddSequences(sequences);
	// this->T->AddNames(names);
	// this->T->AddSitePatternWeights(sitePatternWeights);
	// this->T->SetNumberOfInputSequences(numberOfInputSequences);	
	// this->T->numberOfObservedVertices = numberOfInputSequences;
	
	int largestIdOfVertexInMST = numberOfInputSequences;
	
	bool computeLocalPhylogeneticTree = 1;
	bool numberOfNonSingletonComponentsIsGreaterThanZero = 0;	
	
	while (computeLocalPhylogeneticTree) {
		cout << "Number of vertices in MST is " << this->M->vertexMap->size() << endl;
		this->mstBackboneLogFile << "Number of vertices in MST is " << this->M->vertexMap->size() << endl;
//		cout << "Max vertex degree in MST is " << this->M->maxDegree << endl;
		//----####################################################################---//
		//	2.	Compute the size of the smallest subtree ts = (Vs,Es) of M s.t.		 //
		//		|Vs| > s. Check if |Vm\Vs| > s.								   		 //
		// 		If yes then go to step 3 else go to step 9					   		 //
		// 		Bootstrapped alignments may contain zero-weight edges		   		 //
		//      If so then replace |Vs| with |{non-zero weighted edges in Es}| 		 //
		//      Additionally replace |Vm\Vs| with |{non-zero weighted edges in Es}|  //
		//----####################################################################---//				
		computeLocalPhylogeneticTree = this->M->ShouldIComputeALocalPhylogeneticTree();
//		cout << "2. Checking if local phylogenetic tree should be computed" << endl;
		if (computeLocalPhylogeneticTree) {
			//----####################################################################---//
			//	3.	Extract vertices inducing subtree (Vs), and external vertices (Ve)	 //
			//----####################################################################---//				
//			cout << "3. Extract vertices inducing subtree (Vs), and external vertices (Ve)" << endl;
			this->M->SetIdsOfExternalVertices();
			idsOfExternalVertices = this->M->idsOfExternalVertices;
			idsOfVerticesForSEM = this->M->subtree_v_ptr->idsOfVerticesInSubtree;
			this->numberOfVerticesInSubtree = this->M->subtree_v_ptr->idsOfVerticesInSubtree.size();
			for (int id: idsOfExternalVertices) {
				idsOfVerticesForSEM.push_back(id);
			}
			tie (names, sequences, sitePatternWeights, sitePatternRepetitions) = this->M->GetCompressedSequencesSiteWeightsAndSiteRepeats(idsOfVerticesForSEM);
			//----########################################################---//
			//	4.	Compute local phylogeny t over (Vs U Ve) via SEM      	 //
			//----########################################################---//
//			cout << "4.	Compute local phylogeny t over (Vs U Ve) via SEM" << endl;			
			this->t = new SEM(largestIdOfVertexInMST,this->distance_measure_for_NJ,this->verbose);
			this->t->AddSequences(sequences);
			this->t->SetNumberOfVerticesInSubtree(this->numberOfVerticesInSubtree);
			this->t->SetNumberOfInputSequences(numberOfInputSequences);
			this->t->AddRootVertex();
			this->t->AddNames(names);
			this->t->AddGlobalIds(idsOfVerticesForSEM);
			this->t->AddSitePatternWeights(sitePatternWeights);
			this->t->AddSitePatternRepeats(sitePatternRepetitions);			
			this->t->OptimizeTopologyAndParametersOfGMM();			
//			timeTakenToComputeUnrootedPhylogeny += chrono::duration_cast<chrono::seconds>(t_end_time - t_start_time);
			//----##################################################################---//	
			//  5.	Check if # of non-singleton components of forest f in t that       //
			//		is induced by Vs is greater than zero.							   //
			//		i.e., Does local phylogeny contain vertices/edges of interest?	   //
			//----##################################################################---//
			this->t->SelectIndsOfVerticesOfInterestAndEdgesOfInterest();
			numberOfNonSingletonComponentsIsGreaterThanZero = this->t->IsNumberOfNonSingletonComponentsGreaterThanZero();			
//			cout << "5. Checking if there are any vertices of interest" << endl;
			if (!numberOfNonSingletonComponentsIsGreaterThanZero) {
				//----####################################################---//	
				//  6.	If no then double subtree size and go to step 2 	 //
				//		else reset subtree size and go to to step 7		     //
				//----####################################################---//		
				this->M->DoubleSubtreeSizeThreshold();
//				cout << "6. Doubling subtree size" << endl;
			} else {
				this->M->ResetSubtreeSizeThreshold();		
				//----################################---//
				//  7.	Add vertices/edges in f to T     //
				//----################################---//
//				cout << "7. Adding vertices/edges in f to T" << endl;
				this->t->RenameHiddenVerticesInEdgesOfInterestAndSetIdsOfVerticesOfInterest();
				this->t->SetWeightedEdgesToAddToGlobalPhylogeneticTree();				
				//this->T->AddWeightedEdges(this->t->weightedEdgesToAddToGlobalPhylogeneticTree);
				this->t->SetAncestralSequencesString();
				this->ancestralSequencesString += this->t->ancestralSequencesString;				
				t_start_time = chrono::high_resolution_clock::now();
				this->t->SetEdgeAndVertexLogLikelihoods();				
				// this->T->AddVertexLogLikelihoods(this->t->vertexLogLikelihoodsMapToAddToGlobalPhylogeneticTree);
				// this->T->AddEdgeLogLikelihoods(this->t->edgeLogLikelihoodsToAddToGlobalPhylogeneticTree);
				t_end_time = chrono::high_resolution_clock::now();
				timeTakenToComputeEdgeAndVertexLogLikelihoods += t_end_time - t_start_time;
				// Add vertex logLikelihoods
				// Add edge logLikelihoods
				largestIdOfVertexInMST = this->t->largestIdOfVertexInMST;
				//----##############################---//
				//  8.	Update M and go to step 1	   //
				//----##############################---//
//				cout << "8. Updating MST" << endl;
				this->t->SetInfoForVerticesToAddToMST();				
				this->M->UpdateMSTWithMultipleExternalVertices(t->idsOfVerticesToKeepInMST, t->idsOfVerticesToRemove, t->idAndNameAndSeqTuple, idsOfExternalVertices);								
				delete this->t;
			}			
			computeLocalPhylogeneticTree = this->M->ShouldIComputeALocalPhylogeneticTree();
		}		
		cout << "CPU time used for computing local phylogeny is " << chrono::duration<double>(t_end_time-t_start_time).count() << " second(s)\n";
		this->mstBackboneLogFile << "CPU time used for computing local phylogeny is " << chrono::duration<double>(t_end_time-t_start_time).count() << " second(s)\n";			
	}	
	//----########################################################---//
	//	9.	Compute phylogenetic tree t over vertices in M, and      //
	//		add vertices/edges in t to T							 //
	//----########################################################---//
	cout << "Computing phylogenetic tree over all vertices in MST" << endl;
	this->M->UpdateMaxDegree();
	cout << "Max vertex degree in MST is " << this->M->maxDegree << endl;
	idsOfVerticesForSEM.clear();
	for (pair <int, MST_vertex *> idPtrPair: * this->M->vertexMap) {
		idsOfVerticesForSEM.push_back(idPtrPair.first);
	}
	cout << "Number of vertices in MST is " << idsOfVerticesForSEM.size() << endl;
	cout << "Number of edges in MST is " << this->M->edgeWeightsMap.size() << endl;
	tie (names, sequences, sitePatternWeights, sitePatternRepetitions) = this->M->GetCompressedSequencesSiteWeightsAndSiteRepeats(idsOfVerticesForSEM);
	this->numberOfVerticesInSubtree = sequences.size();
	this->t = new SEM(largestIdOfVertexInMST,this->distance_measure_for_NJ,this->verbose);
	this->t->SetFlagForFinalIterationOfSEM();
	this->t->AddSequences(sequences);
	this->t->SetNumberOfVerticesInSubtree(this->numberOfVerticesInSubtree);
	this->t->SetNumberOfInputSequences(numberOfInputSequences);
	this->t->AddRootVertex();
	this->t->AddNames(names);
	this->t->AddGlobalIds(idsOfVerticesForSEM);
	this->t->AddSitePatternWeights(sitePatternWeights);
	this->t->AddSitePatternRepeats(sitePatternRepetitions);	
	this->t->OptimizeTopologyAndParametersOfGMM();			
	// timeTakenToComputeUnrootedPhylogeny += chrono::duration_cast<chrono::seconds>(t_end_time - t_start_time);
	this->t->SelectIndsOfVerticesOfInterestAndEdgesOfInterest();
	this->t->RenameHiddenVerticesInEdgesOfInterestAndSetIdsOfVerticesOfInterest();
	this->t->SetWeightedEdgesToAddToGlobalPhylogeneticTree();
	this->t->SetAncestralSequencesString();
	this->ancestralSequencesString += this->t->ancestralSequencesString;
	this->t->WriteAncestralSequences();	
	//this->T->AddWeightedEdges(this->t->weightedEdgesToAddToGlobalPhylogeneticTree);	
	t_start_time = chrono::high_resolution_clock::now();
	this->t->SetEdgeAndVertexLogLikelihoods();
	//this->T->AddVertexLogLikelihoods(this->t->vertexLogLikelihoodsMapToAddToGlobalPhylogeneticTree);
	//this->T->AddEdgeLogLikelihoods(this->t->edgeLogLikelihoodsToAddToGlobalPhylogeneticTree);
	t_end_time = chrono::high_resolution_clock::now();
	timeTakenToComputeEdgeAndVertexLogLikelihoods += t_end_time - t_start_time;
	delete this->t;	
	this->mstBackboneLogFile << "CPU time used for computing local phylogeny is " << chrono::duration<double>(t_end_time-t_start_time).count() << " second(s)\n";		
	// assert that T is a tree

}


void MSTBackbone::MSTBackboneWithFullSEMAndMultipleExternalVertices() {
	vector <string> names;
	vector <vector <unsigned char> > sequences;
	vector <int> sitePatternWeights;
	vector <vector <int> > sitePatternRepetitions;	
	vector <int> idsOfVerticesToRemove;
	vector <int> idsOfVerticesToKeep;
	vector <int> idsOfExternalVertices;
	vector <int> idsOfVerticesForSEM;
	vector <tuple <int, string, vector <unsigned char>>> idAndNameAndSeqTupleForVerticesToAdd;	
	//----##############################################################---//
	//	1.	Initialize the global phylogenetic tree T as the empty graph   //
	//----##############################################################---//
//	cout << "Starting MST-backbone" << endl;
//	cout << "1.	Initialize the global phylogenetic tree T as the empty graph" << endl;
	int numberOfInputSequences = (int) this->M->vertexMap->size();		
	current_time = chrono::high_resolution_clock::now();
	// timeTakenToComputeEdgeAndVertexLogLikelihoods = chrono::duration_cast<chrono::seconds>(current_time-current_time);
	
	// Initialize supertree
	idsOfVerticesForSEM.clear();
	for (pair <int, MST_vertex *> vIdAndPtr : * this->M->vertexMap) {
		idsOfVerticesForSEM.push_back(vIdAndPtr.first);
	}
	tie (names, sequences, sitePatternWeights, sitePatternRepetitions) = this->M->GetCompressedSequencesSiteWeightsAndSiteRepeats(idsOfVerticesForSEM);	
	this->T->sequenceFileName = this->sequenceFileName;
	this->T->AddSequences(sequences);
	this->T->AddNames(names);
	this->T->AddSitePatternWeights(sitePatternWeights);
	this->T->SetNumberOfInputSequences(numberOfInputSequences);	
	this->T->numberOfObservedVertices = numberOfInputSequences;
	// add duplicated sequences here
	
	int largestIdOfVertexInMST = numberOfInputSequences;
	
	bool computeLocalPhylogeneticTree = 1;
	bool numberOfNonSingletonComponentsIsGreaterThanZero = 0;	
	
	while (computeLocalPhylogeneticTree) {
		// current_time = chrono::high_resolution_clock::now();
		cout << "Number of vertices in MST is " << this->M->vertexMap->size() << endl;				
		this->mstBackboneLogFile << "Number of vertices in MST is " << this->M->vertexMap->size() << endl;
//		cout << "Max vertex degree in MST is " << this->M->maxDegree << endl;
		//----####################################################################---//
		//	2.	Compute the size of the smallest subtree ts = (Vs,Es) of M s.t.		 //
		//		|Vs| > s. Check if |Vm\Vs| > s.								   		 //
		// 		If yes then go to step 3 else go to step 9					   		 //
		// 		Bootstrapped alignments may contain zero-weight edges		   		 //
		//      If so then replace |Vs| with |{non-zero weighted edges in Es}| 		 //
		//      Additionally replace |Vm\Vs| with |{non-zero weighted edges in Es}|  //
		//----####################################################################---//				
		computeLocalPhylogeneticTree = this->M->ShouldIComputeALocalPhylogeneticTree();
//		cout << "2. Checking if local phylogenetic tree should be computed" << endl;
		if (computeLocalPhylogeneticTree) {
			//----####################################################################---//
			//	3.	Extract vertices inducing subtree (Vs), and external vertices (Ve)	 //
			//----####################################################################---//				
//			cout << "3. Extract vertices inducing subtree (Vs), and external vertices (Ve)" << endl;
			this->M->SetIdsOfExternalVertices();
			idsOfExternalVertices = this->M->idsOfExternalVertices;
			idsOfVerticesForSEM = this->M->subtree_v_ptr->idsOfVerticesInSubtree;
			this->numberOfVerticesInSubtree = this->M->subtree_v_ptr->idsOfVerticesInSubtree.size();
			for (int id: idsOfExternalVertices) {
				idsOfVerticesForSEM.push_back(id);
			}
			tie (names, sequences, sitePatternWeights, sitePatternRepetitions) = this->M->GetCompressedSequencesSiteWeightsAndSiteRepeats(idsOfVerticesForSEM);
			//----########################################################---//
			//	4.	Compute local phylogeny t over (Vs U Ve) via SEM      	 //
			//----########################################################---//
//			cout << "4.	Compute local phylogeny t over (Vs U Ve) via SEM" << endl;			
			this->t = new SEM(largestIdOfVertexInMST,this->distance_measure_for_NJ, this->verbose);
			this->t->AddSequences(sequences);
			this->t->SetNumberOfVerticesInSubtree(this->numberOfVerticesInSubtree);
			this->t->SetNumberOfInputSequences(numberOfInputSequences);
			this->t->AddRootVertex();
			this->t->AddNames(names);
			this->t->AddGlobalIds(idsOfVerticesForSEM);
			this->t->AddSitePatternWeights(sitePatternWeights);
			this->t->AddSitePatternRepeats(sitePatternRepetitions);
			t_start_time = chrono::high_resolution_clock::now();
			this->t->OptimizeTopologyAndParametersOfGMM();
			t_end_time = chrono::high_resolution_clock::now();
			// timeTakenToComputeSubtree = chrono::duration_cast<chrono::seconds>(t_end_time - t_start_time);
			timeTakenToComputeSubtree = t_end_time - t_start_time;
			cout << "CPU time used for computing subtree with " << this->t->numberOfObservedVertices << " leaves is " << timeTakenToComputeSubtree.count() << " seconds\n";
			this->mstBackboneLogFile << "CPU time used for computing subtree with " << this->t->numberOfObservedVertices << " leaves is " << timeTakenToComputeSubtree.count() << " seconds\n";
//			timeTakenToComputeUnrootedPhylogeny += chrono::duration_cast<chrono::seconds>(t_end_time - t_start_time);
			//----##################################################################---//	
			//  5.	Check if # of non-singleton components of forest f in t that       //
			//		is induced by Vs is greater than zero.							   //
			//		i.e., Does local phylogeny contain vertices/edges of interest?	   //
			//----##################################################################---//
			this->t->SelectIndsOfVerticesOfInterestAndEdgesOfInterest();
			numberOfNonSingletonComponentsIsGreaterThanZero = this->t->IsNumberOfNonSingletonComponentsGreaterThanZero();			
//			cout << "5. Checking if there are any vertices of interest" << endl;
			if (!numberOfNonSingletonComponentsIsGreaterThanZero) {
				//----####################################################---//	
				//  6.	If no then double subtree size and go to step 2 	 //
				//		else reset subtree size and go to to step 7		     //
				//----####################################################---//		
				this->M->DoubleSubtreeSizeThreshold();
//				cout << "6. Doubling subtree size" << endl;
			} else {
				this->M->ResetSubtreeSizeThreshold();		
				//----################################---//
				//  7.	Add vertices/edges in f to T     //
				//----################################---//
//				cout << "7. Adding vertices/edges in f to T" << endl;
				this->t->RenameHiddenVerticesInEdgesOfInterestAndSetIdsOfVerticesOfInterest();
				this->t->SetWeightedEdgesToAddToGlobalPhylogeneticTree();				
				this->T->AddWeightedEdges(this->t->weightedEdgesToAddToGlobalPhylogeneticTree);
				this->t->SetAncestralSequencesString();
				this->ancestralSequencesString += this->t->ancestralSequencesString;				
				// t_start_time = chrono::high_resolution_clock::now();
				// this->t->SetEdgeAndVertexLogLikelihoods();				
				// this->T->AddVertexLogLikelihoods(this->t->vertexLogLikelihoodsMapToAddToGlobalPhylogeneticTree);
				// this->T->AddEdgeLogLikelihoods(this->t->edgeLogLikelihoodsToAddToGlobalPhylogeneticTree);
				// // this->T->AddExpectedCountMatrices(t->expectedCountsForVertexPair);
				// t_end_time = chrono::high_resolution_clock::now();
				// timeTakenToComputeEdgeAndVertexLogLikelihoods += chrono::duration_cast<chrono::seconds>(t_end_time - t_start_time);
				// Add vertex logLikelihoods
				// Add edge logLikelihoods
				largestIdOfVertexInMST = this->t->largestIdOfVertexInMST;
				//----##############################---//
				//  8.	Update M and go to step 1	   //
				//----##############################---//
//				cout << "8. Updating MST" << endl;
				this->t->SetInfoForVerticesToAddToMST();				
				this->M->UpdateMSTWithMultipleExternalVertices(t->idsOfVerticesToKeepInMST, t->idsOfVerticesToRemove, t->idAndNameAndSeqTuple, idsOfExternalVertices);								
				delete this->t;
			}			
			computeLocalPhylogeneticTree = this->M->ShouldIComputeALocalPhylogeneticTree();
		}		
//		cout << "CPU time used for computing local phylogeny is " << chrono::duration_cast<chrono::seconds>(t_end_time-t_start_time).count() << " second(s)\n";
//		this->mstBackboneLogFile << "CPU time used for computing local phylogeny is " << chrono::duration_cast<chrono::seconds>(t_end_time-t_start_time).count() << " second(s)\n";			
	}	
	//----########################################################---//
	//	9.	Compute phylogenetic tree t over vertices in M, and      //
	//		add vertices/edges in t to T							 //
	//----########################################################---//
	cout << "Computing phylogenetic tree over all vertices in MST" << endl;
//	this->M->UpdateMaxDegree();
//	cout << "Max vertex degree in MST is " << this->M->maxDegree << endl;
	idsOfVerticesForSEM.clear();
	for (pair <int, MST_vertex *> idPtrPair: * this->M->vertexMap) {
		idsOfVerticesForSEM.push_back(idPtrPair.first);
	}
//	cout << "Number of vertices in MST is " << idsOfVerticesForSEM.size() << endl;
//	cout << "Number of edges in MST is " << this->M->edgeWeightsMap.size() << endl;
	tie (names, sequences, sitePatternWeights, sitePatternRepetitions) = this->M->GetCompressedSequencesSiteWeightsAndSiteRepeats(idsOfVerticesForSEM);
	this->numberOfVerticesInSubtree = sequences.size();
	this->t = new SEM(largestIdOfVertexInMST,this->distance_measure_for_NJ,this->verbose);
	this->t->SetFlagForFinalIterationOfSEM();
	this->t->AddSequences(sequences);
	this->t->SetNumberOfVerticesInSubtree(this->numberOfVerticesInSubtree);
	this->t->SetNumberOfInputSequences(numberOfInputSequences);
	this->t->AddRootVertex();
	this->t->AddNames(names);
	this->t->AddGlobalIds(idsOfVerticesForSEM);
	this->t->AddSitePatternWeights(sitePatternWeights);
	this->t->AddSitePatternRepeats(sitePatternRepetitions);	
	t_start_time = chrono::high_resolution_clock::now();
	this->t->OptimizeTopologyAndParametersOfGMM();
	t_end_time = chrono::high_resolution_clock::now();
	// timeTakenToComputeSubtree = chrono::duration_cast<chrono::seconds>(t_end_time - t_start_time);
	timeTakenToComputeSubtree = t_end_time - t_start_time;
	cout << "CPU time used for computing subtree with " << this->t->numberOfObservedVertices << " leaves is " << timeTakenToComputeSubtree.count() << " seconds\n";
	this->mstBackboneLogFile << "CPU time used for computing subtree with " << this->t->numberOfObservedVertices << " leaves is " << timeTakenToComputeSubtree.count() << " seconds\n";
	// timeTakenToComputeUnrootedPhylogeny += chrono::duration_cast<chrono::seconds>(t_end_time - t_start_time);
	this->t->SelectIndsOfVerticesOfInterestAndEdgesOfInterest();
	this->t->RenameHiddenVerticesInEdgesOfInterestAndSetIdsOfVerticesOfInterest();
	this->t->SetWeightedEdgesToAddToGlobalPhylogeneticTree();
	this->t->SetAncestralSequencesString();
	this->ancestralSequencesString += this->t->ancestralSequencesString;
	this->t->WriteAncestralSequences();	
	this->T->AddWeightedEdges(this->t->weightedEdgesToAddToGlobalPhylogeneticTree);	
	// t_start_time = chrono::high_resolution_clock::now();
	// this->t->SetEdgeAndVertexLogLikelihoods();
	// this->T->AddVertexLogLikelihoods(this->t->vertexLogLikelihoodsMapToAddToGlobalPhylogeneticTree);
	// this->T->AddEdgeLogLikelihoods(this->t->edgeLogLikelihoodsToAddToGlobalPhylogeneticTree);
	// // this->T->AddExpectedCountMatrices(this->t->expectedCountsForVertexPair);
	// t_end_time = chrono::high_resolution_clock::now();
	// timeTakenToComputeEdgeAndVertexLogLikelihoods += chrono::duration_cast<chrono::seconds>(t_end_time - t_start_time);
	delete this->t;
	//	this->mstBackboneLogFile << "CPU time used for computing local phylogeny is " << chrono::duration_cast<chrono::seconds>(t_end_time-t_start_time).count() << " second(s)\n";		
	// assert that T is a tree
	// cout << "Number of vertices in T is " << this->T->vertexMap->size() << endl;
	// cout << "Number of edges in T is " << this->T->edgeLengths.size() << endl;
	// assert(this->T->vertexMap->size() == this->T->edgeLengths.size() + 1);
	// timeTakenToComputeGlobalUnrootedPhylogeneticTree = chrono::duration_cast<chrono::seconds>(current_time-start_time);	
	// timeTakenToComputeGlobalUnrootedPhylogeneticTree -= timeTakenToComputeEdgeAndVertexLogLikelihoods;		
	cout << "Adding duplicated sequences to tree" << endl;
	this->mstBackboneLogFile << "Adding duplicated sequences to tree" << endl;
	this->T->AddDuplicatedSequencesToUnrootedTree(this->M);
	this->T->WriteUnrootedTreeAsEdgeList(this->prefix_for_output_files + ".unrooted_edgeList");
	this->T->RootTreeAtAVertexPickedAtRandom();
	this->T->WriteRootedTreeInNewickFormat(this->prefix_for_output_files + ".unrooted_newick");	
}

void MSTBackbone::RootSuperTree() {
//----##############---//		
	//	10.	Root T via EM  //
	//----##############---//
	// cout << "Fitting a general Markov model GMM to T using reconstructed ancestral sequences" << endl;
	// this->mstBackboneLogFile << "Fitting a general Markov model GMM to T using reconstructed ancestral sequences" << endl;		
	// this->T->RootTreeBySumOfExpectedLogLikelihoods();
	// current_time = chrono::high_resolution_clock::now();
	// timeTakenToRootViaEdgeLoglikelihoods = chrono::duration_cast<chrono::seconds>(current_time-start_time);
	// timeTakenToRootViaEdgeLoglikelihoods -= timeTakenToComputeGlobalUnrootedPhylogeneticTree;
	// cout << "CPU time used for fitting a GM model to fully labeled T is " << timeTakenToRootViaEdgeLoglikelihoods.count() << " second(s)\n";
	// this->mstBackboneLogFile << "CPU time used for fitting a GM model to fully labeled T is " << timeTakenToRootViaEdgeLoglikelihoods.count() << " second(s)\n";
	// cout << "Log likelihood of fitting a GM model to fully labeled T is " << this->T->maxSumOfExpectedLogLikelihoods << endl;
	// this->mstBackboneLogFile << "Log likelihood of fitting a GM model to fully labeled T is " << this->T->maxSumOfExpectedLogLikelihoods << endl;
	// double BIC_full_labeled = -2 * this->T->maxSumOfExpectedLogLikelihoods ;
	// BIC_full_labeled += (3 + 12 * (this->T->numberOfInputSequences -1) * log2(this->T->sequenceLength));
	// cout << "BIC of fitting a GM model to fully labeled T is " << BIC_full_labeled << endl;
	// this->mstBackboneLogFile << "BIC of fitting a GM model to fully labeled T is " << BIC_full_labeled << endl;
	// cout << "Writing rooted tree in edge list format and newick format" << endl;
	// this->T->WriteRootedTreeAsEdgeList(sequenceFileName + ".edgeList_fullyLabeledRooting");
	// this->T->WriteRootedTreeInNewickFormat(sequenceFileName + ".newick_fullyLabeledRooting");

	cout << "Root T by fitting a GMM using EM" << endl;
	this->mstBackboneLogFile << "Root T by fitting a a GMM using EM" << endl;
	t_start_time = chrono::high_resolution_clock::now();
	this->T->RootTreeByFittingAGMMViaEM();
	t_end_time = chrono::high_resolution_clock::now();
	// current_time = chrono::high_resolution_clock::now();
	timeTakenToRootViaRestrictedSEM = t_end_time-t_start_time;
	// timeTakenToRootViaRestrictedSEM = chrono::duration_cast<chrono::seconds>(t_end_time-t_start_time);
	// timeTakenToRootViaRestrictedSEM -= timeTakenToComputeGlobalUnrootedPhylogeneticTree;
	// timeTakenToRootViaRestrictedSEM -= timeTakenToRootViaEdgeLoglikelihoods;
	cout << "CPU time used for rooting T using EM is " << timeTakenToRootViaRestrictedSEM.count() << " second(s)\n";
	this->mstBackboneLogFile << "CPU time used for rooting T using EM is " << timeTakenToRootViaRestrictedSEM.count() << " second(s)\n";
	cout << "Log likelihood under GMM is " << this->T->logLikelihood << endl;
	this->mstBackboneLogFile << "Log likelihood is " << this->T->logLikelihood << endl;
	double BIC = -2 * this->T->logLikelihood + (3 + 12 * (this->T->vertexMap->size() -1)) * log(this->T->sequenceLength);	
	cout << "BIC under GMM is " << BIC << endl;
	this->mstBackboneLogFile << "BIC under GMM is " << BIC << endl;
	cout << "Writing rooted tree in edge list format and newick format" << endl;
	this->T->WriteRootedTreeAsEdgeList(this->prefix_for_output_files + "_unique_seqs_only.directed_edges");
	this->T->WriteRootedTreeInNewickFormat(this->prefix_for_output_files + "_unique_seqs_only.rooted_newick");	
	cout << "Adding duplicated sequences " << endl;
	this->T->AddDuplicatedSequencesToRootedTree(this->M);
		
	// add Identity Matrices for duplicated sequences.
	this->T->WriteRootedTreeAsEdgeList(this->prefix_for_output_files + ".directed_edges");
	this->T->WriteRootedTreeInNewickFormat(this->prefix_for_output_files + ".rooted_newick");
	this->T->WriteParametersOfGMM(this->prefix_for_output_files + ".GMM_parameters"); // only valid for unique sequences, 	
}



void MSTBackbone::MSTBackboneWithRootSEMAndMultipleExternalVertices() {
	vector <string> names;
	vector <vector <unsigned char> > sequences;
	vector <int> sitePatternWeights;
	vector <vector <int> > sitePatternRepetitions;	
	vector <int> idsOfVerticesToRemove;
	vector <int> idsOfVerticesToKeep;
	vector <int> idsOfExternalVertices;
	vector <int> idsOfVerticesForSEM;
	vector <tuple <int, string, vector <unsigned char>>> idAndNameAndSeqTupleForVerticesToAdd;	
	//----##############################################################---//
	//	1.	Initialize the global phylogenetic tree T as the empty graph   //
	//----##############################################################---//
//	cout << "Starting MST-backbone" << endl;
//	cout << "1.	Initialize the global phylogenetic tree T as the empty graph" << endl;
	int numberOfInputSequences = (int) this->M->vertexMap->size();	
	this->T = new SEM(1,this->distance_measure_for_NJ,this->verbose);
	// Initialize global phylogeny
	idsOfVerticesForSEM.clear();
	for (pair <int, MST_vertex *> vIdAndPtr : * this->M->vertexMap) {
		idsOfVerticesForSEM.push_back(vIdAndPtr.first);
	}
	tie (names, sequences, sitePatternWeights, sitePatternRepetitions) = this->M->GetCompressedSequencesSiteWeightsAndSiteRepeats(idsOfVerticesForSEM);	
	this->T->sequenceFileName = this->sequenceFileName;
	this->T->AddSequences(sequences);
	this->T->OpenAncestralSequencesFile();
	this->T->AddNames(names);
	this->T->AddSitePatternWeights(sitePatternWeights);
	this->T->SetNumberOfInputSequences(numberOfInputSequences);	
	this->T->numberOfObservedVertices = numberOfInputSequences;
	
	int largestIdOfVertexInMST = numberOfInputSequences;
	
	bool computeLocalPhylogeneticTree = 1;
	bool numberOfNonSingletonComponentsIsGreaterThanZero = 0;	
	
	while (computeLocalPhylogeneticTree) {
		cout << "Number of vertices in MST is " << this->M->vertexMap->size() << endl;
		this->mstBackboneLogFile << "Number of vertices in MST is " << this->M->vertexMap->size() << endl;
		//----####################################################################---//
		//	2.	Compute the size of the smallest subtree ts = (Vs,Es) of M s.t.		 //
		//		|Vs| > s. Check if |Vm\Vs| > s.								   		 //
		// 		If yes then go to step 3 else go to step 9					   		 //
		// 		Bootstrapped alignments may contain zero-weight edges		   		 //
		//      If so then replace |Vs| with |{non-zero weighted edges in Es}| 		 //
		//      Additionally replace |Vm\Vs| with |{non-zero weighted edges in Es}|  //
		//----####################################################################---//				
		computeLocalPhylogeneticTree = this->M->ShouldIComputeALocalPhylogeneticTree();
//		cout << "2. Checking if local phylogenetic tree should be computed" << endl;
		if (computeLocalPhylogeneticTree) {
			//----####################################################################---//
			//	3.	Extract vertices inducing subtree (Vs), and external vertices (Ve)	 //
			//----####################################################################---//				
//			cout << "3. Extract vertices inducing subtree (Vs), and external vertices (Ve)" << endl;
			this->M->SetIdsOfExternalVertices();
			idsOfExternalVertices = this->M->idsOfExternalVertices;
			idsOfVerticesForSEM = this->M->subtree_v_ptr->idsOfVerticesInSubtree;
			this->numberOfVerticesInSubtree = this->M->subtree_v_ptr->idsOfVerticesInSubtree.size();
			for (int id: idsOfExternalVertices) {
				idsOfVerticesForSEM.push_back(id);
			}
			tie (names, sequences, sitePatternWeights, sitePatternRepetitions) = this->M->GetCompressedSequencesSiteWeightsAndSiteRepeats(idsOfVerticesForSEM);
			//----########################################################---//
			//	4.	Compute local phylogeny t over (Vs U Ve) via SEM      	 //
			//----########################################################---//
//			cout << "4.	Compute local phylogeny t over (Vs U Ve) via SEM" << endl;
			this->t = new SEM(largestIdOfVertexInMST,this->distance_measure_for_NJ,this->verbose);						
			this->t->sequenceFileName = this->sequenceFileName;
			this->t->AddSequences(sequences);			
			this->t->SetNumberOfVerticesInSubtree(this->numberOfVerticesInSubtree);			
			this->t->SetNumberOfInputSequences(numberOfInputSequences);			
//			this->t->AddRootVertex();			
			this->t->AddNames(names);
			this->t->numberOfObservedVertices = sequences.size();			
			this->t->AddGlobalIds(idsOfVerticesForSEM);			
			this->t->AddSitePatternWeights(sitePatternWeights);			
			this->t->AddSitePatternRepeats(sitePatternRepetitions);			
			this->t->ComputeNJTree();
			t_start_time = std::chrono::high_resolution_clock::now();
			this->t->RootTreeByFittingAGMMViaEM();
			t_end_time = std::chrono::high_resolution_clock::now();
			this->t->ComputeMAPEstimateOfAncestralSequencesUsingCliques();
//			this->t->OptimizeTopologyAndParametersOfGMM();		
			//----##################################################################---//	
			//  5.	Check if # of non-singleton components of forest f in t that       //
			//		is induced by Vs is greater than zero.							   //
			//		i.e., Does local phylogeny contain vertices/edges of interest?	   //
			//----##################################################################---//
			this->t->SelectIndsOfVerticesOfInterestAndEdgesOfInterest();
			numberOfNonSingletonComponentsIsGreaterThanZero = this->t->IsNumberOfNonSingletonComponentsGreaterThanZero();			
//			cout << "5. Checking if there are any vertices of interest" << endl;
			if (!numberOfNonSingletonComponentsIsGreaterThanZero) {
				//----####################################################---//	
				//  6.	If no then double subtree size and go to step 2 	 //
				//		else reset subtree size and go to to step 7		     //
				//----####################################################---//		
				this->M->DoubleSubtreeSizeThreshold();
//				cout << "6. Doubling subtree size" << endl;
			} else {
				this->M->ResetSubtreeSizeThreshold();		
				//----################################---//
				//  7.	Add vertices/edges in f to T     //
				//----################################---//
//				cout << "7. Adding vertices/edges in f to T" << endl;
				this->t->RenameHiddenVerticesInEdgesOfInterestAndSetIdsOfVerticesOfInterest();
				this->t->SetWeightedEdgesToAddToGlobalPhylogeneticTree();
				this->t->SetAncestralSequencesString();
				this->t->WriteAncestralSequences();
				this->T->AddWeightedEdges(this->t->weightedEdgesToAddToGlobalPhylogeneticTree);
				this->t->SetEdgeAndVertexLogLikelihoods();
				this->T->AddVertexLogLikelihoods(this->t->vertexLogLikelihoodsMapToAddToGlobalPhylogeneticTree);
				this->T->AddEdgeLogLikelihoods(this->t->edgeLogLikelihoodsToAddToGlobalPhylogeneticTree);
				largestIdOfVertexInMST = this->t->largestIdOfVertexInMST;
				//----##############################---//
				//  8.	Update M and go to step 1	   //
				//----##############################---//
//				cout << "8. Updating MST" << endl;
				this->t->SetInfoForVerticesToAddToMST();
				this->M->UpdateMSTWithMultipleExternalVertices(t->idsOfVerticesToKeepInMST, t->idsOfVerticesToRemove, t->idAndNameAndSeqTuple, idsOfExternalVertices);				
			}			
			computeLocalPhylogeneticTree = this->M->ShouldIComputeALocalPhylogeneticTree();
			cout << "CPU time used for computing local phylogeny is " << chrono::duration<double>(t_end_time-t_start_time).count() << " second(s)\n";
			this->mstBackboneLogFile << "CPU time used for computing local phylogeny is " << chrono::duration<double>(t_end_time-t_start_time).count() << " second(s)\n";			
		}
	}	
	//----########################################################---//
	//	9.	Compute phylogenetic tree t over vertices in M, and      //
	//		add vertices/edges in t to T							 //
	//----########################################################---//
	cout << "Computing phylogenetic tree over remaining vertices" << endl;	
	idsOfVerticesForSEM.clear();
	for (pair <int, MST_vertex *> idPtrPair: * this->M->vertexMap) {
		idsOfVerticesForSEM.push_back(idPtrPair.first);
	}
	cout << "Number of vertices in MST is " << idsOfVerticesForSEM.size() << endl;
	this->mstBackboneLogFile << "Number of vertices in MST is " << this->M->vertexMap->size() << endl;
	tie (names, sequences, sitePatternWeights, sitePatternRepetitions) = this->M->GetCompressedSequencesSiteWeightsAndSiteRepeats(idsOfVerticesForSEM);
	this->numberOfVerticesInSubtree = sequences.size();
	this->t = new SEM(largestIdOfVertexInMST,this->distance_measure_for_NJ,this->verbose);
	this->t->SetFlagForFinalIterationOfSEM();
	this->t->sequenceFileName = this->sequenceFileName;
	this->t->AddSequences(sequences);
	this->t->SetNumberOfVerticesInSubtree(this->numberOfVerticesInSubtree);
	this->t->SetNumberOfInputSequences(numberOfInputSequences);
	this->t->numberOfObservedVertices = sequences.size();
//	this->t->AddRootVertex();
	this->t->AddNames(names);
	this->t->AddGlobalIds(idsOfVerticesForSEM);
	this->t->AddSitePatternWeights(sitePatternWeights);
	this->t->AddSitePatternRepeats(sitePatternRepetitions);
	this->t->ComputeNJTree();
	t_start_time = std::chrono::high_resolution_clock::now();
	this->t->RootTreeByFittingAGMMViaEM();
	t_end_time = std::chrono::high_resolution_clock::now();
	this->t->ComputeMAPEstimateOfAncestralSequencesUsingCliques();
	this->t->SelectIndsOfVerticesOfInterestAndEdgesOfInterest();
	this->t->RenameHiddenVerticesInEdgesOfInterestAndSetIdsOfVerticesOfInterest();
	this->t->SetWeightedEdgesToAddToGlobalPhylogeneticTree();
	this->t->SetAncestralSequencesString();
	this->t->WriteAncestralSequences();
	this->T->AddWeightedEdges(this->t->weightedEdgesToAddToGlobalPhylogeneticTree);	
	this->t->SetEdgeAndVertexLogLikelihoods();
	this->T->AddVertexLogLikelihoods(this->t->vertexLogLikelihoodsMapToAddToGlobalPhylogeneticTree);
	this->T->AddEdgeLogLikelihoods(this->t->edgeLogLikelihoodsToAddToGlobalPhylogeneticTree);
	cout << "CPU time used for computing local phylogeny is " << chrono::duration<double>(t_end_time-t_start_time).count() << " second(s)\n";
	this->mstBackboneLogFile << "CPU time used for computing local phylogeny is " << chrono::duration<double>(t_end_time-t_start_time).count() << " second(s)\n";			
	// assert that T is a tree
//	cout << "Number of vertices in T is " << this->T->vertexMap->size() << endl;
//	cout << "Number of edges in T is " << this->T->edgeLengths.size() << endl;
	assert(this->T->vertexMap->size() == this->T->edgeLengths.size() + 1);
	//----##############---//		
	//	10.	Root T via EM  //
	//----##############---//
	current_time = std::chrono::high_resolution_clock::now();
	cout << "CPU time used for computing unrooted topology is " << chrono::duration<double>(current_time-start_time).count() << " second(s)\n";
	this->mstBackboneLogFile << "CPU time used for computing unrooted topology is " << chrono::duration<double>(current_time-start_time).count() << " second(s)\n";
//	cout << "Rooting T via EM" << endl;
//	this->T->RootTreeByFittingAGMMViaEM();
	cout << "Rooting T by maximizing expected log likelihood" << endl;
	this->T->RootTreeBySumOfExpectedLogLikelihoods();
}


void MSTBackbone::MSTBackboneWithOneExternalVertex() {
	this->T = new SEM(1,this->distance_measure_for_NJ,this->verbose);
//	ofstream edgeListFile;
//	edgeListFile.open(this->sequenceFileName + ".edgeList");
	cout << "Starting MST-backbone" << endl;
//	int numberOfInputSequences = (int) this->M->vertexMap->size();
	bool subtreeExtractionPossible = 1;		
	vector <string> names;
	vector <vector <unsigned char> > sequences;
	vector <int> sitePatternWeights;
	vector <vector <int> > sitePatternRepetitions;
	vector <int> idsOfVerticesForSEM;
	vector <int> idsOfVerticesToRemove;
	vector <int> idsOfVerticesToKeep;
	vector <int> idsOfExternalVertices;
	vector <tuple <int, string, vector <unsigned char>>> idAndNameAndSeqTupleForVerticesToAdd;

	int h_ind = 1;
	
	// Initialize global phylogeny
	idsOfVerticesForSEM.clear();
	for (pair <int, MST_vertex *> vIdAndPtr : * this->M->vertexMap) {
		idsOfVerticesForSEM.push_back(vIdAndPtr.first);
	}
	tie (names, sequences, sitePatternWeights, sitePatternRepetitions) = this->M->GetCompressedSequencesSiteWeightsAndSiteRepeats(idsOfVerticesForSEM);	
	this->T->sequenceFileName = this->sequenceFileName;
	this->T->AddSequences(sequences);
	this->T->AddNames(names);
	this->T->AddSitePatternWeights(sitePatternWeights);
	cout << "Number of leaves is " << this->T->vertexMap->size() << endl;
	
	int numberOfRemainingVertices;
	int numberOfVerticesInSubtree;
	vector <string> weightedEdges;	
	string u_name; string v_name;
	vector <unsigned char> sequenceToAdd;
	string nameOfSequenceToAdd;
//	int totalNumberOfEdges = 0;
//	int vertex_id = numberOfInputSequences;

	vector <unsigned char> seq_u; vector <unsigned char> seq_v;
	vector <unsigned char> compressed_seq_u; vector <unsigned char> compressed_seq_v;
//	int u_ind; int v_ind;
//	int u_id; int v_id;	
	map <int, int> EMVertexIndToPhyloVertexIdMap;	
//	SEM_vertex * v_phylo;		
//	bool resetSubtreeSizeThreshold = 1;
	int subtreeSizeThreshold = this->M->numberOfLargeEdgesThreshold;
	cout << "Subtree size threshold is " << subtreeSizeThreshold << endl;
	// Iterate to completion
	MST_vertex * v_mst;
	tie (subtreeExtractionPossible, v_mst) = this->M->GetPtrToVertexSubtendingSubtree();
	numberOfRemainingVertices = this->M->vertexMap->size() - v_mst->idsOfVerticesInSubtree.size();
	cout << "numberOfRemainingVertices is " << numberOfRemainingVertices << endl;
	if (v_mst->idOfExternalVertex == -1 or numberOfRemainingVertices < 3) {
		subtreeExtractionPossible = 0;
	}
	cout << "Sequence length is " << this->T->sequenceLength << endl;
	while (subtreeExtractionPossible) {
		cout << "Number of vertices in MST is " << this->M->vertexMap->size() << endl;
		this->t = new SEM(h_ind,this->distance_measure_for_NJ,this->verbose);	
		// ids of vertices in subtree
		idsOfVerticesForSEM = v_mst->idsOfVerticesInSubtree;
		numberOfVerticesInSubtree = v_mst->idsOfVerticesInSubtree.size();
		// ids of external vertices
		idsOfVerticesForSEM.push_back(v_mst->idOfExternalVertex);
		tie (names, sequences, sitePatternWeights, sitePatternRepetitions) = this->M->GetCompressedSequencesSiteWeightsAndSiteRepeats(idsOfVerticesForSEM);
		h_ind += sequences.size() -2;
		// Perform SEM
		this->t->SetNumberOfVerticesInSubtree(numberOfVerticesInSubtree);
		this->t->AddSequences(sequences);
		this->t->AddRootVertex();
		this->t->AddNames(names);
		this->t->AddSitePatternWeights(sitePatternWeights);		
		this->t->OptimizeTopologyAndParametersOfGMM(); // use hard EM + soft EM?			
		// Get edges to add
		
		this->T->AddWeightedEdges(this->t->weightedEdgesToAddToGlobalPhylogeneticTree);						
		sequenceToAdd = DecompressSequence(&this->t->compressedSequenceToAddToMST, &sitePatternRepetitions);			
		// edgeListFile << this->t->weightedEdgeListString;			
		// Update MST
		
		this->M->UpdateMSTWithOneExternalVertex(v_mst->idsOfVerticesInSubtree, this->t->nameOfSequenceToAddToMST, sequenceToAdd);
		tie (subtreeExtractionPossible, v_mst) = this->M->GetPtrToVertexSubtendingSubtree();
		numberOfRemainingVertices = this->M->vertexMap->size() - v_mst->idsOfVerticesInSubtree.size();
		if (v_mst->idOfExternalVertex == -1 or numberOfRemainingVertices < 3) {
			subtreeExtractionPossible = 0;
		}
		delete this->t;	
	}		
	cout << "Number of remaining vertices in MST is " << this->M->vertexMap->size() << endl;		
	idsOfVerticesForSEM.clear();
	for (pair <int, MST_vertex *> vIdAndPtr : * this->M->vertexMap) {
		idsOfVerticesForSEM.push_back(vIdAndPtr.first);
	}	
	this->t = new SEM(h_ind,this->distance_measure_for_NJ,this->verbose);
	tie (names, sequences, sitePatternWeights, sitePatternRepetitions) = this->M->GetCompressedSequencesSiteWeightsAndSiteRepeats(idsOfVerticesForSEM);
//	cout << "Number of distinct site patterns is " << sitePatternWeights.size() << endl;
	this->t->SetFlagForFinalIterationOfSEM();
	this->t->numberOfExternalVertices = 1;
	this->t->AddSequences(sequences);
	this->t->AddRootVertex();
	this->t->AddNames(names);
	this->t->AddSitePatternWeights(sitePatternWeights);
	t_start_time = std::chrono::high_resolution_clock::now();
	this->t->OptimizeTopologyAndParametersOfGMM();
	t_end_time = std::chrono::high_resolution_clock::now();
	this->T->AddWeightedEdges(this->t->weightedEdgesToAddToGlobalPhylogeneticTree);
//	if (weightedEdges.size() != this->t->weightedEdgesToAddToGlobalPhylogeneticTree.size()) {
//		cout << "Number of edges print to file is " << weightedEdges.size() << endl;
//		cout << "Number of edges added to global phylogenetic tree is " << this->t->weightedEdgesToAddToGlobalPhylogeneticTree.size() << endl;	
//	}
	delete this->t;	
	// Root T using edge loglikelihoods
//	this->T->PerformModelSelection();
	// Contract leaf incident edges s.t. AIC is minimized
//	this->T->WriteTree();
	// Write model parameters, and rate categories
	// Write tree in edge list and newick format	
//	cout << "Number of vertices in global phylogeny is " << T->vertices.size() << endl;
//	cout << "Number of edges in global phylogeny is " << T->edgeLengths.size() << endl;
	// Write rooted tree in edge list format (analysis on simulated data)
	// Write rooted tree in newick format	
	// Perform model selection
}

string MSTBackbone::GetSequenceListToWriteToFile(map <string, vector <unsigned char>> compressedSeqMap, vector <vector <int> > sitePatternRepetitions) {	
	vector <unsigned char> decompressedSequence;
	string dnaSequence;
	string u_name;
	string listOfVertexNamesAndDNAsequencesToWriteToFile;	
	for (pair <string,vector<unsigned char>> nameSeqPair : compressedSeqMap) {		
		decompressedSequence = DecompressSequence(&(nameSeqPair.second),&sitePatternRepetitions);
		dnaSequence = EncodeAsDNA(decompressedSequence);
		listOfVertexNamesAndDNAsequencesToWriteToFile += nameSeqPair.first + "\t" + dnaSequence + "\n"; 		
	}	
	return (listOfVertexNamesAndDNAsequencesToWriteToFile);
}

int MSTBackbone::GetEdgeIndex (int vertexIndex1, int vertexIndex2, int numberOfVertices){
	int edgeIndex;
	edgeIndex = numberOfVertices*(numberOfVertices-1)/2;
	edgeIndex -= (numberOfVertices-vertexIndex1)*(numberOfVertices-vertexIndex1-1)/2;
	edgeIndex += vertexIndex2 - vertexIndex1 - 1;
	return edgeIndex;
}

int MSTBackbone::ComputeHammingDistance(string seq1, string seq2) {
	int hammingDistance = 0;
	for (unsigned int i=0;i<seq1.length();i++){
		if (seq1[i] != seq2[i]){
			hammingDistance+=1;
		}		
	}
	return (hammingDistance);
};

int MSTBackbone::ComputeHammingDistance(vector<unsigned char> recodedSeq1, vector<unsigned char> recodedSeq2) {
	int hammingDistance = 0;
	float ungappedSequenceLength = 0;
	for (unsigned int i=0;i<recodedSeq1.size();i++) {
		if (recodedSeq1[i] != recodedSeq2[i]) {
			hammingDistance+=1;
		}		
	}	
	return (hammingDistance);
};


#endif



