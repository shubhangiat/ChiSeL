/*** Defines the input graph as a collection of vertices (of vertex.h). Also computes
 *** the output approximate sub-graph matching to the query provided as input.
***/


#ifndef INPUT_GRAPH_H_
 #define INPUT_GRAPH_H_


#include "vertex.h"
#include "query_graph.h"
#include <queue>
#include <unordered_map>



using namespace std;


// Compare operator for min-heap operation
class Compare_min
{
   public:
	bool operator() (Vertex *a, Vertex *b)
	{
		double a_val = a->getChiValue();
		double b_val = b->getChiValue();

		return (a_val > b_val);
	}
};


// The secondary heap in input_graph.cpp uses this comparator
// The edge probability with "parent" considered while inserting in the heap
class Compare_max
{
   public:
	bool operator() (Vertex *a, Vertex *b)
	{
		double a_val = a->getChiValue()*a->getParent().second;
		double b_val = b->getChiValue()*b->getParent().second;

		return (a_val < b_val);
	}
};



class Input
{
   private:
	// Stores the input graph as an unordered_map of vertex ID to the corresponding object
	unordered_map<string, Vertex*> graph;

	// Stores the vertices having the same labels
	unordered_map<string, vector<Vertex*> > vertexLabel;

	// Min-heap structure to keep the vertices with top-k chi-sq values (for greedy approach)
	// The number of elements in the heap is constrained by the ORDER_CONSTANT 
	// Minimum value in the heap is (the value at the top of the heap) is required
	// This value is removed and a new higher value is inserted in the heap when encountered
	// Even with Compare_min as comparator the values in the heap stored are the maximum values
	priority_queue<Vertex*, vector<Vertex*>, Compare_min> heap;

	// Stores unique labels of the input graph
	set<string> uniq_labs;

   public:
	// Creates the input graph from two input files - (1) mapping of node ID to label, (2) list of neighbour labels
	Input(const string, const string);
	~Input();  // Deallocates the constructed graph

	// Returns the vertex IDs of the top-k matching Subgraphs (to the provided query) found in the input graph
	vector<vector<Vertex*> > getSubGraphs(const Query&);

	// Returns the number of vertices in the graph
	unsigned getGraphSize(void) const;

	// Prints the graph characteristics
	void printGraph(void) const;

	// Add label to the unique label set
	void addLabel(string);

	// Get the number of uniqLabels present in the graph, used to compute the chi square values of the vertices
	unsigned getUniqLabsSize(void) const;

	// Perturb the input graph and set the probability of the edges of query subgraph present in input graph as 1
	// Returns the original probability of the edges for those the perturbation was done
	// CRITICAL function, modifies graph
	// The query vertex id must match graph vertex id
	map<pair<string, string>, double> perturb_CRITICAL(const Query&);

	// Unperturb the input graph using the edge probabilities of the original graph
	// CRITICAL function, modifies graph
	void unperturb_CRITICAL(map<pair<string, string>, double>);
};


#endif
