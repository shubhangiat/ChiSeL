/*** Defines the query graph as a collection of vertices (of node.h). Also keeps an
 *** index unordered_mapping the vertex labels to neighbouring labels
***/


#ifndef QUERY_GRAPH_H_
 #define QUERY_GRAPH_H_

#include "const.h"
#include "node.h"
#include <utility>
#include <unordered_map>
#include <vector>
#include <map>

// To sort on the second element of the pair
class Compare_lab
{
	public:
	 bool operator() (const pair<string, string> &a, const pair<string, string> &b) const
	 {
		return(a.second<b.second);
	 }
};

class Query
{
   private:
	unordered_map<string, Node*> graph;  // Stores the query graph as an unordered_map of vertex ID to the corresponding object

	// Provides an index to map a vertex label to the labels of all its neighbouring vertices
	multimap< pair<string, string>, vector<string>, Compare_lab > labelNeighbours;

	// POpulates the above index
	void createLabelNeighbours(void);

   public:
	// Creates the query graph from two files - (1) mapping of node ID to label, (2) list of neighbour labels
	Query(const string, const string);

	~Query();  // Deallocates the graph

	// Returns the labelNeighbours to the vertex label provided as argument minus that of visited query vertices
	const unordered_map<string, vector<string> > getLabelNeighbours(set<string>, string) const;

	// Return the address of the graph of query graph
	// CRITICAL as it could result in accidental modification
	const unordered_map<string, Node*>& getGraph_CRITICAL() const;

	const set<string> getLabels(void) const;  // Returns the unique query graph labels
	unsigned getGraphSize(void) const;  // Returns the number of nodes in the query graph
	void printVertex(void) const;  // Prints the vertex IDs of the query graph
	void printGraph(void) const;   // Prints the graph characteristics
};


#endif
