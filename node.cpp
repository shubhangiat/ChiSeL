/*** Implements the node.h class, for representing the nodes of the Query graph
 *** Is also the base class for Vertex class defined in vertex.h
***/

#include "node.h"
#include <iostream>

// Constructs a vertex with the corresponding characteristics
// Sets the ID and the label of a vertex
Node :: Node(string id, string lab) : ID(id), label(lab)
{
	neighbourIDs.clear();
}


// Deallocates space
Node :: ~Node()
{
	neighbourIDs.clear();
}


// Returns vertex label
const string Node :: getLabel(void) const
{
	return label;
}


// Returns vertex ID
const string Node :: getID(void) const
{
	return ID;
}


// Returns the ID list of the neighbouring vertices
const set<string>& Node :: getNeighbourIDs(void) const
{
	return neighbourIDs;
}


// Add the IDs of adjacent vertices
void Node :: addNeighbours(set<string>& neigh)
{
	set<string>::iterator it = neigh.begin();
	for(; it != neigh.end(); it++)
		neighbourIDs.insert(*it);
}


// Add the ID of one adjacent vertex
void Node :: addNeighbours(string neigh)
{
	neighbourIDs.insert(neigh);
}


// Prints the node characteristics
void Node :: print(void) const
{
	cout<<"ID: "<<ID<<"\tLabel: "<<label<<"\tNeighbourIDs: ";

	set<string>::iterator it = neighbourIDs.begin();
	for(; it!=neighbourIDs.end(); it++)
		cout<<*it<<" ";

	cout<<endl;
}
