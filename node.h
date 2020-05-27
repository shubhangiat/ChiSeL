/*** Defines the basic vertex structure of a graph comprising an unique ID, corresponding label,
 *** and a set of neighbouring node labels. These vertices define the Query graph, and will be
 *** extended for representing the input graph vertices.
***/


#ifndef NODE_H_
 #define NODE_H_


#include <string>
#include <set>

using namespace std;


class Node
{
   protected:
	const string ID;  // Stores the unique ID of the vertex
	const string label;  // Stores the associated label with the vertex

	set<string> neighbourIDs;  // Stores the ID of the adjacent vertices

   public:
	// Constructs a vertex with the corresponding characteristics
	Node(const string, const string);  // Sets the ID and the label of a vertex

	virtual ~Node();  // Deallocates space

	void addNeighbours(set<string>&);  // Add the IDs of adjacent vertices in batch
	void addNeighbours(string);  // Add the ID of an adjacent vertex, one at a time

	const string getLabel(void) const;  // Returns vertex label
	const string getID(void) const;  // Returns vertex ID
	const set<string>& getNeighbourIDs(void) const;  // Returns the ID list of the neighbouring vertices

	virtual void print(void) const;  // Prints the node characteristics
};

#endif
