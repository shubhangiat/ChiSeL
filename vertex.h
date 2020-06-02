/*** Defines the enhanced vertices of the input graph by extending the basic node structure (in node.h)
 *** This provides the vertices with a chi-square value and the associated symbols and incorporates the
 *** uncertainty information associated with edges.
***/


#ifndef VERTEX_H_
 #define VERTEX_H_


#include "const.h"
#include "node.h"
#include "query_graph.h"
#include <vector>
//#include <unordered_map>
#include <algorithm>

using namespace std;

class Edge	// To facilitate storing of edge probabilities
{
   private:
	string nodeID;		// ID of neighbour vertex
	string neighLabel;	// Label of neighbour vertex
    double probability;	// Probability of being neighbour, i.e., existential probability of edge

   public:
   
	//Stores the neighbor with the associated probability
    Edge(const string ID, const string lab, const double pr) : nodeID(ID), neighLabel(lab), probability(pr) {}
    ~Edge() {}  //Deallocates space

    string getID() const { return nodeID; }		// Returns neighbour ID

    string getLabel() const { return neighLabel; }    // Returns neighboring label

    double getProbability() { return probability; } // Returns probability associated with the edge

	// Sets the probability of edge between the two vertices
	// CRITICAL since can modify graph by changing the edge probabilities of graph
    void setProbability_CRITICAL(double pr)	{ probability = pr; }
};

class Vertex : public Node  // publicly inherits the node structure
{
   private:

	double symbolOccurrence[NUM_OF_SYMBOLS];  // Stores the observed values of the symbols s0, s1, s2 (in that order)
	double symbolOccProbability[NUM_OF_SYMBOLS];  // Stores the expected probability of the symbols s0, s1, s2 (in that order)
	double degree;		// expected degree of the vertex

	unsigned long chiLength;  // Stores the number of query triplets in mapped query vertex
	long double chiValue;  // Stores the chi-square value of the vertex
	vector<string> chiLabels;  // Stores the query labels providing the maximum chi sq value

	// Stores the information if the qVertexID value changed from its first computation for primary heap
	// 0 implies no change 1 otherwise
	bool qVertexID_dirty_bit;
	string qVertexID;	// Stores the query vertex ID to which this vertex is mapped

	vector<Edge*> neighbourLabels;  // Stores the adjacent vertices

	pair<string, double> parent;	// Stores parent vertexID and associated edge probability

	unordered_map<string, vector<double> > pr_lx;

	// pr_lx[label][0] Probability that no instance of label l_x exists (= product of non existence of all instances of 'label' as neighbour),
	// pr_lx[label][1] exactly one instance of label l_x exists (= sum of product of existence of exactly one instance of 'label' as neighbour),
	// pr_lx[label][2] atleast one instance of label l_x exists (= (1 - probability that no instance of 'label' exists) = (1 - pr_lx[label][0]) )

   public:

	Vertex(const string, const string);  // Creates a vertex with the characteristics
	~Vertex();  // Deallocates space

	const vector<string> getNeighbourLabels(void) const;  // Returns the list of neighbour labels

	long double getChiValue(void) const;  // Returns the chi-square value associated with the vertex

	const vector<string>& getChiLabels(void) const;  // Returns the chiLabels associated with this vertex

	const bool get_qVertexID_dirty_bit(void) const; // Returns the value of qVertexID_dirty_bit

	void set_qVertexID_dirty_bit(bool); // Sets the value of qVertexID_dirty_bit

	const string get_qVertexID(void) const;	//	Returns the query vertex the vertex is mapped to

	const double* getSymbolOcc(void) const;  // Returns the observed value of occurrences of the symbols s0, s1, s2

	// Computes the chi-square value of the vertex using the symbolOccurrence[]
	void computeChiSqValue(const unordered_map<string, vector<string> >&, bool);

	void computeDegree(void);   // Computes the degree of the vertex given probabilistic neighbours

	// Computes (non)existence probabilities of each unique neighbour label
	// (zero, exactly one and at least one instance exists)
	void compute_pr_lx(void);

	void compute_symOccPr(unsigned); // Computes expected distribution of symbols 

    virtual void print(void) const;  // Prints the vertex characteristics

    void addNeighbourLabels(Edge *);	// Add neighbor labels with associated probabilities

    void setParent(const Vertex *);		// Set parent with parent vertexID and associated edge probability

	const pair<string, double> getParent(void) const;		// Return parent <parent vertexID, associated edge probability> pair

	Edge* getNeighbour(const string) const;		// Return edge with passed neighbourID

};



/*** Inherited data members and methods (from node.h)
        const string ID;
        const string label;  // Stores the associated label with the vertex
        set<string> neighbourIDs;  // Stores the ID of the adjacent vertices

        Node(const string, const string);  // Sets the ID and the label of a vertex
        void addNeighbours(set<string>&);  // Add the IDs of adjacent vertices
        void addNeighbours(string);  // Add the ID of one adjacent vertices
        const string getLabel(void) const;  // Returns vertex label
        const string getID(void) const;  // Returns vertex ID
        const set<string>& getNeighbourIDs(void) const;  // Returns the ID list of the neighbouring vertices
***/


#endif
