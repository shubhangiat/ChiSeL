/*** Implements the query graph (of query_graph.h) and constructs the index as required
***/


#include "query_graph.h"
#include <utility>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <iostream>


// Creates the query graph from two files - (1) mapping of node ID to label, (2) list of neighbour labels
Query :: Query(const string iNodeFile, const string iEdgeFile)
{
	string prevId, id, label, edge;  // Elements to be read from the files
	prevId = id = label = edge = "";

	set<string> neighbours;  // Stores the ID of the neighbours for a vertex

	ifstream iNF;
	iNF.open(iNodeFile.c_str());  // Contains the mapping from vertex IDs to labels

	if(!iNF.is_open())
		cout<<"Unable to open "<<iNodeFile<<" file"<<endl;


	iNF >> id >> label;
	while(!iNF.eof())  // Reads the vertex IDs and the labels
	{
		Node *node = new Node(id, label);
		graph[id] = node;

		iNF >> id >> label;
	}

	iNF.close();


	ifstream iEF;
	iEF.open(iEdgeFile.c_str());  // Contains the edges between the vertices

	if(!iEF.is_open())
		cout<<"Unable to open input "<<iEdgeFile<<" file"<<endl;

	unsigned int file_empty_flag = 1;	// Initialized to true
	char c;

	iEF.get(c);
	while(!iEF.eof())	// Checking for empty file
	{
		if(!isspace(c))
		{
			file_empty_flag = 0;	// Set to false
			break;
		}
	}

	if(!file_empty_flag)	// If file not empty
	{
		iEF.seekg(0, iEF.beg);
		iEF >> id >> edge;  // Reads a new node ID and its label
		while(!iEF.eof())  // Read the neighbours for vertices
		{
			#ifdef DEBUG
				cout<<id<<" "<<edge<<endl;
			#endif
			
			if(prevId.compare(id) == 0)
				neighbours.insert(edge);
			else
			{
				if (prevId != "")
				{
					Node *node = graph[prevId];
					node->addNeighbours(neighbours);
					neighbours.clear();
				}

				neighbours.insert(edge);
				prevId = id;
			}

			// Add the reverse edges for undirected graph
			Node *node = graph[edge];
			node->addNeighbours(id);

			iEF >> id >> edge;  // Reads a new node ID and its label
		}

		Node *node = graph[prevId];
		node->addNeighbours(neighbours);
		neighbours.clear();
	}
	iEF.close();


	createLabelNeighbours();  // Creates the neighbourLabels index
}


// Deallocates the graph
Query :: ~Query()
{
	unordered_map<string, Node*>::iterator it = graph.begin();
	for(; it!=graph.end(); it++)
		delete it->second;

	graph.clear();

	labelNeighbours.clear();
}


// Return the address of the graph
// CRITICAL as it could result in accidental modification
const unordered_map<string, Node*>&  Query :: getGraph_CRITICAL() const
{
	return graph;
}


// Returns the labelNeighbours of the vertex label provided as argument minus that of visited query vertices
const unordered_map<string, vector<string> > Query :: getLabelNeighbours(set<string> visited, string lab) const
{
	unordered_map<string, vector<string> > listOfLabelNeigh;

	// If vertex with same label not present
	// Define dummy pair<empty_label, q_vertex_label> to search in multimap
	if(labelNeighbours.find(make_pair(EMPTY_LABEL, lab)) == labelNeighbours.end())
		return listOfLabelNeigh;

	std::pair< multimap< pair<string,string>, vector<string> >::const_iterator, multimap< pair<string,string>, vector<string> >::const_iterator> ret;
	ret = labelNeighbours.equal_range(make_pair(EMPTY_LABEL, lab));

	for(multimap< pair<string,string>, vector<string> >::const_iterator it=ret.first; it!=ret.second; it++)
	{
		// If query vertex is not in the visited list
		if(find(visited.begin(), visited.end(), (it->first).first) == visited.end())
			listOfLabelNeigh.insert(make_pair((it->first).first, it->second));
	}

	return listOfLabelNeigh;
}


// Returns the unique query graph labels
const set<string> Query :: getLabels(void) const
{
	set<string> labels;
	multimap< pair<string, string>, vector<string> >::const_iterator it = labelNeighbours.begin();

	for(; it!=labelNeighbours.end(); it++)
		labels.insert((it->first).second);

	return labels;
}


// Populates the index variable labelNeighbours, mapping vertices to the labels of the neighbours
void Query :: createLabelNeighbours(void)
{
	unordered_map<string, Node*>::iterator it = graph.begin();

	for(; it!=graph.end(); it++)
	{
		string currID = (it->second)->getID();
		string currLab = (it->second)->getLabel();
		const set<string> neighId = (it->second)->getNeighbourIDs();
		set<string>::const_iterator it1 = neighId.begin();

		vector<string> nLabels;
		nLabels.clear();

		for(; it1!=neighId.end(); it1++)
			nLabels.push_back( (graph[*it1])->getLabel() );

		labelNeighbours.insert( pair< pair<string, string>, vector<string> >(make_pair(currID, currLab), nLabels) );
	}
}


// Returns the number of nodes in the query graph
unsigned Query :: getGraphSize(void) const
{
	return graph.size();
}


// Prints the vertex IDs of the query graph
void Query :: printVertex(void) const
{
	unordered_map<string, Node*>::const_iterator it = graph.begin();

	for(; it!=graph.end(); it++)
		cout<<it->first<<" ("<<(it->second)->getLabel()<<"), ";

	cout<<endl;
}


// Prints the graph characteristics
void Query :: printGraph(void) const
{
	cout<<"\nThe query graph contains the following vertices:"<<endl;

	unordered_map<string, Node*>::const_iterator it = graph.begin();
	for(; it!=graph.end(); it++)
		(it->second)->print();

	cout<<"\nThe index of the vertex and the neighbour labels are:"<<endl;
	multimap< pair<string, string>, vector<string> >::const_iterator it1 = labelNeighbours.begin();
	for(; it1!=labelNeighbours.end(); it1++)
	{
		cout<<"("<<(it1->first).first<<", "<<(it1->first).second<<") has neighbour label: [";
		if(it1->second.size()!=0)
		{
			copy((it1->second).begin(),(it1->second).end()-1, ostream_iterator<string>(cout, ", " ));
			cout<<(it1->second).back();
		}
		cout<<"]"<<endl;
	}
	cout<<endl;
}
