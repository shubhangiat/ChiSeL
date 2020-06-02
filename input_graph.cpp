/*** Implements the class Input (input_graph.h) and its functionalities
***/


#include "input_graph.h"
#include <fstream>
#include <iostream>
#include <cstdlib>


// Creates the input graph from two input files - (1) mapping of node ID to label, (2) list of neighbour labels
Input :: Input(const string iNodeFile, const string iEdgeFile)
{
    // Elements to be read from the files
    string prevId, id, label, edge;
    prevId = id = label = edge = "";
    double pr = 0; //To be read from edge file

	ifstream iNF;
	iNF.open(iNodeFile.c_str());  // Contains the mapping from vertex IDs to labels

	if(!iNF.is_open())
	{
		cout<<"Unable to open "<<iNodeFile<<" file for input"<<endl;
		exit(0);
	}

	iNF >> id >> label;
	while(!iNF.eof())  // Reads the vertex IDs and the labels
	{
		#ifdef DEBUG
				cout<<id<<" "<<label<<endl;
		#endif	// DEBUG
		Vertex *ver = new Vertex(id, label);
		graph[id] = ver;

		// Add label to unique label set
		addLabel(label);

		// Storing which vertices ahre labels
		if(vertexLabel.find(label) == vertexLabel.end())
		{
			vector<Vertex*> v;
			v.push_back(ver);
			vertexLabel[label] = v;
		}
		else
			(vertexLabel[label]).push_back(ver);

		iNF >> id >> label;
	}

	iNF.close();

	#ifdef DEBUG
		cout<<"Node file read"<<endl;
	#endif // DEBUG

	ifstream iEF;
	iEF.open(iEdgeFile.c_str());  // Contains the edges between the vertices

	set<string> neighbours;  // Stores the ID of the neighbours for a vertex

	if(!iEF.is_open())
	{
		cout<<"Unable to open "<<iEdgeFile<<" file for input"<<endl;
		exit(0);
	}

	iEF >> id >> edge >> pr;  // Reads a new node ID and its label
	while(!iEF.eof())  // Read the neighbours for vertices
	{
		#ifdef DEBUG
				cout<<id<<" "<<edge<<" "<<pr<<endl;
		#endif	// DEBUG
		if(prevId.compare(id) == 0)
		{
			// If the previous edge that was read had first vertex same as this edge
			// Then simply add the neighbour vertex to the batch variable 'neighbours'
			neighbours.insert(edge);
		}
		else
		{
			// Add the neighbours to the neighbour list (in batch)
			if (prevId != "")
			{
				Vertex *ver = graph[prevId];
				ver->addNeighbours(neighbours);
				neighbours.clear();
			}

			neighbours.insert(edge);
			prevId = id;
		}

		// Storing neighbour label with associated neighbor probabilities
		Edge *e = new Edge(edge, graph[edge]->getLabel(), pr);
		Vertex *ver = graph[prevId];
		ver->addNeighbourLabels(e);

		// Add the reverse edges for undirected graph
		ver = graph[edge];
		ver->addNeighbours(id);
		e = new Edge(id, graph[id]->getLabel(), pr);
		ver->addNeighbourLabels(e);

		iEF >> id >> edge >> pr;  // Reads a new node ID and its label
	}

	Vertex *ver = graph[prevId];
	ver->addNeighbours(neighbours);
	neighbours.clear();

	iEF.close();

	// Compute number of unique labels
	unsigned card = getUniqLabsSize();

	// Compute vertex degrees, (non)existence probabilities of neighbour labels
	// and expected symbol occurrence probability for each graph vertex
	unordered_map<string, Vertex*>::const_iterator it = graph.begin();
	for(; it!=graph.end(); it++)
	{
		it->second->computeDegree();
		it->second->compute_pr_lx();
		it->second->compute_symOccPr(card);
	}
}


// Deallocates the constructed graph
Input :: ~Input()
{
	unordered_map<string, Vertex*>::iterator it = graph.begin();
	for(; it!=graph.end(); it++)
		delete it->second;

	graph.clear();

	while(!heap.empty())
		heap.pop();

	vertexLabel.clear();
	uniq_labs.clear();
}


// Perturb the input graph and set the probability of the edges of query subgraph present in input graph as 1
// Returns the original probability of the edges for those the perturbation was done
// CRITICAL function, modifies graph
// The query vertex id must match graph vertex id
map<pair<string, string>, double> Input :: perturb_CRITICAL(const Query& qry)
{
    map<pair<string, string>, double> orig_epr;		// to store original edge probabilities

    // Loop over the query graph vertices (to loop over the edges)
    unordered_map<string, Node*>::const_iterator qg_itr = qry.getGraph_CRITICAL().begin();
    for(; qg_itr!=qry.getGraph_CRITICAL().end(); qg_itr++)
    {
        // For each neighbour of qg_itr store the original contents in orig_epr and perturb the probability to 1 in graph
        set<string> neighIDs = qg_itr->second->getNeighbourIDs();
        set<string>::const_iterator nid_itr = neighIDs.begin();
        for(; nid_itr != neighIDs.end(); nid_itr++)
        {
            if(graph[qg_itr->first]->getNeighbour(*nid_itr))
			{
				// If the vertices are connected in the input target graph else ignore
				// Condition is useful for noisy queries
				// store the info in orig_epr
				double pr = graph[qg_itr->first]->getNeighbour(*nid_itr)->getProbability();
				orig_epr[make_pair(qg_itr->first, *nid_itr)] = pr;

				// perturb the graph
				graph[qg_itr->first]->getNeighbour(*nid_itr)->setProbability_CRITICAL(1);
			}
        }
    }

    // Re-compute number of unique labels
	unsigned card = getUniqLabsSize();
	
	// Re-compute vertex degrees, (non)existence probabilities of neighbour labels 
	// and expected symbol occurrence probability
 	for(qg_itr = qry.getGraph_CRITICAL().begin(); qg_itr!=qry.getGraph_CRITICAL().end(); qg_itr++)
    {
		graph[qg_itr->first]->computeDegree();
		graph[qg_itr->first]->compute_pr_lx();
		graph[qg_itr->first]->compute_symOccPr(card);
	}

	return orig_epr;
}


// Unperturb the input graph using the edge probabilities of the original graph
// CRITICAL function, modifies graph
void Input :: unperturb_CRITICAL(map<pair<string, string>, double> epr)
{
	// For each edge in the map reset the edge probability in the graph to original value (as given by the map)
	map<pair<string, string>, double>::const_iterator epr_itr = epr.begin();
	set<string> vset;	// to keep track of vertices being modified
	for(; epr_itr!=epr.end(); epr_itr++)
	{
		vset.insert(epr_itr->first.first);
		vset.insert(epr_itr->first.second);

		graph[epr_itr->first.first]->getNeighbour(epr_itr->first.second)->setProbability_CRITICAL(epr_itr->second);
	}

	// Re-compute number of unique labels
	unsigned card = getUniqLabsSize();
	
	// Re-compute vertex degrees, (non)existence probabilities of neighbour labels 
	// and expected symbol occurrence probability
 	set<string>::iterator vset_itr = vset.begin();
	for(; vset_itr!=vset.end(); vset_itr++)
    {
		graph[*vset_itr]->computeDegree();
		graph[*vset_itr]->compute_pr_lx();
		graph[*vset_itr]->compute_symOccPr(card);
	}
}


// Prints the graph characteristics
void Input :: printGraph(void) const
{
	cout<<"\nThe input graph contains the following vertices:"<<endl;

	unordered_map<string, Vertex*>::const_iterator it = graph.begin();
	for(; it!=graph.end(); it++)
		(it->second)->print();

	cout<<endl<<"The vertices associated with the labels are:"<<endl;
	unordered_map<string, vector<Vertex*> >::const_iterator it1 = vertexLabel.begin();

	for(; it1!=vertexLabel.end(); it1++)
	{
		cout<<it1->first<<": ";
		const vector<Vertex*> vert = it1->second;

		for(unsigned j=0; j<vert.size(); j++)
			cout<<(vert[j])->getID()<<" ";

		cout<<endl;
	}
}


// Returns the number of vertices in the graph
unsigned Input :: getGraphSize(void) const
{
	return graph.size();
}


// Add label to the unique label set
void Input :: addLabel(string label)
{
	uniq_labs.insert(label);
}


// Get the number of uniqLabels present in the graph, used to compute the chi square values of the vertices
unsigned Input :: getUniqLabsSize(void) const
{
	return uniq_labs.size();
}


// Returns the vertex IDs of the top-k matching Subgraphs (to the provided query) found in the input graph
vector<vector<Vertex*> > Input :: getSubGraphs(const Query& qry)
{
	#ifdef DEBUG
			cout<<"\nComputing chi-square values for vertices corresponding to best ip_vertex-q_vertex pair...\n";
	#endif // DEBUG

	vector<vector<Vertex*> > subgraph;  // Stores the top-k approximate matching subgraph
	bool match;  // Checks if the input vertex label is present in the query graph or not
	set<string> uniq_labs;  // Number of uniq labels in the input graph

	// Compare the vertex label and its neighbours with that in the query graph
	// to obtain the symbols for the chi-square computation
	// Only labels from the query graph are selected for comparison
	const set<string> q_labels = qry.getLabels();  // Get the labels in the query graph
	set<string>::const_iterator itq = q_labels.begin();

	// For each of the labels in the query, evaluate the vertices in the input with the same label
	set<string> q_visited;		// Store visited query vertices
	for(; itq!=q_labels.end(); itq++)
	{
		string label = *itq;

		match = false;

		vector<Vertex*> cand;   // Store all the input vertices with the same label

		if(vertexLabel.find(label) != vertexLabel.end())
		{
			cand = vertexLabel[label];
			match = true;
		}

		if(!match)  // 	Modify here for structure matching here
			continue;

		// Get the label neighbourhood in the query
		// No vertex is visited yet, so q_visited empty
		// In this case all possible label neighbourhoods are 
		unordered_map<string, vector<string> > qlabel_neig = qry.getLabelNeighbours(q_visited, label);

        // For each candidate vertex
		for(unsigned i=0; i<cand.size(); i++)
		{
			Vertex *ver = cand[i];
			ver->computeChiSqValue(qlabel_neig, true);

			// Keep the top-k vertices with the maximum chi value as candidates
			if(heap.size() < (ORDER_CONSTANT*TOPK))
			{
				heap.push(ver);
			}
			else
			{
				const Vertex *v = heap.top();

				// If the minimum chiValue vertex has chi value less than candidate vertex ver, discard the min chiValue vertex
				if(v->getChiValue() < ver->getChiValue())
				{
					heap.pop();
					heap.push(ver);
				}
			}
		}
	}   // Populating min heap complete, heap now contains ORDER_CONSTANT*TOPK vertex pairs with maximum chi-sq values

	// Initialize answer variable subgraph with the candidates present in the heap
	// Each vertex pair is inserted as an independent candidate answer
	// Note that heap is a min-heap, so vertex pairs are inserted lowest chi-sq first
	// All at the beginning of subgraph vector, after the loop in vector varable subgraph
	// Chi-sq value of vertex pairs decreases as the index position increases
	// Thus, subgraph serves as Primary heap, during exploration and as answer variable at the end of it
	while(!heap.empty())
	{
		vector<Vertex*> cand;
		cand.push_back(heap.top());
		heap.pop();
		subgraph.insert(subgraph.begin(), cand);
	}

	#ifdef DMEASURE
		// For measuring heap sizes
		size_t primHeap_size = subgraph.size();		// Number of entries in Primary Heap
		cout<<"\nSize of Primary heap = "<<primHeap_size;
		vector<size_t> max_secondary_heapSize;
	#endif // DMEASURE

	set<Vertex*> duplicate; //Store visited vertex

	#ifdef DEBUG
		// Check vertex pairs
		if(subgraph.size()>0)
		{
			cout<<"\n\Printing each candidate ip_vertex-q_vertex pair: ";
			for(unsigned dbg_i=0; dbg_i<subgraph.size(); dbg_i++)
				cout<<"("<<subgraph[dbg_i][0]->getID()<<", "<<subgraph[dbg_i][0]->get_qVertexID()<<") ";
			cout<<endl;
		}
		else
		{
			cout<<"\n\nNO CANDIDATE VERTEX PAIR FOUND!!!";
		}
	#endif // DEBUG

	// Run the greedy approach using the candidate vertices obtained
	unsigned max_vertex = qry.getGraphSize();
	for(unsigned i=0; i<subgraph.size(); i++)   // For each candidate subgraph explore (grow) subgraph
	{
		// Max-heap for exploring subgraph, stores neighbour ranked on (chi square*neighbour probability)
		priority_queue<Vertex*, vector<Vertex*>, Compare_max> secondary_heap;
		secondary_heap.push(subgraph[i][0]); // Push subgraph vertex in primary heap

		(subgraph[i]).clear();

		// To ensure that for each iteration all query vertices are mapped and injectively
		q_visited.clear();

		#ifdef DEBUG
			cout<<"\n\nExploring subgraph "<<i<<":";
		#endif // DEBUG

		#ifdef DMEASURE
			size_t max_heapSize = 0;
		#endif // DMEASURE

		while( ((subgraph[i]).size() < max_vertex) && (!secondary_heap.empty()) )
		{
			#ifdef DEBUG
				cout<<"\nPrimary Heap: ";
				priority_queue<Vertex*, vector<Vertex*>, Compare_max> copy_ph = secondary_heap;
				while(!copy_ph.empty()) { cout<<copy_ph.top()->getID()<<" "; copy_ph.pop();}
			#endif // DEBUG
			#ifdef DMEASURE
				size_t cur_heapSize = secondary_heap.size();		// Maximum number of elements stored in the Secondary Heap
				if(max_heapSize < cur_heapSize)
				{
					max_heapSize = cur_heapSize;
				}
			#endif //DMEASURE

			// Get the vertex with the maximum chi-sq value and insert in heap
			Vertex *cand = secondary_heap.top();
			secondary_heap.pop();

			if(duplicate.find(cand) != duplicate.end() )
			{
				// Vertex *cand is already visited in some iteration, so discard
				continue;
			}

			if(q_visited.empty() && cand->get_qVertexID_dirty_bit())
			{
				// Dirty bit is set implies that the candidate vertex was explored as a Secondary heap candidate
				// And the chi-values or labels may have now changed from the initial computation done as a Primary heap candidate
				// Condition q_visited.empty() ensures that the condition is being checked for only the Primary heap candidate pairs
				unordered_map<string, vector<string> > qlabel_neig = qry.getLabelNeighbours(q_visited, cand->getLabel());
				cand->computeChiSqValue(qlabel_neig, true);
				cand->set_qVertexID_dirty_bit(0);

			}

			// To avoid visiting a pair with visited query vertex
			// Situation will arise if (v_i, q_i) was added to Secondary heap after (v_j, q_i)
			// (this addition could be across iterations of exploration)
			// And (v_i,q_i) has higher chi square value than (v_j, q_i)
			// This leads to q_i being marked visited and chi-sq for v_j is now rquired to be computed again
			if(q_visited.find(cand->get_qVertexID())!=q_visited.end())
			{
				string neig_label = cand->getLabel();

				// Get possible neighbor query labels for the neighbor vertex label
				unordered_map<string, vector<string> > q_labs = qry.getLabelNeighbours(q_visited, neig_label);

				bool match = true;
				if(q_labs.size() == 0)  // If no query neighborhood is found for the label of neighbor vertex (of cand)
				{
					match = false;
					continue;           // Discard and move on
				}

				cand->computeChiSqValue(q_labs, match);
			}

			// candidate vertex explored, mark visited
			(subgraph[i]).push_back(cand);
			duplicate.insert(cand);

			q_visited.insert(cand->get_qVertexID());

			// Find neighbours of the candidate vertex
			// Add to Secondary heap based on intersection with chi labels
			set<string> neighID = cand->getNeighbourIDs();
			vector<string> chiLab = cand->getChiLabels();

			set<string>::iterator it = neighID.begin();
			
			for(; it!=neighID.end(); it++)  // For each neighbor vertex
			{
				Vertex *view = graph[*it];  // Neighbor vertex of cand

				if(duplicate.find(view) != duplicate.end())		// If view is visited then discard
				{
					continue;
				}

				// For unvisited vertex or vertices mapped to a visited query vertex
				// (Re)compute chiValue and find the best query neighbourhood based on the labels
				if(view->getChiValue() == 0.0 || q_visited.find(view->get_qVertexID())!=q_visited.end())
				{
					if(view->getChiValue() != 0.0)
					{
						// The chi square value for vertex view is being recomputed
						// Set the dirty bit, indicating there has been a change in the uery vertex mapping
						// Useful if view is a Primary heap candidate and is not visited/ chosen
						// Until explored as a Primary heap vertex pair
						view->set_qVertexID_dirty_bit(1);
					}

					string neig_label = view->getLabel();

					// Get possible query label neighborhoods for chi-sq computation of view
					unordered_map<string, vector<string> > q_labs = qry.getLabelNeighbours(q_visited, neig_label);

					bool match = true;
					if(q_labs.size() == 0)  // If no query neighborhood is found for view (neighbor vertex of cand)
					{
						match = false;
						continue;           // Discard and move on
					}

					view->computeChiSqValue(q_labs, match);
				}

				// Check if the label of view is present in the chiLabels of the cand
				if(find(chiLab.begin(), chiLab.end(), view->getLabel())!=chiLab.end())
				{
					view->setParent(cand);
					secondary_heap.push(view);		// Insertion in heap based on chisq value and the edge probability of edge connecting "parent"
					
					// Note: Since the edge probability is required only once during insertion
					// Even if the parent is over-written in some iteration it doesn't affect the computation
				}

			}	// Exploring neighbors of cand
		}	// While loop for exploring a subgraph ends

		// The Primary heap vertex was already visited in a previous iteration and hence not explored
		if(subgraph[i].empty())
		{
            subgraph.erase(subgraph.begin()+i);
            i--;
		}

		#ifdef DEBUG
			// Check the subgraph computed
			if(!subgraph[i].empty())
			{
				cout<<endl;
				for(unsigned j=0; j<subgraph[i].size(); j++)
					cout<<subgraph[i][j]->getID()<<" ("<<subgraph[i][j]->get_qVertexID()<<", "<<subgraph[i][j]->getChiValue()<<"), ";

			}
		#endif //DEBUG
		#ifdef DMEASURE
			if(!subgraph[i].empty())
			{
				max_secondary_heapSize.push_back(max_heapSize);
			}
		#endif // DMEASURE
	}

	#ifdef DMEASURE
		if(subgraph.size()>0)
		{
			size_t max_secHS=0;
			float avg_secHS=0;
			for(unsigned int dm_i = 0; dm_i<max_secondary_heapSize.size(); dm_i++)
			{
				//cout<<"\nheap "<<max_secondary_heapSize[dm_i];
				if(max_secHS < max_secondary_heapSize[dm_i])
				{
					max_secHS = max_secondary_heapSize[dm_i];
				}
				avg_secHS += max_secondary_heapSize[dm_i];
			}
			avg_secHS /= max_secondary_heapSize.size();
			cout<<"\nMax size of secondary heap was = "<<max_secHS;
			cout<<"\nAverage size of secondary heap was = "<<avg_secHS;
			cout<<"\nSize of Vertex* = "<<sizeof(Vertex*);
		}
		else
		{
			cout<<"\nNO MATCHING SUBGRAPH FOUND!!!";
		}
	#endif // DMEASURE


	return subgraph;
}

