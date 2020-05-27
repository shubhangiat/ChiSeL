/*** The main file to reads the input and query graphs and subsequently obtains the
 *** approximate matching sub-graphs in the input graph to the query graph provided.
 *** Multiple query graphs can be given. Searching for multiple query graphs on the
 *** same  input graph  in one  execution  doesn't affect the  results, as for each
 *** query graph the chi square values are recomputed irrespective.
 ***/


#include "query_graph.h"
#include "input_graph.h"
#include "const.h"
#include <iostream>
#include <algorithm>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <unordered_map>


bool desc_sort(std::pair<vector<Vertex*>, double> a, std::pair<vector<Vertex*>, double> b)
{
	return (a.second > b.second);
}


bool vert_sort(Vertex *a, Vertex *b)
{
	return (a->getID() < b->getID());
}


int main(int argc, char **argv)
{
	clock_t begin, end;

	if(argc!=4)
	{
		cout<<"USAGE: ./subgraph <input graph vertex label file> <input graph edge file> <query arg file>";
        cout<<"\n"<<"Please go through the readMe file\n";
        return 0;
	}

	// Create the input graph
	string i_label_file = argv[1];
	string i_edge_file = argv[2];

	cout<<"Reading Input Graph files and Indexing Input Structure...\n";
	fflush(stdout);
	begin = clock();
	Input inp(i_label_file, i_edge_file);
	end = clock();
	cout<<"Read in "<<((end-begin)*1.0)/CLOCKS_PER_SEC<<" sec."<<endl;

#ifdef DEBUG
	// Check the input graph obtained
	inp.printGraph();
#endif // DEBUG

	// Create the query graph
	ifstream qfile;
	qfile.open(argv[3]);

	string f1, f2;
	qfile>>f1>>f2;	// Reading path of query vertex and graph edge file
	while(!qfile.eof())
	{
		cout<<"Query files:\t"<<f1<<"\t"<<f2<<endl;
		string q_label_file =f1;	//= argv[1];
		string q_edge_file = f2;	//argv[2];

		cout<<"Reading Query Graph files and Indexing Query Structure...\n";
		begin = clock();
		Query qry(q_label_file, q_edge_file);
		end = clock();
		cout<<"Done in "<<((end-begin)*1.0)/CLOCKS_PER_SEC<<" sec."<<endl;
		#ifdef DEBUG
				// Check the query graph obtained
				qry.printGraph();
		#endif // DEBUG

		#ifdef PERTURB_CRIT
				// perturb the input graph
				map<pair<string, string>, double> orig_epr = inp.perturb_CRITICAL(qry);
				#ifdef DEBUG
					// Check the input graph obtained
					inp.printGraph();
				#endif // DEBUG
		#endif	// PERTURB_CRIT

		// Compute the chi square values for all the vertices in the input graph. Then
		// find the approximate sub-graph(s) in the input graph matching the query graph
		// Also find the amount of time spent to find the match
		cout<<endl<<"Computing the Approximate Matching Subgraph based on Statistical Significance. Please wait...\n";
		fflush(stdout);
		begin = clock();
		vector<vector<Vertex*> > subgraphs = inp.getSubGraphs(qry);
		end = clock();
		cout<<"\nDone in = "<<((end-begin)*1.0)/CLOCKS_PER_SEC<<" sec."<<endl;

		#ifdef DEBUG 
			// Check the query glraph post processing for any accidental changes
			cout<<endl<<"The Query graph with Vertex IDs (and labels) and edges are: "<<endl;
			qry.printGraph();
		#endif // DEBUG

		vector<std::pair<vector<Vertex*>, long double> > topk;  // Stores the sorted final subgraphs

		cout<<endl<<"The top-"<<TOPK<<" subgraphs (with vertex IDs, query Vertex IDs, and Chi-square value) found are:"<<endl;
		for(unsigned i=0; i<subgraphs.size(); i++)		// Calculate total chi-square value of the answer subgraphs
		{
			sort(subgraphs[i].begin(), subgraphs[i].end(), vert_sort);

			double total_chi = 0.0;
			for(unsigned j=0; j<subgraphs[i].size(); j++)
				total_chi += (subgraphs[i][j])->getChiValue();

			topk.push_back(std::make_pair(subgraphs[i], total_chi));
		}

		cout<<endl;

		subgraphs.clear();

		sort(topk.begin(), topk.end(), desc_sort);

		// Print the answer graphs sorted on total chi-sq
		unsigned size = (topk.size() > TOPK)? TOPK : topk.size();
		for(unsigned i=0; i<size; i++)
		{
			for(unsigned j=0; j<((topk[i]).first).size(); j++)
				cout<<((topk[i]).first)[j]->getID()<<" ("<<((topk[i]).first)[j]->get_qVertexID()<<"), ";

			//cout<<" ==> "<<topk[i].second<<endl;
			printf(" ==> %Lf\n", topk[i].second);
		}
		#ifdef PERTURB_CRIT
				// unperturb the input graph
				inp.unperturb_CRITICAL(orig_epr);
		#endif	// PERTURB_CRIT
				cout<<"Input graph:"<<endl;

#ifdef DEBUG
	// Check the input graph obtained
	inp.printGraph();
#endif // DEBUG

		cout<<endl;
		qfile>>f1>>f2;	//	Reading next query vertex and graph edge file
	}

	qfile.close();

	return 0;
}
