/*** Implements the vertex structure for the input graph (in vertex.h) and
 *** computes the chi-square value of the vertices based on the occurrence
 *** counts of the symbols
***/


#include "vertex.h"
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>

// C++ template to print vector container elements
template <typename T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
    os << "[";
    for (unsigned i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i != v.size() - 1)
            os << ", ";
    }
    os << "]";
    return os;
}

// Creates a vertex with the characteristics
Vertex :: Vertex(const string id, const string lab) : Node(id, lab)
{
	chiValue = 0.0L;
	chiLength = 0;
	qVertexID_dirty_bit = 0;

	neighbourLabels.clear();
	chiLabels.clear();

	degree = -1;

	// initializes the occurrence counts of the symbols
	for(int i=0; i<NUM_OF_SYMBOLS; i++)
	{
		memset( &(symbolOccurrence[i]), 0, sizeof(unsigned long) );
	}
}


// Deallocates space
Vertex :: ~Vertex()
{
	neighbourLabels.clear();
	chiLabels.clear();
}


// Add neighbor labels with associated probabilities
void Vertex :: addNeighbourLabels(Edge *e)
{
	neighbourLabels.push_back(e);
}


// Returns the chi-square value associated with the vertex
long double Vertex :: getChiValue(void) const
{
	return chiValue;
}


// Returns the chiLabels associated with this vertex
const vector<string>& Vertex :: getChiLabels(void) const
{
	return chiLabels;
}


//Returns the value of qVertexID_dirty_bit
const bool Vertex :: get_qVertexID_dirty_bit(void) const
{
	return qVertexID_dirty_bit;
}


//Sets the value of qVertexID_dirty_bit
void Vertex :: set_qVertexID_dirty_bit(bool b)
{
	qVertexID_dirty_bit = b;
}


//	Returns the query vertex the vertex is mapped to
const string Vertex :: get_qVertexID(void) const
{
	return qVertexID;
}


// Returns the count of occurrences of the symbols s0, s1, s2
const double* Vertex :: getSymbolOcc(void) const
{
	return symbolOccurrence;
}


// Returns the list of neighbour labels
const vector<string>& Vertex :: getNeighbourLabels(void) const
{
    vector<string> nLabels;
    vector<Edge*>::const_iterator it = neighbourLabels.begin();

    for(; it!=neighbourLabels.end(); it++)
        nLabels.push_back((*it)->getLabel());

	return nLabels;
}


// Return edge with passed neighbourID
Edge* Vertex :: getNeighbour(const string vID) const
{
	vector<Edge*>::const_iterator it = neighbourLabels.begin();
	Edge *e = NULL;

    for(; it!=neighbourLabels.end(); it++)
    {
		if((*it)->getID().compare(vID) == 0)
		{
			e = *it;
			break;
		}
    }

    return e;
}


// Get parent with parent vertexID and associated edge probability
const pair<string, double> Vertex :: getParent(void) const
{
	return parent;
}


// Set parent with parent vertexID and associated edge probability
void Vertex :: setParent(const Vertex * p)
{
	parent.first = p->getID();
	parent.second = getNeighbour(parent.first)->getProbability();
}


// Prints the vertex characteristics
void Vertex :: print(void) const
{
	cout<<"ID: "<<ID<<"\tLabel: "<<label<<"\tNeighbourIDs: ";

	set<string>::iterator it = neighbourIDs.begin();
	for(; it!=neighbourIDs.end(); it++)
		cout<<*it<<" ";

	cout<<"\tDegree: "<<degree;

	cout<<"\t(NeighbourLabels, probability): ";
	for(unsigned i=0; i<neighbourLabels.size(); i++)
		cout<<"("<<neighbourLabels[i]->getLabel()<<","<<neighbourLabels[i]->getProbability()<<") ";

	cout<<"\nP(l_x does not exist): [";
	unordered_map<string, vector<double> >::const_iterator iter = pr_lx.begin();
	for( ; iter != pr_lx.end(); iter++)
	{
		cout<<"("<<iter->first<<","<<iter->second[0]<<"),";
	}
	cout<<"]";

	cout<<"\tP(l_x exactly 1 exists): [";
	for( iter = pr_lx.begin(); iter != pr_lx.end(); iter++)
	{
		cout<<"("<<iter->first<<","<<iter->second[1]<<"),";
	}
	cout<<"]";

	cout<<"\tP(l_x >=1 exists): [";
	for( iter = pr_lx.begin(); iter != pr_lx.end(); iter++)
	{
		cout<<"("<<iter->first<<","<<iter->second[2]<<"),";
	}
	cout<<"]\n";

	if(chiLength!=0)
	{
		cout<<"\tMapped to query vertex: "<<qVertexID;

		cout<<"\tChi-square string length: "<<chiLength;

		cout<<"\tChi-Sq Value: "<<chiValue<<endl;

//		cout<<"\tChi labels: ";
//		for(unsigned i=0; i<chiLength; i++)
//			cout<<chiLabels[i]<<" ";
	}

	cout<<endl;
}


// Computes the degree of the vertex given probabilistic neighbours
void Vertex :: computeDegree(void)
{
    degree = 0;

    for(unsigned i=0; i<neighbourLabels.size(); i++)
    {
        degree+=(neighbourLabels[i]->getProbability());
    }

	#ifdef DEBUG
		cout<<"\nDegree computed for "<<this->ID<<": "<<degree<<endl;
	#endif // DEBUG
}


// Computes (non)existence probabilities of each unique neighbour label (zero, exactly one and at least one instance exists)
void Vertex :: compute_pr_lx(void)
{
	vector<Edge*> vlabel = neighbourLabels;
	unordered_map<string, vector<double> > label_pr;	// for each label, a vector of probability of existence for each instance of it

	// iterating over neighbours to get unique labels and keep track of probability of existence for each instance of it
	for(unsigned v1=0; v1<vlabel.size() && vlabel.size()>1; v1++)
	{
		string l = vlabel[v1]->getLabel();		// Label being explored
		if(label_pr.find(l) == label_pr.end())
		{
			// Initialise placeholder for key l in label_pr
			vector<double> v;
			v.push_back(vlabel[v1]->getProbability());
			label_pr[l] = v;
		}
		else
		{
			(label_pr[l].push_back(vlabel[v1]->getProbability()));
		}
	}

	// Computing (non)existence probabilities of each unique neighbour label (zero, exactly one and at least one instance exists)
	// Loop over unique neighbourhood labels
	unordered_map<string, vector<double> >::iterator iter = label_pr.begin();
	for( ; iter != label_pr.end(); iter++)
	{
		string l = iter->first;
		double pr_lx_ni=1;	// Probability that no instance of l_x exists
		double pr_lx_e1=0;	// Probability that exactly one instance of l_x exists
		double pr_lx_al1=1;	// Probability that at least one instance of l_x exists

		vector<double> lx_inst_pr = iter->second;

		// Calculate probability that no instance of l_x exists
		for(unsigned int i=0; i<lx_inst_pr.size(); i++)
		{
			pr_lx_ni *= (1 - lx_inst_pr[i]);
		}

		// Calculate probability that exactly one instance of l_x exists
		if(pr_lx_ni!=0)	// To check if an instance exists with probability 1 (avoid division with 0 and miscalculation)
		{
			for(unsigned int i=0; i<lx_inst_pr.size(); i++)
			{
				pr_lx_e1 += (pr_lx_ni*lx_inst_pr[i]/(1 - lx_inst_pr[i]));
			}
		}
		else
		{
			// Calculate from scratch
			for(unsigned int i=0; i<lx_inst_pr.size(); i++)
			{
				double pr_lx_e1_subterm = 1;

				for(unsigned int j=0; j<lx_inst_pr.size(); j++)
				{
					if(j!=i)
						pr_lx_e1_subterm *= (1 - lx_inst_pr[j]);
				}

				pr_lx_e1 += lx_inst_pr[i] * pr_lx_e1_subterm;
			}
		}

		// Calculate probability that at least one instance of l_x exists
		pr_lx_al1 = (1 - pr_lx_ni);

		vector<double> v;
		v.push_back(pr_lx_ni);
		v.push_back(pr_lx_e1);
		v.push_back(pr_lx_al1);

		pr_lx[l] = v;
	} // Loop over unique neighbourhood labels ends
}


// Computes expected distribution of symbols s0, s1 and s2 (in that order)
void Vertex :: compute_symOccPr(unsigned card)
{
	symbolOccProbability[0] = pow( pow( ( 1 - (1.0/card) ), degree), 2);
	symbolOccProbability[1] = 2 * (1 - pow( 1-(1.0/card), degree)) * pow( 1-(1.0/card), degree);
	symbolOccProbability[2] = pow(1 - pow( ( 1 - (1.0/card) ), degree), 2);
}


// Computes the chi-square value of the vertex using the symbolOccurrence[]
void Vertex :: computeChiSqValue(const unordered_map<string, vector<string> >& qry_label, bool match)
{
    string best_q;
    double maxChi = -1;
	#ifdef DEBUG
		double best_expected[3];
		double best_symOcc[3];
	#endif // DEBUG

	// Compute degree for vertex, if somehow missed
	if(degree<0)
		computeDegree();
	
	// For each query neighborhood (qry_label contains label neighbourhood (value) of the query vertex (key))
    unordered_map<string, vector<string> >::const_iterator it = qry_label.begin();

    for(; it!=qry_label.end(); it++)
	{
		vector<string> qlabel = it->second;

		//reinitialise symbol occurrences
		for(unsigned j=0; j<NUM_OF_SYMBOLS; j++)
		{
			symbolOccurrence[j] = 0;
		}

		// For query vertex having none or only one neighbor label
		while(qlabel.size()<2)
            qlabel.push_back(EMPTY_LABEL);

        // Comparing two label pairs from query neighborhood with two labels of vertex label
        for(unsigned q1=0; q1<qlabel.size(); q1++)
        {
            for(unsigned q2=q1+1; q2<qlabel.size(); q2++)
            {
				double symbol_count[NUM_OF_SYMBOLS]={0};

				if(pr_lx.find(qlabel[q1])==pr_lx.end() && pr_lx.find(qlabel[q2])==pr_lx.end())
				{
					//CASE IV: both of the query triplet labels don't match any neighbour labels of the vertex
					symbol_count[0] = 1;												// s0
					symbol_count[1] = 0;												// s1
					symbol_count[2] = 0;												// s2
				}
				else if(pr_lx.find(qlabel[q1])==pr_lx.end())
				{
					//CASE III-a: qlabel[q1] doesn't match any neighbour labels of the vertex, all values depend on qlabel[q2]
					symbol_count[0] = pr_lx[qlabel[q2]][0];								// s0
					symbol_count[1] = pr_lx[qlabel[q2]][2];								// s1
					symbol_count[2] = 0;												// s2
				}
				else if(pr_lx.find(qlabel[q2])==pr_lx.end())
				{
					//CASE III-b: qlabel[q2] doesn't match any neighbour labels of the vertex, all values depend on qlabel[q1]
					symbol_count[0] = pr_lx[qlabel[q1]][0];								// s0
					symbol_count[1] = pr_lx[qlabel[q1]][2];								// s1
					symbol_count[2] = 0;												// s2
				}
				else if(qlabel[q1]!=qlabel[q2])
				{
					//CASE I: qlabel[q1]!=qlabel[q2]
					symbol_count[0] = pr_lx[qlabel[q1]][0] * pr_lx[qlabel[q2]][0];		// s0
					symbol_count[2] = pr_lx[qlabel[q1]][2] * pr_lx[qlabel[q2]][2];		// s2
					symbol_count[1] = 1 - symbol_count[0] - symbol_count[2];			// s1
				}
				else
				{
					//CASE II: qlabel[q1]=qlabel[q2]
					symbol_count[0] = pr_lx[qlabel[q1]][0];								// s0
					symbol_count[1] = pr_lx[qlabel[q2]][1];								// s1
					symbol_count[2] = 1 - symbol_count[0] - symbol_count[1];			// s2
				}

				// Computing symbolOccurrence, summation of symbol_count over all possible query triplets
                for(unsigned j=0; j<NUM_OF_SYMBOLS; j++)
                {
					symbolOccurrence[j] += symbol_count[j];
                }

            }   // For loop for q2 ends
        }   // For loop for q1 ends

        // Computing chi-square value for given vertex and query vertex
        chiLength = (qlabel.size()*(qlabel.size()-1))/2;	// Number of possible query triplets
		chiValue = 0;	// Clear previous chiValue

		#ifdef DEBUG_CHI_COMPUTATION
			cout<<endl<<this->getID()<<", "<<this->getLabel()<<":\t"<<it->first;
			cout<<endl<<symbolOccurrence[0]<<", "<<symbolOccurrence[1]<<", "<<symbolOccurrence[2];
			cout<<endl<<symbolOccProbability[0]*chiLength<<", "<<symbolOccProbability[1]*chiLength<<", "<<symbolOccProbability[2]*chiLength;
		#endif // DEBUG_CHI_COMPUTATION

		for(unsigned j=0; j<NUM_OF_SYMBOLS; j++)
		{
			double expected = symbolOccProbability[j]*chiLength + LAPLACIAN_BIAS;
			if(expected == 0.0)
			{
				// cout<<"Expected symbol occurrence found to be 0 for vertex ID: "<<ID<<endl;
				chiValue += 0.0;
				break;
			}

			chiValue += ( pow(expected - symbolOccurrence[j], 2.0) / expected);
		}

		if(maxChi < chiValue)
		{
			maxChi = chiValue;
			best_q = it->first;
			#ifdef DEBUG
				for(int iter=0; iter<3; iter++)
				{
					best_symOcc[iter] = symbolOccurrence[iter];
					best_expected[iter] = symbolOccProbability[iter]*chiLength;
				}
			#endif // DEBUG
		}

	}   // For loop iterating over query neighborhoods ends

	// Restoring best chi value and corresponding query neighbor set
	chiValue = maxChi;
	qVertexID = best_q;
	chiLabels = qry_label.at(best_q);
	chiLength = qry_label.at(best_q).size()*(qry_label.at(best_q).size()-1)/2;	// Number of query triplets
	if(chiLength==0)
		chiLength = 1;		// When query vertex has only one neighbour
		
	#ifdef DEBUG
		cout<<"\nChi square was computed for ("<<this->ID<<","<<qVertexID<<") as: "<<chiValue
			<<"\tsymOcc: ["<<best_symOcc[0]<<", "<<best_symOcc[1]<<", "<<best_symOcc[2]<<"]"
			<<"\t\texpected: ["<<best_expected[0]<<", "<<best_expected[1]<<", "<<best_expected[2]<<"]";
	#endif // DEBUG
}
