# *ChiSeL*

## Introduction
*ChiSeL* is a novel algorithmic framework that uses the idea of
statistical significance for approximate subgraph matching on
uncertain graphs. For each candidate matching vertex in the target
graph that matches a query vertex, it computes its statistical
significance using the chi-squared statistic. The search algorithm
then proceeds in a greedy manner by exploring the vertex neighbors
having the largest chi-square score. The algorithm works well on
large real-life graphs with millions of vertices and billions of
edges and can compute answer seconds after reading the input graph.

Please cite our paper, if you use our source code.
* "ChiSeL: Graph Similarity Search using Chi-Squared Statistics in Large Probabilistic Graphs. VLDB'20"

## How to run the binary file?

Compile all the header files and dependent cpp files through make command and then run the binary file with appropriate arguments.

```
make  
./subgraph <input graph vertex-label file> <input graph edge file> <list of query graph files>
```

For instance,

```
./subgraph ip_v15_label.txt ip_v15_edge.txt qry_graphs.txt
```

## Parameters

The _k_ for top-_k_ matching subgraphs can be set in the _const.h_ header file.

### Flags

During compilation flags can be set in makefile.

* __-DDMEASURE__  
To measure heap sizes.

* __-DPERTURB_CRIT__  
To perturb the original graph at runtime. The probability of graph edges
that match with the query graph edges are set to 1 before the subgraph
search and are reset to old value after the computation of answers.  
__Note:__ All query graph vertex-ids must match the input graph vertex-ids for this.

## Argument and Graph files format:

### Vertex file 
File format: vid label

e.g.
```
v1 l1  
v2 l2  
v3 l1  
```

### Edge file

File format: vid1 vid2 edge_probability

e.g.
```
v1 v2 pr12
v2 v3 pr23
v3 v1 pr31
```

This format is followed for both input target graph.
For the query graph the edge file has only two values, the vertex ids of an edge.


### Format of third argument (list of query graph files)

```
<query graph vertex-label file 1> <query graph edge file 1>
<query graph vertex-label file 2> <query graph edge file 2>
```
