How to run the binary file?
===========================

Compile all the header files and dependent cpp files through make command and then run the binary file with appropriate arguments.

make
./subgraph \<input graph vertex-label file\> \<input graph edge file\> \<list of query graph files\>


Parameters:
===========

The k for top-k matching subgraphs can be set in the const.h header file.


Graph files format:
===================

Vertex file: vid label
------------

e.g.

v1 l1
v2 l2
v3 l1


Edge file: vid1 vid2 edge_probability
----------

e.g.

v1 v2 pr12
v2 v3 pr23
v3 v1 pr31


This format is followed for both input target graph.
For the query graph the edge file has only two values, the vertex ids of an edge.



Format of 3rd argument (list of query graph files):
===================================================

<query graph vertex-label file 1> <query graph edge file 1>
<query graph vertex-label file 2> <query graph edge file 2>
