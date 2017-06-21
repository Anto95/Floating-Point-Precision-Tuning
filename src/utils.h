/*
 * utils.h
 *
 *  Created on: 30/mag/2017
 *      Author: antonio
 */
#include<iostream>
#include<fstream>
#include <string>
#include <vector>
#include <algorithm>
#include<list>
#include <stdio.h>
#include<stdlib.h>
#include <mpfr.h>
using namespace std;

#ifndef UTILS_H_
#define UTILS_H_

struct vertex {
	string valueRead;
	bool visited;
	int precision;
	double value;
	mpfr_t mpfr_var;
};
// Graph class represents a directed graph using adjacency list representation

class Graph {
	int V;    // No. of vertices
	int N;    //Number of input variables
	double* input; // Input of the graph
	list<int> *adj;    // Adjacency lists
	list<int> *dep_path; //Dependency paths lists
	vertex *prop;		//Property lists
public:
	Graph();
	Graph(int V);   // Constructor
	Graph(string name); //Constructor from file
	int get_V(); //Return number of vertexes
	int get_n_input(); //get number of inputs
	double get_value(int v); //get value of vertex v
	int get_precision(int v); //get precision value of vertex v
	void add_edge(int v, int w);   // function to add an edge to graph
	void add_properties(string* properties); //bind properties with vertexes
	void print_adj_lists(); //Print graph as an adjacency list
	bool DFS_print(int v);    // DFS traversal of the vertices reachable from v
	void print_results(); //Print graph results after compute
	void init_mpfr(); //Initialize vertexes' mpfr variables
	void init_input(double* values); //Initialize input array
	void set_precision(int v, int precision); //Set precision of a vertex
	void set_value(int v, double); //Set value of a vertex
	void set_precisions(int* precisions); //Initialize vertexes' mpfr variables' precision
	void set_input(); //set graph with input values
	void fill_dep_path(list<int> root_path, int root, int v); //Helper function for set_paths
	void find_dependency_paths(); //Initializes vertexes' dependencies paths
	void increase_path_precisions(int vertex, int* actual_precisions); //Increase the precision of the path of a given vertex
	void DFS_compute_ret(int v, mpfr_t* final); // DFS compute DFG starting from v and returns final
	void DFS_compute(int v); // DFS compute DFG starting from v
	void DFS_compute_all(); //Call DFS_compute for each unvisited vertex
	//int countVariables(int v);// return the number of variables which will be created to compute the DFG

};
class Precision {
private:
	Graph graph; //Graph analyzed
	double* max_prec_values; //Actual values contained in the graph
	int* old_precisions; //Old precision array (needed to detect infinite loops
	int* actual_precisions; //Actual precision array
	int* MWL; //Minimum World Length array
	int N; //Number of input variables
	double* input; // Input of the graph
	double error_threshold; //error_threshold for the precision tuning
	double current_error; //used to avoid infinite loops
public:
	int* init_precisions; //Initial precision array
	Precision(); //Default constructor
	Precision(string graph_file_path); //Constructor with name of graph_file
	Graph get_graph(); //Return graph object
	int* get_MWL(); //Return MWL
	int* get_actual_precisions(); //get actual precisions array
	double get_threshold(); //get threshold value
	void initialize_precisions(); //Set initial precisions
	void initialize_precisions(int* precisions); //Set initial precisions to precisions
	void set_max_precisions(); //Set each mpfr variable to a precision of 53 bit
	void initialize_MWL(); //Set initial MWL
	void set_error_threshold(double value); //Set error threshold to value
	void set_input(double* input, int N); //Set input of the graph
	void set_max_prec_values(); //Set actual values of the graph
	void compute_graph(int* precision); //Compute the DFG with given precision
	void compute_graph(); //Compute the graph as it is actually configured
	double compute_error(int* precision); //Compute error obtained computing the graph with a given precision
	int binary_search_MWL(int* trial_prec, int v, int upper_bound,
			int lower_bound); //Helper function for compute_MWL
	void compute_MWL(); //Compute Minimum World Length array
	int bit_sum(int* precisions, int N);//Sum number of bits of a precision result
	int find_best_path(); //Find the most significant path
	void increase_path_precision(int best_path); //Increase of one bit precisions of path path
	void best_prec_increase(); //Increase of one bit the precision of the best path
	bool check_convergence(); //Check the convergence between MWL and actual_precision
	void grouped_upward(); //execute best_prec_increase untill the error is below the threshold
	void iterative_search(double* input, int N); //Compute iterative search over an input
	void restore(); //Restore precision values to default
	void end(); //Print results

};

class Statistics {
private:
	Precision p;
	int V; //Number of vertexes
	int N; //Number of input vertexes
	int M; //Dimension of training set
	double* sqnr_array; //sqnr computed for each iteration of average refinement
	double target_sqnr;
	double current_sqnr;
	bool percentile_guided;//True by default, if false target_sqnr= average_sqnr
	int* p_refined;

public:
	double* worst_input;
	Statistics();
	Statistics(Precision p, int M, bool percentile_guided);
	void find_worst_input();
	void update_p_refined();
	Precision get_p();
	void iterative_search(); //Launch iterative search on worst_input
	bool check_end();

};

void operation(string op, mpfr_t left_operand, mpfr_t right_operand,
		mpfr_t* result);
void print_array(int* array, int N); //Print array elements
#endif /* UTILS_H_ */
