/*
 * util.cpp
 *
 *  Created on: 30/mag/2017
 *      Author: antonio
 */
#include<iostream>
#include<fstream>
#include <vector>
#include <algorithm>
#include <list>
#include <stdio.h>
#include<stdlib.h>
#include <iomanip>
#include <mpfr.h>
#include <array>
#include <math.h>
#include "utils.h"
#include <time.h>
using namespace std;
const int upper_bit_bound = 53;
const int lower_bit_bound = 4;

void sleepcp(int milliseconds) // cross-platform sleep function
		{
	clock_t time_end;
	time_end = clock() + milliseconds * CLOCKS_PER_SEC / 1000;
	while (clock() < time_end) {
	}
}
void print_array(int *array, int N) {
	cout << endl;
	for (unsigned int n = 0; n < N; n++) {
		cout << array[n] << "-";
	}
	cout << endl;
}

void operation(string op, mpfr_t left_operand, mpfr_t right_operand,
		mpfr_t result) {
	//Bind string op to the correct mpfr_operation
	if (op == "*") {
		mpfr_mul(result, left_operand, right_operand, MPFR_RNDZ);
	} else if (op == "/") {
		if (mpfr_get_d(right_operand, MPFR_RNDZ) != 0) {
			mpfr_div(result, left_operand, right_operand, MPFR_RNDZ);
		} else {
			cout << "ERROR_Division_by_zero, return random value" << endl;
			mpfr_set_d(right_operand,17.0,MPFR_RNDZ);
			mpfr_div(result, left_operand, right_operand, MPFR_RNDZ);
			//throw exception();
		}
	} else if (op == "-") {
		mpfr_sub(result, left_operand, right_operand, MPFR_RNDZ);
	} else if (op == "+") {
		mpfr_add(result, left_operand, right_operand, MPFR_RNDZ);
	} else {
		cout << "Unknown operator, please add \"" + op << "\"" << endl;
		throw exception();
	}
	/* Debug code
	 cout << "Operation: ";
	 mpfr_printf("%.15Rf", left_operand);
	 cout << op;
	 mpfr_printf("%.15Rf", right_operand);
	 cout << "----";
	 cout << "result: ";
	 mpfr_printf("%.15Rf", result);
	 cout << endl;*/
}

/*
 * Graph is a class created to create and compute
 * a Data Flow Graph.
 * Each vertex contains a mpfr variable, its
 * precision value, and the initialization value.
 * The graph is directed and can be computed
 * with DFS with post order.
 * Each vertex which is not a terminal vertex is an
 * operator. The results of operations are saved in
 * operator vertexes.
 * The DFS is implemented to visit each vertex even if
 * the graph is not connected.
 */

Graph::Graph() {
}
Graph::Graph(int V) {
	this->V = V;
	adj = new list<int> [V];
	dep_path = new list<int> [V + 1];
	prop = new vertex[V];
}
Graph::Graph(string name) {
	vector<string> vertexs;
	vector<string> edges;
	ifstream dotfile(name.c_str());
	string line;
	bool first = true;
	if (dotfile.is_open()) {
		while (!dotfile.eof()) {
			//First line contains "Digraph{" then is ignored
			getline(dotfile, line);
			if (first) {
				first = false;
				continue;
			}
			//read vertexes and edges
			if (line[line.size() - 1] == ';') {
				vertexs.push_back(line);
			} else {
				edges.push_back(line);
			}
		}
		edges.pop_back(); //Last line is a }
		//Graph initialization
		V = vertexs.size();
		adj = new list<int> [V];
		dep_path = new list<int> [V + 1];
		prop = new vertex[V];
		/*Build the graph
		 * Add nodes
		 */
		for (size_t i = 0; i < V; i++) {
			line = vertexs[i];
			size_t first = line.find("\"") + 1;
			size_t second = line.find_last_of("\"");
			prop[i].valueRead = line.substr(first, second - first);
		}
		//Add edges
		size_t nEdges = edges.size();
		for (size_t i = 0; i < nEdges; i++) {
			line = edges[i];
			size_t index = line.find("->");
			int first = atoi(line.substr(0, index).c_str());
			int second = atoi(line.substr(index + 2, line.size() - 1).c_str());
			this->add_edge(first, second);
		}
		print_adj_lists();
		N = get_n_input();
		init_mpfr();		//Initialize mpfr variables
		find_dependency_paths();		//Set dependencies path
	} else {
		cout << "Error- Input graph file can't be opened" << endl;
	}
}
int Graph::get_V() {
	return V;
}
int Graph::get_precision(int v) {
	return mpfr_get_prec(prop[v].mpfr_var);
}
double Graph::get_value(int v) {
	return prop[v].value;
}
void Graph::print_adj_lists() {
	/*
	 * Print the adj lists in graphviz readable format
	 */
	cout << "digraph G {" << endl;
	for (int a = 0; a < V; a++) {
		cout << a << "[label=\"" << prop[a].valueRead << "\"];" << endl;
	}
	list<int>::iterator i;
	for (int a = 0; a < V; a++) {
		for (i = adj[a].begin(); i != adj[a].end(); ++i) {
			cout << a << "->" << *i << endl;
		}

	}
	cout << "}" << endl;
}
void Graph::print_results() {
	/*
	 * Print the adj lists in graphviz readable format
	 * Show the actual values of each vertex
	 */
	cout << "digraph G {" << endl;
	for (int a = 0; a < V; a++) {
		cout << a << "[label=\"";
		mpfr_printf("%.5Rf", prop[a].mpfr_var);
		cout << "\"];" << endl;
	}

	list<int>::iterator i;
	for (int a = 0; a < V; a++) {
		for (i = adj[a].begin(); i != adj[a].end(); ++i) {
			cout << a << "->" << *i << endl;
		}

	}
	cout << "}" << endl;
}
void Graph::add_edge(int v, int w) {
	// Add w to v’s list.
	adj[v].push_back(w);
}
void Graph::add_properties(string* properties) {
	/*
	 * Save properties readen from file
	 * in the Graph element
	 */
	for (int i = 0; i < V; i++) {
		prop[i].valueRead = properties[i];
	}
}
bool Graph::DFS_print(int v) {
	/*
	 * Print graph vertexes visited with
	 * DFS in post order
	 */
	if (adj[v].size() == 0) {
		return true;
	} else {
		list<int>::iterator i;
		for (i = adj[v].begin(); i != adj[v].end(); ++i) {
			if (DFS_print(*i)) {
				cout << *i;
			}
		}
		return true;
	}
}
int Graph::get_n_input() {
	/*
	 * Get the number of terminal vertexes
	 * in the graph. This number coincides
	 * with the number of inputs
	 */
	int N = 0;
	for (unsigned int i = 0; i < V; i++) {
		if (adj[i].size() == 0) {
			N++;
		}
	}
	return N;
}
void Graph::init_mpfr() {
	/*
	 * Init the mpfr variables contained in each vertex
	 * and set their precision to upper_bit_bound bits
	 */
	for (unsigned int i = 0; i < V; i++) {
		prop[i].precision = upper_bit_bound;
		mpfr_init2(prop[i].mpfr_var, prop[i].precision);
	}
}
void Graph::set_precision(int v, int precision) {
	//Set precision of a vertex
	prop[v].precision = precision;
	int p = precision;
	if (!(p >= 2
			&& p <= ((mpfr_prec_t) ((mpfr_uprec_t) (~(mpfr_uprec_t) 0) >> 1)))) {
		cout << "ERROR PRECISION: " << p << endl;
		cout << "Possible missing precision value..." << endl;
	}
	mpfr_set_prec(prop[v].mpfr_var, precision);
}
void Graph::set_value(int v, double value) {
	//Set value of a vertex
	prop[v].value = value;
	mpfr_set_d(prop[v].mpfr_var, value, MPFR_RNDZ);
}
void Graph::set_precisions(int* precisions) {
	/*
	 * Modify precision of mpfr variables according
	 * to precision array.
	 */
	for (unsigned int i = 0; i < V; i++) {
		set_precision(i, precisions[i]);
	}
}
void Graph::fill_dep_path(list<int> root_path, int root, int v) {
	root_path.push_back(v);
	list<int>::iterator i;
	//cout<<"v: "<<v<<endl;
	for (i = root_path.begin(); i != root_path.end(); i++) {
		if (find(dep_path[v].begin(), dep_path[v].end(), *i)
				== dep_path[v].end()) {
			/*cout<<"v: "<<v<<"  -->";
			cout<<"inserted: "<<*i<<endl;*/
			dep_path[v].push_back(*i);
		}
	}

	root = v;
	if (adj[v].size() != 0) {
		for (i = adj[v].begin(); i != adj[v].end(); ++i) {
			fill_dep_path(dep_path[v], root, *i);
		}
	}
}
void Graph::find_dependency_paths() {
	int roots[V];
	for (unsigned int v = 0; v < V; v++) {
		roots[v] = v;
	}
	list<int>::iterator i;
	list<int>::iterator i2;
	std::list<int> roots_list(roots, roots + V);
	for (unsigned int v = 0; v < V; v++) {
		for (i = adj[v].begin(); i != adj[v].end(); ++i) {
			roots_list.remove(*i);
		}
	}
	for (i = roots_list.begin(); i != roots_list.end(); ++i) {
		dep_path[*i].push_back(*i);
		for (i2 = adj[int(*i)].begin(); i2 != adj[int(*i)].end(); ++i2) {
			cout << *i<<" : "<<*i2 << endl;
			fill_dep_path(dep_path[*i], *i, *i2);
		}
	}
	cout << endl << "Dependency paths";
	for (unsigned int v = 0; v < V; v++) {
		for (i = dep_path[v].begin(); i != dep_path[v].end(); ++i) {
			cout << *i << "-";
		}
	}
	cout << endl << endl;
	//Create a path with all the vertexes to be used to avoid infinite loops
	for (unsigned int v = 0; v < V; v++) {
		dep_path[V].push_back(v);
	}
}
void Graph::increase_path_precisions(int vertex, int* actual_precisions) {
	/*Increase of one bit precisions
	 *of variables included in the dependency
	 *path of vertex
	 */
	list<int>::iterator i;
	for (i = dep_path[vertex].begin(); i != dep_path[vertex].end(); ++i) {
		if (actual_precisions[*i] < upper_bit_bound) {
			actual_precisions[*i] += 1;
		}
	}
}
void Graph::init_input(double* values) {
	/*
	 * Initialize input array and vertexes value
	 * Inputs must be given as first elements
	 * in the graph.
	 */
	input = new double[N];
	for (unsigned int i = 0; i < N; i++) {
		input[i] = values[i];
		set_value(i, values[i]);
	}
}
void Graph::set_input() {
	/*
	 * Set input values in vertexes
	 */
	for (unsigned int i = 0; i < N; i++) {
		set_value(i, input[i]);
	}
}
void Graph::DFS_compute_ret(int v, mpfr_t* final) {
	/*
	 * Compute the DFG starting from a vertex v
	 * Value of vertex v is saved in final
	 *
	 */
	if (adj[v].size() != 0 and !prop[v].visited) {
		list<int>::iterator i = adj[v].begin();
		list<int>::iterator left = i;
		++i;
		list<int>::iterator right = i;
		DFS_compute_ret(*left, &prop[*left].mpfr_var);
		DFS_compute_ret(*(right), &prop[*(right)].mpfr_var);
		i = adj[v].begin();
		operation(prop[v].valueRead, prop[*left].mpfr_var,
				prop[*(right)].mpfr_var, prop[v].mpfr_var);
	}
	mpfr_set(*final, prop[v].mpfr_var, MPFR_RNDZ);
}
void Graph::DFS_compute(int v) {
	/*
	 * Compute the DFG starting from a vertex v
	 */
	if (adj[v].size() != 0 and !prop[v].visited) {
		list<int>::iterator i = adj[v].begin();
		list<int>::iterator left = i;
		list<int>::iterator right = ++i;
		DFS_compute(*left);
		DFS_compute(*(right));
		i = adj[v].begin();
		operation(prop[v].valueRead, prop[*left].mpfr_var,
				prop[*(right)].mpfr_var, prop[v].mpfr_var);
		prop[v].value = mpfr_get_d(prop[v].mpfr_var, MPFR_RNDZ);
	}
}
void Graph::DFS_compute_all() {
	/*
	 * Visit all the graph without starting vertex
	 * For each vertex compute the subgraph
	 */
	for (unsigned int i = 0; i < V; i++) {
		prop[i].visited = false;
	}
	for (unsigned int i = 0; i < V; i++) {
		if (prop[i].visited == false) {
			DFS_compute(i);
			prop[i].visited = true;
		}
	}
}

/*
 * Precision class is used to implement the algorithm
 * and manage the elaboration of final precision array
 * and intermediary steps.
 */

Precision::Precision() {

}

Precision::Precision(string graph_file_path) {
	/*
	 * Precision constructor initialize graph variable
	 * with a new graph readen fom graph_file_path
	 * Initialize input array length to the number of
	 * inputs of the graph and the precisions vector to
	 * the number of vertexes of the graph
	 */
	graph = Graph(graph_file_path);
	N = graph.get_n_input();
	input = new double[N];
	max_prec_values = new double[graph.get_V()];
	initialize_precisions();
	initialize_MWL();
}
Graph Precision::get_graph() {
	/*
	 * Retun the graph object
	 */
	return graph;
}
int* Precision::get_MWL() {
	return MWL;
}
void Precision::initialize_precisions() {
	/*
	 * Initialize precision values
	 */
	int n_vertexes = graph.get_V();
	actual_precisions = new int[n_vertexes];
	init_precisions = new int[n_vertexes];
	old_precisions = new int[n_vertexes];
	fill(old_precisions, old_precisions + n_vertexes, 0);
	fill(actual_precisions, actual_precisions + n_vertexes, upper_bit_bound);
	fill(init_precisions, init_precisions + n_vertexes, upper_bit_bound);
}
void Precision::initialize_precisions(int* precisions) {
	/*
	 * Initialize precision values
	 */
	for (unsigned int v = 0; v < graph.get_V(); v++) {
		init_precisions[v] = precisions[v];
	}
}
void Precision::set_max_precisions() {
	initialize_precisions();
	graph.set_precisions(actual_precisions);
}
void Precision::initialize_MWL() {
	/*
	 * Initialize Minimum World Length array
	 */
	int n_vertexes = graph.get_V();
	MWL = new int[n_vertexes];
	fill(MWL, MWL + n_vertexes, lower_bit_bound);
}
void Precision::set_error_threshold(double value) {
	if (value >= 0 and value <= 1) {
		error_threshold = value;
	} else {
		cout << "Error, invalid threshold. Insert a value between 0 and 1"
				<< endl;
	}
}
void Precision::set_input(double* input, int N) {
	/*
	 * Initialize input values according to input array
	 */
	if (this->N == N) {
		for (unsigned int i = 0; i < N; i++) {
			this->input[i] = input[i];
		}
		graph.init_input(input);
	} else {
		cout << "Error, wrong number of inputs!";
		throw exception();
	}

}
void Precision::set_max_prec_values() {
	compute_graph(init_precisions);
	for (unsigned int v = 0; v < graph.get_V(); v++) {
		max_prec_values[v] = graph.get_value(v);
	}
}
int* Precision::get_actual_precisions() {
	//get actual precisions array
	return actual_precisions;
}

double Precision::get_threshold() {
	//get threshold value
	return error_threshold;
}
void Precision::compute_graph(int* precisions) {
	/*
	 * Compute the DFG with given precision and input
	 */
	graph.set_precisions(precisions);
	graph.set_input();
	graph.DFS_compute_all();
}
void Precision::compute_graph() {
	/*
	 * Compute DFG with its actual precision and inputs
	 */
	graph.set_input();
	graph.DFS_compute_all();
}
double Precision::compute_error(int* trial_prec) {
	int n_vertexes = graph.get_V();
	compute_graph(trial_prec);
	double signal_sqr = 0;
	double error_sqr = 0;
	for (unsigned int v = 0; v < n_vertexes; v++) {
		/*cout<<v<<"-->"<<(max_prec_values[v] - graph.get_value(v))
		 <<"    ("<<max_prec_values[v]<<","<<graph.get_value(v)<<")"<<endl;*/
		signal_sqr += pow(max_prec_values[v], 2);
		error_sqr += pow((max_prec_values[v] - graph.get_value(v)), 2);
	}
	double sqnr;
	if (error_sqr != 0.00) {
		sqnr = signal_sqr / error_sqr;
		if (sqnr != 0) {
			return 1.0 / sqnr;
		} else {
			return 0.00;
		}
	}
}
int Precision::binary_search_MWL(int* trial_prec, int v, int upper_bound,
		int lower_bound) {
	int trial = (upper_bound + lower_bound) / 2;
	double error;
	//cout << endl << "trial: " << trial;
	if (trial == lower_bound) {
		trial_prec[v] = trial;
		error = compute_error(trial_prec);
		trial_prec[v] = init_precisions[v];
		if (error < error_threshold) {
			/*cout << endl << "vertex " << v << ": prec: " << lower_bound
			 << " error: " << error;*/
			return lower_bound;
		} else {
			/*cout << endl << "vertex " << v << ": prec: " << upper_bound
			 << " error: " << error;*/
			return upper_bound;
		}
	}
	trial_prec[v] = trial;
	error = compute_error(trial_prec);
	//cout << " - error: " << error << endl;
	if (error < error_threshold) {
		binary_search_MWL(trial_prec, v, trial, lower_bound);
	} else {
		binary_search_MWL(trial_prec, v, upper_bound, trial);
	}
}
void Precision::compute_MWL() {
	/*
	 * Compute the Minimum World Length
	 */
	//cout << endl << "Computing MWL...";
	int n_vertexes = graph.get_V();
	int trial_prec[n_vertexes];
	int final_prec[n_vertexes];
	initialize_precisions(actual_precisions);
	for (unsigned int v = 0; v < n_vertexes; v++) {
		trial_prec[v] = init_precisions[v];
	}/*
	 cout << endl << "Binary search, minimum precision value for vertexes..."
	 << endl;*/
	for (unsigned int v = 0; v < n_vertexes; v++) {
		final_prec[v] = binary_search_MWL(trial_prec, v, init_precisions[v],
				MWL[v]);
	}
	for (unsigned int v = 0; v < n_vertexes; v++) {
		MWL[v] = final_prec[v];
	}
	cout << endl << "MWL: ";
	for (unsigned int v = 0; v < n_vertexes; v++) {
		cout << MWL[v] << "-";
		actual_precisions[v] = MWL[v];
	}
	cout << " error: " << compute_error(actual_precisions) << endl;
}
int Precision::bit_sum(int* precisions, int N) {
	int sum = 0;
	for (int n = 0; n < N; n++) {
		sum += precisions[n];
	}
	return sum;
}
int Precision::find_best_path() {
	/*
	 * Find the most significant group
	 */
	int n_vertexes = graph.get_V();
	double errors[n_vertexes];
	double costs[n_vertexes];
	double error_min = 1.0;
	double cost_min = 53 * n_vertexes;
	int best_path = 0;
	for (unsigned int v = 0; v < n_vertexes; v++) {
		graph.increase_path_precisions(v, actual_precisions);
		errors[v] = compute_error(actual_precisions);
		costs[v] = bit_sum(actual_precisions, n_vertexes);
		//cout << endl << v << ": " << errors[v];
		if (errors[v] < error_min) {
			error_min = errors[v];
			cost_min = costs[v];
			best_path = v;
		} /*else if (errors[v] == error_min && costs[v] < cost_min) {
		 error_min = errors[v];
		 cost_min=costs[v];
		 cout<<v<<" instead of "<<best_path<<endl;
		 best_path = v;
		 }*/
		for (unsigned int v = 0; v < n_vertexes; v++) {
			actual_precisions[v] = init_precisions[v];
		}
	}
	if (error_min == current_error) {		//Inifinite loop detected
		best_path = n_vertexes;
	}
	if (error_min == 1.0) {
		cout << "ERORR, no path found" << endl;
	}
	return best_path;
}
void Precision::best_prec_increase() {
//Increase of one bit the precision of the best path
	initialize_precisions(actual_precisions);
	int best_path = find_best_path();
	graph.increase_path_precisions(best_path, actual_precisions);//CHECK PASS ARRAY
	if (best_path == graph.get_V()) {
		cout << " infinite loop detected,";
	}
	cout << " increased path: " << best_path << " -->";
	/*cout << endl << "path: " << best_path << endl;
	 for (unsigned int v = 0; v < graph.get_V(); v++) {
	 cout << actual_precisions[v] << "--";
	 }
	 cout << endl << "ERROR: " << compute_error(actual_precisions) << endl;
	 */
}
bool Precision::check_convergence() {
	bool sum_convergence = false, MWL_convergence = true, infinite_loop = true;
	int sum_bit_actual = 0, sum_bit_old = 0;
	for (unsigned int v = 0; v < graph.get_V(); v++) {
		sum_bit_actual += actual_precisions[v];
		sum_bit_old += old_precisions[v];
		old_precisions[v] = actual_precisions[v];
		if (MWL[v] != actual_precisions[v]) {
			MWL_convergence = false;
		}
	}
	if (sum_bit_actual == sum_bit_old) {
		sum_convergence = true;
	}
	/*
	 * Strict rule
	 * or (unsigned int v = 0; v < graph.get_V(); v++) {
	 if (old_precisions[v] != actual_precisions[v]) {
	 infinite_loop = false;
	 old_precisions[v] = actual_precisions[v];
	 }

	 if (MWL[v] != actual_precisions[v]) {
	 convergence = false;
	 }
	 }
	 if (infinite_loop) {
	 cout << "Infinite_loop detected! Increasing each precision of 1 bit"
	 << endl;

	 /*for (unsigned int v=0;v<graph.get_V();v++){
	 MWL[v]=MWL[v]+1;
	 }
	 graph.increase_path_precisions(graph.get_V(), actual_precisions);
	 }*/
	return sum_convergence or MWL_convergence;
}
void Precision::grouped_upward() {
	//execute best_prec_increase untill the error is below the threshold
	cout << "Grouped Upward..." << endl;
	while (compute_error(actual_precisions) > error_threshold) {
		cout << "prec: ";
		for (int v = 0; v < graph.get_V(); v++) {
			cout << actual_precisions[v] << "-";
		}
		current_error = compute_error(actual_precisions);
		cout << " error: " << current_error;
		//sleepcp(500);
		best_prec_increase();
		cout << endl;
	}
	cout << "prec: ";
	for (int v = 0; v < graph.get_V(); v++) {
		cout << actual_precisions[v] << "-";
	}
	fill(init_precisions, init_precisions + graph.get_V(), upper_bit_bound);//CHANGE
	cout << " error: " << compute_error(actual_precisions);
	cout << endl;
}
void Precision::end() {
	cout << endl << "Iterative Search Finished" << endl;
	cout << "Precision: ";
	int sum = 0;
	for (unsigned int v = 0; v < graph.get_V(); v++) {
		cout << actual_precisions[v];
		sum += actual_precisions[v];
		if (v != graph.get_V() - 1) {
			cout << "-";
		}
	}
	cout << "  bit_sum: " << sum << endl;
	cout << "Error: " << compute_error(actual_precisions) << endl;
	//cout << "Threshold: " << error_threshold << endl;
}
void Precision::iterative_search(double *input, int N) {
	set_input(input, N);		//Set input variables
	compute_graph();		//Compute graph with initial precisions
	set_max_prec_values();		//Set maximum precision values
	/*cout << "Maximum precision values: " << endl;
	 get_graph().print_results();*/
	while (!check_convergence()) {		//Untill MWL==actual_precisions
		compute_MWL();		//ISOLATED DOWNWARD
		grouped_upward();		//GROUPED UPWARD
	}
	end();
}

Statistics::Statistics(Precision p, int M, bool percentile_guided = true) {
	this->p = p;
	this->M = M;
	this->V = p.get_graph().get_V();
	this->N = p.get_graph().get_n_input();
	this->percentile_guided = percentile_guided;
	this->worst_input = new double[N];
	this->sqnr_array = new double[M];
	this->current_sqnr = 0;
	p_refined = new int[V];
	for (unsigned int v = 0; v < V; v++) {
		p_refined[v] = 0;
	}
}
Precision Statistics::get_p() {
	return p;
}
void Statistics::find_worst_input() {
	double error_avg = 0;
	int worst_seed = 0;
	double max_error = 0.0;
	double errors[M];
	for (unsigned int m = 0; m < M; m++) {
		srand(m);
		double trial_input[N];
		for (unsigned int n = 0; n < N; n++) {
			trial_input[n] = (double) rand() / (1e+07);
		}
		p.set_input(trial_input, N);
		p.set_max_prec_values();
		double error = p.compute_error(p_refined);
		//cout<<m<<"--> "<<error<<endl;
		errors[m] = error;
		error_avg += error / M;
		sqnr_array[m] = error;
		if (error > max_error) {
			max_error = error;
			worst_seed = m;
		}
	}
	//cout << "worst: " << worst_seed << endl;
	if (percentile_guided == false) {
		target_sqnr = error_avg;
	} else {
		sort(errors, errors + M);
		int percentile_95 = M * 95 / 100;
		target_sqnr = errors[percentile_95];
	}
	srand(worst_seed);
	for (unsigned int n = 0; n < N; n++) {
		worst_input[n] = (double) rand() / (1e+07);
	}
	//p.set_input(worst_input, N);
}
void Statistics::update_p_refined() {
	for (unsigned int v = 0; v < V; v++) {
		p_refined[v] = max(p_refined[v], p.get_actual_precisions()[v]);
		p.get_actual_precisions()[v] = p_refined[v];
	}
}
bool Statistics::check_end() {
	cout << endl << "Actual conditions" << endl;
	if (percentile_guided) {
		cout << "95° percentile: ";
	} else {
		cout << "Average sqnr: ";
	}
	int sum = 0;
	cout << target_sqnr << endl;
	cout << "p_refined: ";
	for (int v = 0; v < V; v++) {
		sum += p_refined[v];
		cout << p_refined[v] << "-";
	}
	int max_bit_sum=V*53;
	cout <<endl<< " bit_sum: " << sum <<" out of "<<max_bit_sum<< endl;
	cout << "threshold= " << p.get_threshold() << endl;
	/*cout << "worst_input: ";
	 for (int v = 0; v < N; v++)
	 cout << worst_input[v] << ",";
	 cout << endl;*/
	if (target_sqnr <= p.get_threshold()) {
		cout<<"It is possible to represent this program with: "<<((double)sum/max_bit_sum)*100<<"%"<<" of originary memory"<<endl;
		return true;
	} else {
		if (target_sqnr == current_sqnr) {
			for (int v = 0; v < V; v++) {
				p_refined[v] += 1;
			}
		} else {
			current_sqnr = target_sqnr;
		}
		p.set_max_precisions();
		p.initialize_MWL();
		return false;
	}
}
void Statistics::iterative_search() {
	p.iterative_search(worst_input, N);
}
