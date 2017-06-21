/*
 * Graph and its methods are defined in "utils.h" and implemented in "utils.cpp"
 *
 */
#include<iostream>
#include<fstream>
#include <vector>
#include <algorithm>
#include<list>
#include<iomanip>
#include <stdio.h>
#include<stdlib.h>
#include <mpfr.h>
#include "utils.h"
#include "genetic.h"
using namespace std;

#define THRESHOLD 			1e-5
#define INPUT_SIZE  		10
#define	INPUT_GRAPH			"graph_100.txt"
#define TRAINING_SET_SIZE	1000
#define PERCENTILE_GUIDED	true

int main() {
	/*
	 * Graphs input:
	 * graph_0.txt { 3.4325, 2.5435, 3.98876, 2.12345 }
	 * graph_1.txt { 3.4325, 2.5435, 3.98876, 2.12345, 3.98876, 2.12345 }
	 * graph_2.txt { 3.4325, 2.5435, 3.98876, 2.12345, 3.98876, 2.12345,
	 * 3.4325, 2.5435, 3.98876, 2.12345 }
	 * double input[INPUT_SIZE] = { 0.4325, 0.5435, 0.98876, 0.12345, 0.98876,
	 * 0.12345, 0.12345, 0.4325, 0.98876, 0.12345 };
	 *  */
	/*Precision p(INPUT_GRAPH); //Read graph
	p.set_error_threshold(THRESHOLD); //Set Error threshold
	double input[INPUT_SIZE] = { 12.32, 14.18, 16.21, 2.9214, 1.05, 4.02 };
	p.set_input(input, INPUT_SIZE);
	p.set_max_prec_values();
	int precisions[12] = { 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 };
	p.get_graph().set_precisions(precisions);
	p.set_input(input,INPUT_SIZE);
	cout << p.compute_error(precisions);
	p.get_graph().print_results();*/
	 Precision p(INPUT_GRAPH); //Read graph
	 p.set_error_threshold(THRESHOLD); //Set Error threshold
	 int counter = 0;
	 Statistics avg(p, TRAINING_SET_SIZE, PERCENTILE_GUIDED);
	 srand(time(0));
	 for (int i = 0; i < INPUT_SIZE; i++) {
	 avg.worst_input[i] = (double) (rand() / (1e+07));
	 }
	 do {
	 cout << endl
	 << "##########################################################################################"
	 << endl;
	 cout << "Iteration: " << counter << endl;
	 counter += 1;
	 cout << endl << "--->Iterative search";
	 avg.iterative_search();
	 avg.update_p_refined();
	 avg.find_worst_input();
	 } while (!avg.check_end());
	 cout << endl << endl << "Finished with " << counter << " iterations"
	 << endl;
	 /*
	 //Try Genetic
	 double input[INPUT_SIZE] = { 0.4325, 0.5435, 0.98876, 0.12345, 0.98876,
	 0.12345, 0.12345, 0.4325, 0.98876, 0.12345 };
	 p.set_input(input, INPUT_SIZE);
	 p.set_max_prec_values();
	 genetic(p);
	 p.iterative_search(input, INPUT_SIZE);
	 //p.iterative_search(input,INPUT_SIZE);
	 * */
	return 0;
}
