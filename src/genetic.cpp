//
//  code to illustrate the use of a genetic algorithm to solve the problem described
//  at
//
//  by Mat Buckland aka fup
//
//-----------------------------------------------------------------------------------------------
#include <string>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <time.h>
#include <math.h>
#include <cstring>
#include <algorithm>
#include "utils.h"
#include "genetic.h"

using std::string;
using namespace std;

#define CROSSOVER_RATE            0.7
#define MUTATION_RATE             0.001
#define POP_SIZE                  100           //must be an even number
#define CHROMO_LENGTH             126
// int from 4 to 64 is represented by 6 bits, condition: 4<x<53
#define GENE_LENGTH               6
#define MAX_ALLOWABLE_GENERATIONS   400

//returns a double between 0 & 1
#define RANDOM_NUM      ((double)rand()/(RAND_MAX))

//----------------------------------------------------------------------------------------
//
//  define a data structure which will define a chromosome
//

//-------------------------------main--------------------------------------------------
//
//-------------------------------------------------------------------------------------
int bit_cost = 0;
vector<int> total_costs;
double all[POP_SIZE];
int* genetic(Precision p) {
	//seed the random number generator
	srand((int) time(NULL));
	int total_counter = 0;
	//just loop endlessly until user gets bored :0)
	while (total_counter < 10) {
		total_counter++;
		//storage for our population of chromosomes.
		chromo_typ Population[POP_SIZE];

		//get a target number from the user. (no error checking)
		double Target;
		double threshold;
		threshold = 1e-15;
		Target = 1 / threshold;
		cout << endl << endl;

		//first create a random population, all with zero fitness.
		for (int i = 0; i < POP_SIZE; i++) {
			Population[i].bits = GetRandomBits(CHROMO_LENGTH);
			Population[i].fitness = 0.0;
		}
		int GenerationsRequiredToFindASolution = 0;

		//we will set this flag if a solution has been found
		bool bFound = false;

		//enter the main GA loop
		while (!bFound) {
			//this is used during roulette wheel sampling
			double TotalFitness = 0.0;

			// test and update the fitness of every chromosome in the
			// population
			for (int i = 0; i < POP_SIZE; i++) {
				Population[i].fitness = AssignFitness(Population[i].bits,
						threshold, p);
				TotalFitness += Population[i].fitness;

			}
			// check to see if we have found any solutions (fitness will be 999)
			for (int i = 0; i < POP_SIZE; i++) {
				all[i] = Population[i].fitness;
				if (Population[i].fitness == 999999999.0) {
					cout << "\nSolution found in "
							<< GenerationsRequiredToFindASolution
							<< " generations!" << endl << endl;
					;

					PrintChromo(Population[i].bits,p);

					bFound = true;

					break;
				}
			}
			// create a new population by selecting two parents at a time and creating offspring
			// by applying crossover and mutation. Do this until the desired number of offspring
			// have been created.

			//define some temporary storage for the new population we are about to create
			chromo_typ temp[POP_SIZE];

			int cPop = 0;

			//loop until we have created POP_SIZE new chromosomes
			while (cPop < POP_SIZE) {
				// we are going to create the new population by grabbing members of the old population
				// two at a time via roulette wheel selection.
				string offspring1 = RankSelection(Population);
				string offspring2 = RankSelection(Population);

				//add crossover dependent on the crossover rate

				Crossover(offspring1, offspring2);
				//now mutate dependent on the mutation rate
				Mutate(offspring1);
				Mutate(offspring2);

				//add these offspring to the new population. (assigning zero as their
				//fitness scores)
				temp[cPop++] = chromo_typ(offspring1, 0.0);
				temp[cPop++] = chromo_typ(offspring2, 0.0);

			}        //end loop
					 //copy temp population into main population array
			for (int i = 0; i < POP_SIZE; i++) {
				Population[i] = temp[i];
			}

			++GenerationsRequiredToFindASolution;

			// exit app if no solution found within the maximum allowable number
			// of generations
			if (GenerationsRequiredToFindASolution > MAX_ALLOWABLE_GENERATIONS) {
				sort(all, all + POP_SIZE);
				cout << "No solutions found this run! First ten: " << endl;
				for (int i = 0; i < 10; i++) {
					cout << all[POP_SIZE - i] << " - ";
				}
				bFound = true;
			}

		}

		cout << "\n\n\n";

	}          //end while
	sort(total_costs.begin(),total_costs.end());
	cout<<"Best: "<<total_costs[0]<<endl;
	int array[2] = { 1, 2 };
	return array;
}

//---------------------------------GetRandomBits-----------------------------------------
//
//  This function returns a string of random 1s and 0s of the desired length.
//
//-----------------------------------------------------------------------------------------
string GetRandomBits(int length) {
	string bits;
	int sum = 0;
	for (int i = 0; i < length; i++) {
		if (RANDOM_NUM > 0.5)

			bits += "1";

		else

			bits += "0";
	}
	return bits;
}

//---------------------------------BinToDec-----------------------------------------
//
//  converts a binary string into a decimal integer
//
//-----------------------------------------------------------------------------------
int BinToDec(string bits) {
	int val = 0;
	int value_to_add = 1;

	for (int i = bits.length(); i > 0; i--) {

		if (bits.at(i - 1) == '1')

			val += value_to_add;

		value_to_add *= 2;

	}          //next bit

	return val;
}
string DecToBin(int dec) {
	string r;
	while (dec != 0) {
		r = (dec % 2 == 0 ? "0" : "1") + r;
		dec /= 2;
	}
	return r;
}
string checkGene(string bits) {
	int this_gene = BinToDec(bits);
	if (this_gene > 53)
		this_gene = 53;
	else if (this_gene < 4)
		this_gene = 4;
	bits = DecToBin(this_gene);
	while (bits.size() < GENE_LENGTH) {
		bits.insert(0, "0");
	}
	return bits;
}

//---------------------------------ParseBits------------------------------------------
//
// Given a chromosome this function will step through the genes one at a time and insert
// the decimal values of each gene (which follow the operator -> number -> operator rule)
// into a buffer. Returns the number of elements in the buffer.
//------------------------------------------------------------------------------------
int ParseBits(string bits, int* buffer) {
	bit_cost = 0;
	//counter for buffer position
	int cBuff = 0;

	// step through bits a gene at a time until end and store decimal values
	// of valid operators and numbers. Don't forget we are looking for operator -
	// number - operator - number and so on... We ignore the unused genes 1111
	// and 1110

	//flag to determine if we are looking for an operator or a number
	bool bOperator = true;

	//storage for decimal value of currently tested gene
	int this_gene = 0;

	for (int i = 0; i < CHROMO_LENGTH; i += GENE_LENGTH)
	{
		//convert the current gene to decimal
		string checked = checkGene(bits.substr(i, GENE_LENGTH));
		bits.replace(i, GENE_LENGTH, checked);
		this_gene = BinToDec(checked);
		buffer[cBuff] = this_gene;
		bit_cost += this_gene;
		cBuff++;
		/*find a gene which represents an operator
		 if (bOperator)
		 {
		 if ( (this_gene < 10) || (this_gene > 13) )

		 continue;

		 else
		 {
		 bOperator       = false;
		 buffer[cBuff++] = this_gene;
		 continue;
		 }
		 }

		 //find a gene which represents a number
		 else
		 {
		 if (this_gene > 9)

		 continue;

		 else
		 {
		 bOperator       = true;
		 buffer[cBuff++] = this_gene;
		 continue;
		 }
		 }

		 }//next gene

		 //  now we have to run through buffer to see if a possible divide by zero
		 //  is included and delete it. (ie a '/' followed by a '0'). We take an easy
		 //  way out here and just change the '/' to a '+'. This will not effect the
		 //  evolution of the solution
		 for (int i=0; i<cBuff; i++)
		 {
		 if ( (buffer[i] == 13) && (buffer[i+1] == 0) )

		 buffer[i] = 10;
		 */}

	return cBuff;
}

//---------------------------------AssignFitness--------------------------------------
//
//  given a string of bits and a target value this function will calculate its
//  representation and return a fitness score accordingly
//------------------------------------------------------------------------------------
double AssignFitness(string bits, double target_value, Precision p) {
	//holds decimal values of gene sequence
	int buffer[(int) (CHROMO_LENGTH / GENE_LENGTH)];
	int num_elements = ParseBits(bits, buffer);
	/*
	 for (int v = 0; v < num_elements; v++) {
	 cout << buffer[v] << endl;
	 }
	 cout << endl;*/

	// ok, we have a buffer filled with valid values of: operator - number - operator - number..
	// now we calculate what this represents.
	double result = 0.0;

	for (int i = 0; i < num_elements; i++) {
		//result+=buffer[i];
		result = p.compute_error(buffer);
		/* switch (buffer[i])
		 {
		 case 10:

		 result += buffer[i+1];
		 break;

		 case 11:

		 result -= buffer[i+1];
		 break;

		 case 12:

		 result *= buffer[i+1];
		 break;

		 case 13:

		 result /= buffer[i+1];
		 break;

		 }//end switch
		 */
	}

	// Now we calculate the fitness. First check to see if a solution has been found
	// and assign an arbitarily high fitness score if this is so.

	//Add soft condition of termination
	if (result <= target_value) {
		//if (bit_cost < 85) {
		return 999999999.0;
	} else {
		//CONTROLLA E SISTEMA
		//double value = 1 / (double) fabs((double) (target_value - result));
		double value=target_value/(double)result;
		value = (double) value / pow(bit_cost,2);
		return value;
	}
	//  return result;
}

//---------------------------------PrintChromo---------------------------------------
//
// decodes and prints a chromo to screen
//-----------------------------------------------------------------------------------
void PrintChromo(string bits, Precision p) {
	//holds decimal values of gene sequence
	int buffer[(int) (CHROMO_LENGTH / GENE_LENGTH)];

	//parse the bit string
	int num_elements = ParseBits(bits, buffer);
	bit_cost = 0;
	for (int i = 0; i < num_elements; i++) {
		bit_cost += buffer[i];
		PrintGeneSymbol(buffer[i]);
	}
	total_costs.push_back(bit_cost);
	double error=p.compute_error(buffer);
	cout << "   bit cost: " << bit_cost;
	cout<<"  error: "<<error<<endl;
	return;
}

//--------------------------------------PrintGeneSymbol-----------------------------
//
//  given an integer this function outputs its symbol to the screen
//----------------------------------------------------------------------------------
void PrintGeneSymbol(int val) {
	cout << val << "-";
	/*if (val < 10 )

	 cout << val << " ";

	 else
	 {
	 switch (val)
	 {

	 case 10:

	 cout << "+";
	 break;

	 case 11:

	 cout << "-";
	 break;

	 case 12:

	 cout << "*";
	 break;

	 case 13:

	 cout << "/";
	 break;

	 }//end switch

	 cout << " ";
	 }
	 */
	return;
}

//------------------------------------Mutate---------------------------------------
//
//  Mutates a chromosome's bits dependent on the MUTATION_RATE
//-------------------------------------------------------------------------------------
void Mutate(string &bits) {    //Check to create values between 4 and 53
	for (int i = 0; i < bits.length(); i++) {
		if (RANDOM_NUM < MUTATION_RATE) {
			if (bits.at(i) == '1')

				bits.at(i) = '0';

			else

				bits.at(i) = '1';
		}
	}

	return;
}

//---------------------------------- Crossover ---------------------------------------
//
//  Dependent on the CROSSOVER_RATE this function selects a random point along the
//  lenghth of the chromosomes and swaps all the  bits after that point.
//------------------------------------------------------------------------------------
void Crossover(string &offspring1, string &offspring2) {
	//dependent on the crossover rate
	if (RANDOM_NUM < CROSSOVER_RATE) {
		//create a random crossover point
		int crossover = (int) (RANDOM_NUM * CHROMO_LENGTH);
		string t1 = offspring1.substr(0, crossover)
				+ offspring2.substr(crossover, CHROMO_LENGTH);
		string t2 = offspring2.substr(0, crossover)
				+ offspring1.substr(crossover, CHROMO_LENGTH);
		offspring1 = t1;
		offspring2 = t2;
	}
}

//--------------------------------Roulette-------------------------------------------
//
//  selects a chromosome from the population via roulette wheel selection
//------------------------------------------------------------------------------------
string Roulette(int total_fitness, chromo_typ* Population) {
	//generate a random number between 0 & total fitness count
	double Slice = (double) (RANDOM_NUM * total_fitness);

	//go through the chromosones adding up the fitness so far
	double FitnessSoFar = 0.0;

	for (int i = 0; i < POP_SIZE; i++) {
		FitnessSoFar += Population[i].fitness;

		//if the fitness so far > random number return the chromo at this point
		if (FitnessSoFar >= Slice)

			return Population[i].bits;
	}

	return "";
}
string RankSelection(chromo_typ* Population) {
	double array[POP_SIZE];
	for (int i = 0; i < POP_SIZE; i++) {
		array[i] = Population[i].fitness;
	}
	sort(array, array + POP_SIZE);
	double max = array[POP_SIZE - 1];
	for (int i = 0; i < POP_SIZE; i++) {
		if (Population[i].fitness == max) {
			return Population[i].bits;
		}
	}
	cout << "MEGAERROR" << endl;
	return "";
}
//- See more at: https://www.codemiles.com/c-examples/genetic-algorithm-example-t7548.html#sthash.3PFpVcvp.dpuf
