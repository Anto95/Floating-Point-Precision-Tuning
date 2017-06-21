#include <string>
#include "utils.h"

struct chromo_typ
{
    //the binary bit string is held in a std::string
  string    bits;

    double     fitness;

    chromo_typ(): bits(""), fitness(0.0f){};
    chromo_typ(string bts, double ftns): bits(bts), fitness(ftns){}
};


/////////////////////////////////prototypes/////////////////////////////////////////////////////
int* 	genetic(Precision p);
void    PrintGeneSymbol(int val);
string  GetRandomBits(int length);
int     BinToDec(string bits);
string 	DecToBin(int dec);
double   AssignFitness(string bits, double target_value,Precision p);
void    PrintChromo(string bits,Precision p);
void    PrintGeneSymbol(int val);
int     ParseBits(string bits, int* buffer);
string  Roulette(int total_fitness, chromo_typ* Population);
string  RankSelection(chromo_typ* Population);
void    Mutate(string &bits);
void    Crossover(string &offspring1, string &offspring2);
