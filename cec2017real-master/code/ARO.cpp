/* 
 * File:   main.cpp
 * Author: carlos
 *
 * Created on 24 de junio de 2023, 12:35
 */

extern "C" {
#include "cec17.h"
}
#include <iostream>
#include <vector>
#include <random>
#include <cmath>

// From https://github.com/effolkronium/random
#include "random.hpp"

using namespace std;

// get base random alias which is auto seeded and has static API and internal state
// it is not threads secure, you can also use ::random_thread_local
using Random = effolkronium::random_static;

const double upperBound = 100.0, lowerBound = -100.0;

std::uniform_real_distribution<> dis(lowerBound, upperBound);
std::normal_distribution<double> norm_dist(0.0, 33.3);

int seed = 42;
int dim = 10;

//init. seed
std::mt19937 gen(seed); 

void Mutate(vector<double> & v){
	// a random number following a normal distribution is added to every component
	for(int i = 0; i < dim; i++){
		v[i] += norm_dist(gen);
		
		// values are truncated if they are out of bounds
		if(v[i] > upperBound){
			v[i] = upperBound;
		} 
		else if (v[i] < lowerBound){
			v[i] = lowerBound;
		}
	}
}

//Crossover as random selection from parent and larva
/*void Combine(const vector<double> & parent, vector<double> & larva){
	// larva's components are randomly taken from parent with 50% probability
	for(int i = 0; i < dim; i++){
		if( Random::get<bool>() ){
			larva[i] = parent[i];
		}
	}
}
*/

//Crossover as BLX-alpha crossover operator
/*void Combine(const vector<double> & parent, vector<double> & larva){
	double alpha = 0.3;
	
	for(int i = 0; i < dim; i++){
		double min = larva[i] < parent[i] ? larva[i] : parent[i];
		double max = larva[i] > parent[i] ? larva[i] : parent[i];
		double interval = max - min;
		
		larva[i] = Random::get(min - alpha*interval, max + alpha*interval);
	}
		
}*/

//Crossover as random arithmetic crossover operator
void Combine(const vector<double> & parent, vector<double> & larva){
	double alpha = Random::get(0.0, 1.0);
	
	for(int i = 0; i < dim; i++){
		larva[i] = parent[i] * alpha + (1-alpha)*larva[i];
	}
		
}


vector<double> Reproduce(const vector<double> & parent){
	vector<double> bud = parent;
	
	Mutate(bud);
	
	Combine(parent, bud);
	
	return bud;
}

int main() {

	for (int funcid = 1; funcid <= 30; funcid++) {
		cec17_init("ARO", funcid, dim);
		int evals = 0;
		
		//cerr <<"Warning: output by console, if you want to create the output file you have to comment cec17_print_output()" <<endl;
		//cec17_print_output(); // Comment to generate the output file
		
		//init. parent
		vector<double> parent(dim);
		for(int i = 0; i < dim; i++){
			parent[i] = dis(gen);
		}
		
		// fitness of parent is calculated
		double fitness = cec17_fitness(&parent[0]);
		evals++;
		
		while (evals < 10000*dim) {
			// parent reproduces a bud
			vector<double> bud = Reproduce(parent);
			
			// fitness of bud is calculated
			double fit_bud = cec17_fitness(&bud[0]);
			evals++;
			// if bud is better than parent, it replaces parent, otherwise it is discarded
			if(fit_bud < fitness) {
				fitness = fit_bud;
				parent = bud;
			}
		}

		cout <<"Best ARO[F" <<funcid <<"]: " << scientific <<cec17_error(fitness) <<endl;
	}
    
}
