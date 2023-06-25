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


// Mutation as the the sum of values following a normal distribution
/*
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
*/


// Mutation with chaotic sequences for exploration and exploitation
void Mutate(vector<double> & v){

	// variable for chaotic sequence
	static double r = Random::get(-1.0, 1.0);
	
	// iteration control
	static int iter =  0;
	static const int k = 0.48 * dim * 10000;
	
	// Exploration performed at early stages
	if(iter < k){	
		for(int i = 0; i < dim; i++){
			r = sin(20/r);
			double max = upperBound - v[i] > v[i] - lowerBound ? upperBound - v[i] : v[i] - lowerBound;
			v[i] += r * max;
		}
	}
	// Exploitation performed at late stages
	else{
		for(int i = 0; i < dim; i++){
			r = sin(20/r);
			double a = pow((10000 * dim - iter) / 10000 * dim, 2.0); 
			v[i] += r * (upperBound - lowerBound) * a / 2;
		}
	}
	
	// values are truncated if they are out of bounds
	for(int i = 0; i < dim; i++){
		if(v[i] > upperBound){
			v[i] = upperBound;
		} 
		else if (v[i] < lowerBound){
			v[i] = lowerBound;
		}
	}
	
	iter = (iter + 1) % (10000 * dim);
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
/*void Combine(const vector<double> & parent, vector<double> & larva){
	double alpha = Random::get(0.0, 1.0);
	
	for(int i = 0; i < dim; i++){
		larva[i] = parent[i] * alpha + (1-alpha)*larva[i];
	}
		
}
*/

//Crossover as chaotic aritmetic operator
void Combine(const vector<double> & parent, vector<double> & larva){
	//variable for chaotic sequence
	static double alpha = Random::get(-1.0, 1.0);
	
	for(int i = 0; i < dim; i++){
		alpha = 0.5 + 0.5 * sin(20 / alpha);
		larva[i] = parent[i] * alpha + (1-alpha)*larva[i];
	}
		
}


pair<vector<double>, vector<double>> Reproduce(const vector<double> & parent){
	
	pair<vector<double>, vector<double>> result;
	vector<double> bud = parent;
	
	Mutate(bud);
	
	//first element is the "larva" (only mutation)
	result.first = bud;
	
	Combine(parent, bud);
	
	//second element is the "bud" (mutation + crossover)
	result.second = bud;
	
	return result;
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
			pair<vector<double>, vector<double>> result = Reproduce(parent);
			
			// fitness of larva is calculated
			double fit_larva = cec17_fitness(&result.first[0]);
			evals++;
			
			// fitness of bud is calculated
			double fit_bud = cec17_fitness(&result.second[0]);
			evals++;
			
			// best of larva, bud and parent stays, the other two are discarded
			if(fit_larva < fitness) {
				fitness = fit_larva;
				parent = result.first;
			}
			
			// best of larva, bud and parent stays, the other two are discarded
			if(fit_bud < fitness) {
				fitness = fit_bud;
				parent = result.second;
			}
			
			
		}

		cout <<"Best ARO[F" <<funcid <<"]: " << scientific <<cec17_error(fitness) <<endl;
	}
    
}
