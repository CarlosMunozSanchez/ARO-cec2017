/* 
 * File:   main.cpp
 * Author: carlos
 *
 * Created on 26 de junio de 2023, 12:35
 */

extern "C" {
#include "cec17.h"
}
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <map>

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
		
		// initial population
		const int POPULATION_SIZE = 50;
		vector<vector<double>> population(POPULATION_SIZE);
		
		for(int i = 0; i < POPULATION_SIZE; i++){
			for(int j = 0; j < dim; j++){
				population[i].push_back(dis(gen));
			}
		}
		
		// individuals' indexes sorted by fitness
		multimap<double, int> ranking;
		// fitness for each individual
		vector<double> fitness_values(POPULATION_SIZE);
		
		// fitness for each individual is calculated
		for(int i = 0; i < POPULATION_SIZE; i++){
			double f = cec17_fitness(&population[i][0]);
			evals++;
			pair<double, int> p(f, i);
			
			ranking.insert(p);
			fitness_values[i] = f;
		}
		
		while(evals < 10000*dim){
			//cout << evals << endl;
			// Selection process: binary tournament
			int p1 = Random::get(0, POPULATION_SIZE - 1);
			int p2;
			do{
				p2 = Random::get(0, POPULATION_SIZE - 1);
			}while(p1 == p2);
			
			pair<vector<double>, vector<double>> result;
			  
			// reproduction: tournament winner reproduces
			if(fitness_values[p1] < fitness_values[p2]){
				result = Reproduce(population[p1]);
			} else{
				result = Reproduce(population[p2]);
				// swap p1 and p2 so that p1 contains the index of the individual that has reproduced
				swap(p1, p2);
			}
			
			// replacement: new individuals challenge worst individuals
			
			// larva and bud fitness are calculated
			double f1 = cec17_fitness(&result.first[0]);
			evals++;
			double f2 = cec17_fitness(&result.second[0]);
			evals++;
			
			// f1 will contain best from larva and bud
			if(f2 < f1){
				swap(f1, f2);
				swap(result.first, result.second);
			}
			
			// obtain two worst individuals
			multimap<double, int>::iterator worst, second_worst = ranking.end();
			second_worst--;
			worst = second_worst--;
			
			// new individuals challenge worst older individuals
			// best new vs 2nd worst old
			if(f1 < second_worst->first){
				//worst old is out, best new is in
				population[worst->second] = result.first;
				fitness_values[worst->second] = f1;
				ranking.emplace(f1, worst->second);
				// second best new vs worst old
				if(f2 < worst->first){
					// second worst old is also out, second best new is also in
					population[second_worst->second] = result.second;
					fitness_values[second_worst->second] = f2;
					ranking.emplace(f2, second_worst->second);
					ranking.erase(second_worst);
				}
				ranking.erase(worst);
			// if best new is not better than second worst old, best new vs worst old
			} else if(f1 < worst->first){
				// worst old is out, best new is in
				//worst old is out, best new is in
				population[worst->second] = result.first;
				fitness_values[worst->second] = f1;
				ranking.emplace(f1, worst->second);
				ranking.erase(worst);
			}			
			// else -> no new one is in
			
			
			//replacement: best one among the new individuals challenges the older one
			/*
			// larva and bud fitness are calculated
			double f1 = cec17_fitness(&result.first[0]);
			evals++;
			double f2 = cec17_fitness(&result.second[0]);
			evals++;
			
			// f1 will contain best from larva and bud
			if(f2 < f1){
				swap(f1, f2);
				swap(result.first, result.second);
			}
			
			// best new one vs his progenitor
			if(f1 < fitness_values[p1]){
				// if new is better, erase old
				ranking.erase(ranking.find(fitness_values[p1]));
				fitness_values[p1] = f1;
				population[p1] = result.first;
				ranking.emplace(f1, p1);
			}
			else{ // if new is worse, there is a chance it is still accepted
				double p = 1.0 - (evals / (10000.0*dim));
				p = p * p;
				if(Random::get<bool>(p)){
					ranking.erase(ranking.find(fitness_values[p1]));
					fitness_values[p1] = f1;
					population[p1] = result.first;
					ranking.emplace(f1, p1);
				}			
			}
			*/
		}
		
		double best = ranking.begin()->first;
		cout <<"Best ARO[F" <<funcid <<"]: " << scientific <<cec17_error(best) <<endl;
	}
    
}
