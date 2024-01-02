//
//  kernel.hpp
//  signal
//
//  Created by Tom Lenaerts on 03/02/2020.
//  Copyright © 2020 Tom Lenaerts. All rights reserved.
//

#ifndef kernel_hpp
#define kernel_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>

#include "cpiegame.hpp"
#include "strategy.hpp"
#include "rangen.h"
#include "cpiedata.hpp"

#include <gsl/gsl_matrix.h>

using namespace std;

class Kernel {
public:
    Kernel (unsigned psize, unsigned ssize, CpieGame& game):_psize(psize), _ssize(ssize), _game(game) {
        _initialised=false;
        _fprobs = NULL;
        _trans = NULL;
        _stationary=NULL;
        _payoffs=NULL;
    };

    
    ~Kernel(){
        if(_fprobs != NULL)
            gsl_matrix_free(_fprobs);
        if(_trans != NULL)
            gsl_matrix_free(_trans);
        if(_stationary != NULL)
            gsl_vector_free(_stationary);
        if(_payoffs != NULL)
            gsl_matrix_free(_payoffs);
    }
    
    bool setGame(CpieGame& other) {
        _game=other;
        return true;
    }
    
    bool execute(StrategySpace* strats, double betas,double mut, double eps, unsigned repeats, RanGen* ran, double scaling);
    bool calcPayoffs(StrategySpace* strategies,double eps, unsigned repeats, RanGen* ran, double cost);
    bool execute(StrategySpace* strats, double betas, double mut,double scaling);
    
    void showEvoRobustStrategies(StrategySpace* strats, double mult);
    void showDomStrategies(StrategySpace* strats, double mult);
    void showStrategyNetwork(StrategySpace* strats, double mult);
    double calcLevelTransition(StrategySpace* strats, unsigned numbeliefs, unsigned from, unsigned to);
    double calcBeliefTransition(StrategySpace* strats, unsigned numlevels, unsigned from, unsigned to);
    void printFixation(StrategySpace* strats, double threshold, ofstream& of);
    void printStationary(StrategySpace* strats, double threshold, ofstream& of);

    bool levelKdistribution(StrategySpace* strats, double v1, unsigned levels, ofstream& of);
    bool stepsdistribution(StrategySpace* strats, double v1, double eps, unsigned steps, unsigned repeats, RanGen* ran, ofstream& of);
    bool decisionDistribution(StrategySpace* strats, double v1, double eps, unsigned repeats, RanGen* ran, ofstream& of);
    bool decperkDistribution(StrategySpace* strats, double v1, unsigned levels, double eps, unsigned repeats, RanGen* ran, ofstream& of);
    bool misbeliefperkDistribution(StrategySpace* strats, double v1, unsigned levels, double eps, unsigned repeats, RanGen* ran, ofstream& of);

    bool beliefDistribution(StrategySpace* strats, double v1, ofstream& of);  // same as full stationary distribution
    double averageBelief(StrategySpace* strats);
    double averageLevel(StrategySpace* strats);
    double averageFitness(StrategySpace* strats);

    double averageDecision(StrategySpace* strats, unsigned repeats, double eps, RanGen* ran);
    double calcMSE(StrategySpace* strats, CpieEntry* elm, double eps, unsigned repeats, RanGen* ran);
    friend ostream & operator<<(ostream &o, Kernel& k){return k.displayStationary(o);}

protected:
    bool initialize();
    
    double calcAvgPayoff(Strategy* first, Strategy* second, double eps, unsigned repeats, RanGen* ran);

    void calcPairwiseFitness(unsigned num_inv, Strategy* res, Strategy* inv, long double& res_fit, long double& inv_fit,double eps, unsigned repeats, RanGen* ran);
    void calcPairwiseFitness(unsigned num_inv, unsigned res, unsigned inv, long double& res_fit, long double& inv_fit);

    long double fermiFunction(bool sign_beta, long double first, long double second, double beta);
  
    void probIncreaseDecrease(unsigned num_inv, Strategy* res, Strategy* inv, long double& increase, long double& decrease, double betas, double mut, double eps, unsigned repeats, RanGen* ran);
    void probIncreaseDecrease(unsigned num_inv, unsigned res, unsigned inv, long double& increase, long double& decrease, double betas, double mut);

    long double fixationProbability(Strategy* res, Strategy* inv, double betas, double mut,double eps, unsigned repeats, RanGen* ran);
    long double fixationProbability(unsigned res, unsigned inv, double betas, double mut);

    void transitionMatrix(StrategySpace* strats, double betas, double mut, double eps, unsigned repeats, RanGen* ran, double scaling);
    void transitionMatrix(StrategySpace* strats, double betas, double mut, double scaling);

    bool createStationarySDistribution();
    bool createStationaryNSDistribution();

    long double clip(long double n, long double lower, long double upper);
    bool isEqual(double value1, double value2, double precision);

    void printMatrix(gsl_matrix* mat);
    std::ostream& displayStationary(std::ostream& os) const ;

    unsigned _psize;
    unsigned _ssize;
    bool _initialised;
    CpieGame _game;
    
    gsl_matrix* _payoffs;
    gsl_matrix* _fprobs;
    gsl_matrix* _trans;
    gsl_vector* _stationary;

};

#endif /* kernel_hpp */
