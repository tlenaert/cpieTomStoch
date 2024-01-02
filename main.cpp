//
//  main.cpp
//  signalling
//
//  Created by Tom Lenaerts on 23/01/2020.
//  Copyright © 2020 Tom Lenaerts. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "kernel.hpp"
#include "rangen.h"
#include "cpiegame.hpp"
#include "strategy.hpp"
#include <sstream>
#include <fstream>
#include <chrono>
#include <gsl/gsl_matrix.h>
#include "cpiedata.hpp"

using namespace std;
using namespace std::chrono;
#define FIXED_FLOAT(x) std::fixed <<setprecision(9)<<(x)


void runNetwork(unsigned psize, double cost, unsigned levels, unsigned maxlevel, double betas, double mut, double eps, unsigned repeats, double scaling, CpieGame& game, RanGen& ran, CpieData& data, ofstream& nodef, ofstream& edgef, ofstream& belf, ofstream& levf){
    
    StrategySpace strategies;
    strategies.createSymmetricUTStrategies(game, levels, maxlevel);
    cout << strategies << endl;

    high_resolution_clock::time_point start = high_resolution_clock::now();
    Kernel run(psize, strategies.size(), game);
    run.calcPayoffs(&strategies, eps, repeats, &ran, cost);
    run.execute(&strategies, betas, mut, scaling);
    high_resolution_clock::time_point stop = high_resolution_clock::now();
    duration<double> duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by function: "
         << duration.count() << " s" << endl;

    //belief transitions
    for(unsigned first=0; first <= game.length(); first++){
        for(unsigned second=0; second <= game.length(); second++){
            double trans = run.calcBeliefTransition(&strategies,5, first, second);
            belf << first << "\t" << second << "\t" << trans <<endl;
        }
    }
    //level transitions
    for(unsigned first=0; first <= levels; first++){
        for(unsigned second=0; second <= levels; second++){
            double trans = run.calcLevelTransition(&strategies,5, first, second);
            levf << first << "\t" << second << "\t" << trans <<endl;
        }
    }
    
    //nodes output
    run.printStationary(&strategies, 0.0, nodef); // include all SD
    //edges output (fixation probabilities)
    run.printFixation(&strategies,0.002,edgef); // links bigger than 0.002 (1/Z0
    
}


void runAnalytical(unsigned psize, double cost, unsigned levels, unsigned maxlevel, double betas, double mut, double eps, unsigned repeats, double scaling, CpieGame& game, RanGen& ran, CpieData& data, ofstream& of, ofstream& lf, ofstream& sf, ofstream& ef, ofstream& df, ofstream& levf, ofstream& mbf){
    
    StrategySpace strategies;
    strategies.createSymmetricUTStrategies(game, levels, maxlevel);
//    strategies.createSymmetricReducedStrategies(game, levels, maxlevel);
    cout << strategies << endl;

    CpieEntry* elm1 = data["MKP6cteavg"];


    high_resolution_clock::time_point start = high_resolution_clock::now();
    Kernel run(psize, strategies.size(), game);
    run.calcPayoffs(&strategies, eps, repeats, &ran, cost);
    run.execute(&strategies, betas, mut, scaling);
    high_resolution_clock::time_point stop = high_resolution_clock::now();
    duration<double> duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by function: "
         << duration.count() << " s" << endl;
    
    run.levelKdistribution(&strategies, betas, levels+1, lf);
    run.stepsdistribution(&strategies, betas, eps, levels+1, repeats, &ran, sf);
    run.decperkDistribution(&strategies, betas, levels+1, eps, repeats, &ran, levf);
    run.misbeliefperkDistribution(&strategies, betas, levels+1, eps, repeats, &ran, mbf);

    run.beliefDistribution(&strategies, betas, of);
    run.decisionDistribution(&strategies, betas, eps, repeats, &ran, df);
    //compare to actual data
    double MSE1 = run.calcMSE(&strategies,elm1, eps, repeats, &ran);
//    double MSE2 = run.calcMSE(&strategies,elm2, eps, repeats, &ran);
    ef << FIXED_FLOAT(betas) << "\t" << FIXED_FLOAT(eps) << "\t" << FIXED_FLOAT(MSE1) << "\t" << FIXED_FLOAT(sqrt(MSE1));
//    ef << "\t" << FIXED_FLOAT(MSE2) << "\t" << FIXED_FLOAT(sqrt(MSE2));
    ef << endl;
    
}


void runBeta(unsigned psize, double cost, unsigned levels, unsigned maxlevel, double mut, double eps, unsigned repeats, double scaling, CpieGame& game, RanGen& ran, CpieData& data, ofstream& of, ofstream& lf, ofstream& sf, ofstream& ef, ofstream& df, ofstream& levf){
    StrategySpace strategies;
    strategies.createSymmetricUTStrategies(game, levels, maxlevel);
//    strategies.createSymmetricReducedStrategies(game, levels, maxlevel);
    cout << strategies << endl;
    double tmp=0.0001;
    double basebeta = tmp;
    double beta=basebeta;
    unsigned step=0;
    unsigned interval = 1;
    CpieEntry* elm = data["MKP6cteavg"];

    Kernel run(psize, strategies.size(), game);
    run.calcPayoffs(&strategies, eps, repeats, &ran, cost);

    while(beta<=10.0){
        cout << interval << " - " << step << " - " << FIXED_FLOAT(beta) << endl;
        high_resolution_clock::time_point start = high_resolution_clock::now();
        run.execute(&strategies, beta, mut, scaling);

        run.levelKdistribution(&strategies, beta, levels+1, lf);
//        run.stepsdistribution(&strategies, beta, eps, levels+1, repeats, &ran, sf);
//        run.decisionDistribution(&strategies, beta, eps, repeats, &ran, df);
//        run.beliefDistribution(&strategies, beta, of);
        run.showEvoRobustStrategies(&strategies, 1.0);
        //compare to actual data
        double MSE = run.calcMSE(&strategies,elm, eps, repeats, &ran);

        ef << FIXED_FLOAT(beta) << "\t" << FIXED_FLOAT(eps) << "\t" << FIXED_FLOAT(MSE) << "\t" << FIXED_FLOAT(sqrt(MSE)) << endl;
        high_resolution_clock::time_point stop = high_resolution_clock::now();
        duration<double> duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by run function: "
             << fixed << duration.count() << " seconds" << endl;

        interval +=1;
        beta = (basebeta*interval);
        if(interval==10){
            interval=1;
            step +=1;
            basebeta = (tmp * pow(10, step));
        }
    }
}

void runEpsilon(unsigned psize, double cost, unsigned levels, unsigned maxlevel, double mut, double beta, unsigned repeats, double scaling, CpieGame& game, RanGen& ran, ofstream& of, ofstream& lf, ofstream& sf){
    StrategySpace strategies;
//    strategies.createEquilibriumStrategies(game);
    strategies.createSymmetricUTStrategies(game, levels, maxlevel);
    cout << strategies << endl;
    
    for(double eps=0.0 ; eps < 1.01  ; eps+=0.05){
        Kernel run(psize, strategies.size(), game);
        cout << FIXED_FLOAT(eps) << " - " << FIXED_FLOAT(beta) << endl;
        high_resolution_clock::time_point start = high_resolution_clock::now();
        run.calcPayoffs(&strategies, eps, repeats, &ran, cost);
        run.execute(&strategies, beta, mut, scaling);
        high_resolution_clock::time_point stop = high_resolution_clock::now();
        duration<double> duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by run function: "
             << fixed << duration.count() << " seconds" << endl;
        run.levelKdistribution(&strategies, eps, levels+1, lf);
//        run.beliefDistribution(&strategies, eps, of);
//        run.payoffResults(&strategies, eps, sf);
//        run.stepsdistribution(&strategies, beta, eps, levels+1, repeats, &ran, sf);

    }
}

void runBoth(unsigned psize, double cost, unsigned levels, unsigned maxlevel, double mut, unsigned repeats, double scaling, CpieGame& game, RanGen& ran, CpieData& data, ofstream& ef){
    StrategySpace strategies;
    strategies.createSymmetricUTStrategies(game, levels, maxlevel);
//    strategies.createSymmetricReducedStrategies(game, levels, maxlevel);
    cout << strategies << endl;
    CpieEntry* elm1 = data["MKP6cteavg"];
//    CpieEntry* elm2 = data["KT12"];

    for(double eps=0.5; eps <= 1.01; eps+=0.02){
        Kernel run(psize, strategies.size(), game);
        run.calcPayoffs(&strategies, eps, repeats, &ran,cost); //use the same payoff matrix for all beta
        
        for(double beta=0.0; beta <= 0.31; beta+=0.02){
            cout << "BETA=" << beta << "EPS=" << eps << endl;
            high_resolution_clock::time_point start = high_resolution_clock::now();
            run.execute(&strategies, beta, mut,scaling);
            high_resolution_clock::time_point stop = high_resolution_clock::now();
            duration<double> duration = duration_cast<microseconds>(stop - start);
            cout << "Time taken by function: "
                 << duration.count() << " s" << endl;

            //compare to actual data
            double MSE1 = run.calcMSE(&strategies,elm1, eps, repeats, &ran);
//            double MSE2 = run.calcMSE(&strategies,elm2, eps, repeats, &ran);
            double avgb = run.averageBelief(&strategies);
            double avgl = run.averageLevel(&strategies);
            double avgf = run.averageFitness(&strategies);
            double avgd = run.averageDecision(&strategies, repeats, eps, &ran);

            ef << FIXED_FLOAT(beta) << "\t" << FIXED_FLOAT(eps) << "\t" << FIXED_FLOAT(MSE1) << "\t" << FIXED_FLOAT(sqrt(MSE1));
//            ef << "\t" << FIXED_FLOAT(MSE2) << "\t" << FIXED_FLOAT(sqrt(MSE2));
            ef << "\t" << avgb << "\t" << avgl << "\t" << avgd << "\t" << (avgb-avgd) << "\t" << avgf << endl;
        }
    }
}


int main(int argc, char * argv[]) {
 
    //game specific for constant pie game
    unsigned length = 6;

    //constant pie
    double first = 1.6;
    double second = 1.6;
    double factor= 0.25;

    //strategy specific
    unsigned levels = 6;// or 6
    
    //reasoning stochasticity
    double epsilon=0.76;
    unsigned repeats=50000;

    
    //evolution specific
    unsigned psize = 500;
    double betas=0.02;
    double mut = 0.0;
    
//
    string lfname("ctplevelsLall.txt");
    string dfname("ctpbeliefsLall.txt");
    string sfname("ctpstepsLall.txt");
    string efname("ctperrorLall.txt");
    string decfname("ctpdecisionLall.txt");
    string levfname("ctpldecperKLall.txt");
    string mbfname ("ctpmbperKLall.txt");
    string belfname("ctpltransitbeliefs.txt");
    string nodefname("ctpnodes.txt");
    string edgefname("ctpedges.txt");

    
    cout << "Settings: ";
    cout << "psize(" << psize << "),  beta-s(" << betas << "), ";
    cout << "Length(" << length <<"), Factor(" << factor <<"), First(" << first <<"), Second(" << second <<"), Levels(" << levels <<"), " ;
    cout <<  "mut(" << mut << "), eps(" << epsilon << "), repeats("<<repeats<<"), " << endl;

    
    CpieGame game(first, second, factor, length);
    
    cout << game << endl;
    
    CpieData data;
    CpieEntry MKP6cteavg(6, {0.59, 0.33, 0.07, 0.007, 0.003, 0.0, 0.0});
    data.add("MKP6cteavg", &MKP6cteavg);

    cout << data;
    
    ofstream errf(efname);
    ofstream lf(lfname);
    ofstream sf(sfname);
    ofstream df(dfname);
    ofstream decf(decfname);
    ofstream levf(levfname);
    ofstream mbf(mbfname);
    ofstream belf(belfname);
    ofstream nodef(nodefname);
    ofstream edgef(edgefname);

    //
    
    RanGen ran;
    betas=0.04;
    epsilon=0.76;
    unsigned maxlevel=6;
    double scaling = 1.0;
    double cost = 0.0;
  
    
    
//    runNetwork(psize, cost, levels, maxlevel, betas, mut, epsilon, repeats, scaling, game, ran, data, nodef, edgef, belf, levf);
//    runAnalytical(psize, cost, levels, maxlevel, betas, mut, epsilon, repeats, scaling, game, ran, data, df, lf, sf, errf, decf, levf,mbf);
//    runBeta(psize, cost, levels, maxlevel, mut, epsilon, repeats, scaling, game, ran, data, df, lf, sf, errf, decf,levf);
    
    runEpsilon(psize, cost, levels, maxlevel, mut, betas, repeats, scaling, game, ran, df, lf, sf);

//    runBoth(psize, cost, levels, maxlevel, mut, repeats, scaling, game, ran, data, errf);

    errf.close();
    lf.close();
    sf.close();
    df.close();
    decf.close();
    levf.close();
    mbf.close();
    belf.close();
    nodef.close();
    edgef.close();
    return 0;
}

