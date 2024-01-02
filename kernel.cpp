//
//  kernel.cpp
//  signal
//
//  Created by Tom Lenaerts on 03/02/2020.
//  Copyright © 2020 Tom Lenaerts. All rights reserved.
//

#include "kernel.hpp"
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_randist.h>
#include <iomanip>
#include <sstream>
#include <limits.h>

#define FIXED_FLOAT(x) std::fixed << setprecision(6)<<(x)
#define SHORT_FIXED_FLOAT(x) std::fixed << setprecision(3)<<(x)
using namespace std;


long double Kernel::clip(long double n, long double lower, long double upper) {
  return std::max(lower, std::min(n, upper));
}

long double Kernel::fermiFunction(bool sign_beta, long double first, long double second, double beta){
    long double result = NAN;
    if (sign_beta)
        result = 1.0 / (1.0 + gsl_sf_exp((-beta) *(first-second)));
    else result = 1.0 / (1.0 + gsl_sf_exp((beta) *(first-second)));
    return result;
}

double Kernel::calcAvgPayoff(Strategy* first, Strategy* second, double eps, unsigned repeats, RanGen* ran){
    long double stochpayoff(0);
    for(unsigned i=0; i<repeats; i++){
        first->stochasticInferDecision(&_game, eps, ran);
        second->stochasticInferDecision(&_game, eps, ran);
//        cout << *first << "\t" << *second << "\t";
        long double tmp = _game.payoff(first->decisionR1(), first->decisionR2(), second->decisionR1(), second->decisionR2());
//        cout << tmp << endl;
        stochpayoff += tmp;
    }
    stochpayoff /=  double(repeats);
    return stochpayoff;
}

bool Kernel::calcPayoffs(StrategySpace* strategies,double eps, unsigned repeats, RanGen* ran, double cost){
    if(_payoffs== NULL){
        _payoffs = gsl_matrix_calloc(_ssize, _ssize);
    }
    gsl_matrix_set_zero(_payoffs);
    for(unsigned i =0; i < strategies->size(); i++){
        Strategy* one = (*strategies)[i];
        double lev = one->level();
        for(unsigned j =0; j < strategies->size(); j++){
            Strategy* two = (*strategies)[j];
            long double payoff = calcAvgPayoff(one, two, eps, repeats,ran);
            gsl_matrix_set(_payoffs, i, j, (payoff- (cost*lev)));
        }
    }
//    printMatrix(_payoffs);
    return true;
}



void Kernel::calcPairwiseFitness(unsigned num_inv, Strategy* res, Strategy* inv, long double& res_fit, long double& inv_fit,double eps, unsigned repeats, RanGen* ran){
    long double avg_inv_inv = calcAvgPayoff(inv, inv, eps, repeats,ran);
    long double avg_inv_res = calcAvgPayoff(inv, res, eps, repeats,ran);
    long double avg_res_inv = calcAvgPayoff(res, inv, eps, repeats,ran);
    long double avg_res_res = calcAvgPayoff(res, res, eps, repeats,ran);

    inv_fit= ((num_inv-1) * avg_inv_inv + (_psize - num_inv)* avg_inv_res) /  double(_psize-1);
    res_fit= (num_inv * avg_res_inv + (_psize - num_inv - 1)* avg_res_res) /  double(_psize-1);
}

void Kernel::calcPairwiseFitness(unsigned num_inv, unsigned res, unsigned inv, long double& res_fit, long double& inv_fit){
    long double avg_inv_inv = gsl_matrix_get(_payoffs, inv, inv);
    long double avg_inv_res = gsl_matrix_get(_payoffs, inv, res);
    long double avg_res_inv = gsl_matrix_get(_payoffs, res, inv);
    long double avg_res_res = gsl_matrix_get(_payoffs, res, res);

    inv_fit= ((num_inv-1) * avg_inv_inv + (_psize - num_inv)* avg_inv_res) /  double(_psize-1);
    res_fit= (num_inv * avg_res_inv + (_psize - num_inv - 1)* avg_res_res) /  double(_psize-1);
}

void Kernel::probIncreaseDecrease(unsigned num_inv, Strategy* res, Strategy* inv, long double& increase, long double& decrease, double betas, double mut, double eps, unsigned repeats, RanGen* ran){
    long double res_fit(0),inv_fit(0);
    calcPairwiseFitness(num_inv, res, inv, res_fit, inv_fit, eps, repeats, ran);
    increase = (((_psize - num_inv)/double(_psize))*(num_inv/double(_psize))*fermiFunction(false, res_fit, inv_fit, betas))*(1.0 - mut);
    increase += mut*((_psize - num_inv)/double(_psize));
    decrease = (((_psize - num_inv)/double(_psize))*(num_inv/double(_psize))*fermiFunction(true, res_fit, inv_fit, betas))*(1.0 - mut);
    decrease += mut*(num_inv/double(_psize));
}

void Kernel::probIncreaseDecrease(unsigned num_inv, unsigned res, unsigned inv, long double& increase, long double& decrease, double betas, double mut){
    long double res_fit(0),inv_fit(0);
    calcPairwiseFitness(num_inv, res, inv, res_fit, inv_fit);
    increase = (((_psize - num_inv)/double(_psize))*(num_inv/double(_psize))*fermiFunction(false, res_fit, inv_fit, betas))*(1.0 - mut);
    increase += mut*((_psize - num_inv)/double(_psize));
    decrease = (((_psize - num_inv)/double(_psize))*(num_inv/double(_psize))*fermiFunction(true, res_fit, inv_fit, betas))*(1.0 - mut);
    decrease += mut*(num_inv/double(_psize));
}


long double Kernel::fixationProbability(Strategy* res, Strategy* inv, double betas, double mut, double eps, unsigned repeats, RanGen* ran){
    long double result=0;
    for(unsigned i=0; i < _psize; i++ ){
        long double sub =1.0;
        for(unsigned j=1; j<(i+1) ; j++){
            long double increase(0),decrease(0);
            probIncreaseDecrease(j, res, inv, increase, decrease, betas, mut, eps, repeats, ran);
            sub *= (decrease/increase);
        }
        result += sub;
    }
    return (1.0/result); //clip((1.0/result),0.0, 1.0);
}

long double Kernel::fixationProbability(unsigned res, unsigned inv, double betas, double mut){
    long double result=0;
    for(unsigned i=0; i < _psize; i++ ){
        long double sub =1.0;
        for(unsigned j=1; j<(i+1) ; j++){
            long double increase(0),decrease(0);
            probIncreaseDecrease(j, res, inv, increase, decrease, betas, mut);
            sub *= (decrease/increase);
        }
        result += sub;
    }
    return (1.0/result); //clip((1.0/result),0.0, 1.0);
}

void Kernel::transitionMatrix(StrategySpace* strats, double betas, double mut, double eps, unsigned repeats, RanGen* ran, double scaling){
    for(unsigned i=0; i < _ssize ; i++){
        Strategy* res = (*strats)[i];
        double total=0;
        for(unsigned j=0 ; j < _ssize; j++){
            Strategy* inv = (*strats)[j];
            if(i != j){
                long double fp= fixationProbability(res, inv, betas, mut, eps, repeats, ran);
//                if(fp == 0){
//                    cout << "fix prob "<< *res << " to " << *inv << " is zero: " << scientific << fp << endl;
//                }
                if(res->level() != inv->level()){
                    fp = fp*scaling; //reduce the prob imitation of levels
                }
                gsl_matrix_set(_fprobs, i, j, fp);
                long double div =_ssize-1;
                gsl_matrix_set(_trans, i, j, (fp/div));
                total += gsl_matrix_get(_trans, i, j);
            }
        }
        gsl_matrix_set(_trans, i, i, (1.0 - total));

    }
//    printMatrix(_fprobs);
//    printMatrix(_trans);
}


void Kernel::transitionMatrix(StrategySpace* strats, double betas, double mut, double scaling){
    for(unsigned i=0; i < _ssize ; i++){
        double total=0;
        Strategy* res = (*strats)[i];
        for(unsigned j=0 ; j < _ssize; j++){
            Strategy* inv = (*strats)[j];
            if(i != j){
                long double fp= fixationProbability(i, j, betas, mut);
//                if(isEqual(fp, 0.0, 10)){
//                    cout << "fix prob "<< i << " to " << j << " is zero: " << scientific << fp << endl;
//                }
                if(res->level() != inv->level()){
                    fp = fp*scaling; //reduce the prob imitation of levels
                }
                gsl_matrix_set(_fprobs, i, j, fp);
                long double div =_ssize-1;
                gsl_matrix_set(_trans, i, j, (fp/div));
                total += gsl_matrix_get(_trans, i, j);
            }
        }
        gsl_matrix_set(_trans, i, i, (1.0 - total));
    }
//    printMatrix(_fprobs);
//    printMatrix(_trans);
}


bool Kernel::execute(StrategySpace* strats, double betas,double mut, double eps, unsigned repeats, RanGen* ran,double scaling){
    if(initialize()){
        //Calculate fixation probabilities and the transition matrix
        transitionMatrix(strats, betas, mut, eps, repeats, ran, scaling);
        //Determine now the stationary distribution for the transition matrix
        createStationaryNSDistribution();
        return true;
    }
    return false;
}


bool Kernel::execute(StrategySpace* strats, double betas,double mut, double scaling){
    if(initialize()){
        if(_payoffs != NULL){ //this function requires precalculation of payoffs
            //Calculate fixation probabilities and the transition matrix
            transitionMatrix(strats, betas, mut, scaling); // 
            //Determine now the stationary distribution for the transition matrix
            createStationaryNSDistribution();
    //        printStationary(_stationary);
    //        printMatrix(_fprobs);
            return true;
        }
    }
    return false;
}

void Kernel::showEvoRobustStrategies(StrategySpace* strats, double mult){
    double neutral = mult*(1.0/(_psize));
    for(unsigned i=0; i < _ssize ; i++){
        Strategy* strat = (*strats)[i];
        bool is_robust=true;
        for(unsigned j=0 ; j < _ssize; j++){
            if(i!=j && gsl_matrix_get(_fprobs, i, j) > neutral){ // no outgoing links
                is_robust = false;
//                cout << "from "<< *strat << " to " << *(*strats)[j] << " " << gsl_matrix_get(_fprobs, i, j) << endl;
                break;
            }
        }
        if(is_robust){
            cout << "Strategy "<< *strat << "is robust" << endl;
        }
    }
}

void Kernel::showStrategyNetwork(StrategySpace* strats, double mult){
    double neutral = mult*(1.0/(_psize));
    for(unsigned i=0; i < _ssize ; i++){
        Strategy* strat = (*strats)[i];
        for(unsigned j=0 ; j < _ssize; j++){
//            cout << "test " << i << "\t" << j << "\t"<< gsl_matrix_get(_fprobs, i, j) << endl;
            if(i!=j && gsl_matrix_get(_fprobs, i, j) > neutral){ // no outgoing links
                cout << "from "<< *strat << " to " << *(*strats)[j] << " " << gsl_matrix_get(_fprobs, i, j) << endl;
            }
        }
    }
}


void Kernel::showDomStrategies(StrategySpace* strats, double mult){
    double neutral = mult*(1.0/(_psize));
    for(unsigned j=0 ; j < _ssize; j++){
        Strategy* strat = (*strats)[j];
        bool is_dominated=true;
        for(unsigned i=0; i < _ssize ; i++){
            if(i!=j && gsl_matrix_get(_fprobs, i, j) > neutral){ // no incoming links
                is_dominated = false;
                break;
            }
        }
        if(is_dominated){
            cout << "Strategy "<< *strat << "is dominated" << endl;
        }
    }
}


void Kernel::printMatrix(gsl_matrix* mat){
    cout << "\t\t\t";
    for (unsigned i = 0; i < mat->size2; i++) {
        cout << fixed << i << "\t";
    }
    cout << endl;
    for (unsigned i = 0; i < mat->size1; i++) {
        cout << fixed << i << "\t";
        for (size_t j = 0; j < mat->size2; j++) {
            cout << scientific << gsl_matrix_get(mat, i, j);
            if( j < (mat->size2 - 1))
                cout<<"\t";
        }
        cout << endl;
    }
}

bool Kernel::isEqual(double value1, double value2, double precision){
    return std::fabs(value1 - value2) < std::pow(10, -precision);
}

bool Kernel::createStationarySDistribution(){ // symmetric transistion matrix
    gsl_matrix *transposed = gsl_matrix_calloc(_ssize, _ssize);
    gsl_matrix_transpose_memcpy(transposed, _trans);
    
    //transition matrix is not symmetric
    gsl_vector *eval = gsl_vector_alloc(_ssize);
    gsl_matrix *evec = gsl_matrix_alloc (_ssize,_ssize);
    
    
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (_ssize);
    gsl_eigen_symmv(transposed, eval, evec, w);

//    gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);

    int best_pos = -1;
    for(unsigned i=0; i < _ssize; i++){
//        cout << "Check :" << gsl_vector_get(eval,i) << endl;
        double tmp =fabs(gsl_vector_get(eval,i)) - 1.0;  //scaling paramter
        if (isEqual(tmp,0.0,6) && best_pos == -1){
            best_pos=i;
//            cout << "pos found "<< tmp << "\t" << best_pos << endl;
            break;
        }
//        gsl_vector_view test=gsl_matrix_column(evec, i);
//        for(unsigned j=0; j < _ssize; j++){
//            cout << "V" << j << " " << abs(gsl_vector_get(&test.vector, j)) << endl;
//        }
    }
//    cout << "column : " << best_pos << endl;
    
    gsl_vector_view selected=gsl_matrix_column(evec, best_pos);
    
    //normalize data in eigenvector
    double sum=0;
    for(unsigned i=0; i < _ssize; i++){
        gsl_vector_set(_stationary, i, abs(gsl_vector_get(&selected.vector, i)));
        sum+= gsl_vector_get(_stationary,i);
    }
    for(unsigned i=0; i < _ssize; i++){
        gsl_vector_set(_stationary, i,(gsl_vector_get(_stationary,i) / sum));
    }
    
    gsl_eigen_symmv_free (w);
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_matrix_free(transposed);
    return true;
}



bool Kernel::createStationaryNSDistribution(){
    gsl_matrix *transposed = gsl_matrix_calloc(_ssize, _ssize);
    gsl_matrix_transpose_memcpy(transposed, _trans);
    
    //transition matrix is not symmetric
    gsl_vector_complex *eval = gsl_vector_complex_alloc(_ssize);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc (_ssize,_ssize);
    
    
    gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (_ssize);
    gsl_eigen_nonsymmv(transposed, eval, evec, w);

    gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);

    int best_pos = 0;
    for(unsigned i=0; i < _ssize; i++){
        double tmp =fabs(GSL_REAL(gsl_vector_complex_get(eval,i)) - 1.0);
        if (isEqual(tmp,0.0,10) ){
            cout << "Found :" << GSL_REAL(gsl_vector_complex_get(eval,i)) << endl;
            if(best_pos == -1)
                best_pos=i;
            break;
        }
    }

    gsl_vector_complex_view selected=gsl_matrix_complex_column(evec, best_pos);
    
    //normalize data in eigenvector
    double sum=0;
    for(unsigned i=0; i < _ssize; i++){
        gsl_vector_set(_stationary, i, abs(GSL_REAL(gsl_vector_complex_get(&selected.vector, i))));
        sum+= gsl_vector_get(_stationary,i);
    }
    for(unsigned i=0; i < _ssize; i++){
        gsl_vector_set(_stationary, i,(gsl_vector_get(_stationary,i) / sum));
    }
    
    gsl_eigen_nonsymmv_free (w);
    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);
    gsl_matrix_free(transposed);
    return true;
}



std::ostream& Kernel::displayStationary(std::ostream& os) const {
    for (unsigned i = 0; i < _stationary->size; i++) {
        os <<  FIXED_FLOAT(gsl_vector_get(_stationary, i)) ;
        if (i < (_stationary->size - 1))
            os << "\t";
    }
    return os;
}


bool Kernel::initialize(){
    if(_initialised){
        gsl_matrix_free(_fprobs);
        gsl_matrix_free(_trans);
        gsl_vector_free(_stationary);
    }
    _fprobs = gsl_matrix_calloc(_ssize, _ssize);
    gsl_matrix_set_zero(_fprobs);
    _trans = gsl_matrix_calloc(_ssize, _ssize);
    gsl_matrix_set_identity(_trans);
    _stationary = gsl_vector_calloc(_ssize);
    gsl_vector_set_zero(_stationary);
    _initialised=true;
    return _initialised;
}

bool Kernel::levelKdistribution(StrategySpace* strats, double v1, unsigned levels, ofstream& of){
    map<unsigned,vector<double>> results;

    for(unsigned i=0; i < strats->size();i++){
        double tmp = gsl_vector_get(_stationary, i);
        Strategy* elm = (*strats)[i];
        map<unsigned,vector<double>>::iterator found =results.find(elm->level());
        if(found != results.end()){
            vector<double> data = found->second;
            data[elm->beliefR1()] += tmp;  //symmetric beliefs
            found->second = data;
        }
        else {
            vector<double> data(levels, 0.0);
            data[elm->beliefR1()] += tmp;
            results[elm->level()] = data;
        }
    }

    map<unsigned,vector<double>>::iterator start = results.begin();
    map<unsigned,vector<double>>::iterator stop = results.end();

    while (start!=stop){
        unsigned lev = start->first;
        of << v1 << "\t" << lev << "\t";
        double total = 0;
        vector<double> data = start->second;
        for(unsigned i=0; i <  data.size(); i++){
            of << data[i];
            total += data[i];
            if (i < (data.size()-1))
                of << "\t";
        }
        of << "\t" << total << endl;
        start++;
    }
    of << endl;
    return true;
}

bool Kernel::decperkDistribution(StrategySpace* strats, double v1, unsigned levels, double eps, unsigned repeats, RanGen* ran, ofstream& of){
    map<unsigned,vector<double>> results;

    for(unsigned i=0; i < strats->size();i++){
        double tmp = gsl_vector_get(_stationary, i);
        Strategy* elm = (*strats)[i];
        map<unsigned,vector<double>>::iterator found =results.find(elm->level());
        vector<double> collect(levels,0.0);
        for(unsigned j = 0; j < repeats; j++ ){
            Strategy first = *elm;
            Strategy second = *elm;
            first.stochasticInferDecision(&_game, eps, ran);
            second.stochasticInferDecision(&_game, eps, ran);
            unsigned whenp1=first.decisionR1();
            unsigned whenp2=second.decisionR2();
            map<unsigned,double>::iterator found;
            unsigned index=whenp1;
            if(whenp1 > whenp2)
                index=whenp2;
            collect[index] += tmp;
        }
        if(found != results.end()){
            vector<double> data = found->second;
            for(unsigned iter=0; iter < levels; iter++){
                double val = collect[iter]/double(repeats);
                data[iter] += val;
            }
            found->second = data;
        }
        else {
            for(unsigned iter=0; iter < levels; iter++){
                collect[iter]/=double(repeats);
            }
            results[elm->level()] = collect;
        }
    }

    map<unsigned,vector<double>>::iterator start = results.begin();
    map<unsigned,vector<double>>::iterator stop = results.end();

    while (start!=stop){
        unsigned lev = start->first;
        of << v1 << "\t" << lev << "\t";
        double total = 0;
        vector<double> data = start->second;
        for(unsigned i=0; i <  data.size(); i++){
            of << data[i];
            total += data[i];
            if (i < (data.size()-1))
                of << "\t";
        }
        of << "\t" << total << endl;
        start++;
    }
    of << endl;
    return true;
}

bool Kernel::misbeliefperkDistribution(StrategySpace* strats, double v1, unsigned levels, double eps, unsigned repeats, RanGen* ran, ofstream& of){
    map<unsigned,vector<double>> results;
    unsigned size = ((2*levels)-1);
    unsigned middle = (unsigned)floor(size/2);

    for(unsigned i=0; i < strats->size();i++){
        Strategy* elm = (*strats)[i];
        double tmp = gsl_vector_get(_stationary, i);
        map<unsigned,vector<double>>::iterator found =results.find(elm->level());
        vector<double> collect(size,0.0);
        for(unsigned j = 0; j < repeats; j++ ){
            Strategy first = *elm;
            Strategy second = *elm;
            first.stochasticInferDecision(&_game, eps, ran);
            second.stochasticInferDecision(&_game, eps, ran);
            unsigned whenp1=first.decisionR1();
            unsigned whenp2=second.decisionR2();
            map<unsigned,double>::iterator found;
            unsigned index=whenp1;
            if(whenp1 > whenp2)
                index=whenp2;
            int misbelief = (elm->beliefR1() - index); // symmetric belief
            collect[middle+misbelief] += tmp;
        }
        double test = 0;
        if(found != results.end()){
            vector<double> data = found->second;
            for(unsigned iter=0; iter < data.size(); iter++){
                double val = collect[iter]/double(repeats);
                data[iter] += val;
                test+=val;
            }
            found->second = data;
        }
        else {
            for(unsigned iter=0; iter < collect.size(); iter++){
                collect[iter]/=double(repeats);
                test+=collect[iter];
            }
            results[elm->level()] = collect;
        }
//        cout << *elm << "\t" << tmp << "\t" << test << endl;
    }

    map<unsigned,vector<double>>::iterator start = results.begin();
    map<unsigned,vector<double>>::iterator stop = results.end();
    while (start!=stop){
        unsigned lev = start->first;
        of << v1 << "\t" << lev << "\t";
        double pos(0), correct (0), neg(0), scaled(0), sum(0);
        vector<double> data = start->second;
        for(unsigned i=0; i <  data.size(); i++){
            sum+=data[i];
            if(i < middle) neg += data[i];
            if(i > middle) pos += data[i];
            if (i== middle) correct+=data[i];
            if(i < middle)
                scaled -= (data[i] * (levels-i-1));
            else scaled += (data[i] * (i-levels+1));
        }

        of << (correct/sum) << "\t" << (neg/sum) << "\t" << (pos/sum) << "\t" << scaled << endl;  //prints now the fraction not scaled to stationary distribution
        start++;
    }
    of << endl;
    return true;
}

bool Kernel::decisionDistribution(StrategySpace* strats, double v1, double eps, unsigned repeats, RanGen* ran, ofstream& of){
    map<unsigned,double> results;
    //first make all possible decisions.
    for(unsigned i=0; i < (_game.length()+1); i++){
        results[i] = 0;
    }
    for(unsigned i=0; i < strats->size();i++){
        double tmp = gsl_vector_get(_stationary, i);
        Strategy* elm = (*strats)[i];
        for(unsigned j = 0; j < repeats; j++ ){
            Strategy first = *elm;
            Strategy second = *elm;
            first.stochasticInferDecision(&_game, eps, ran);
            second.stochasticInferDecision(&_game, eps, ran);
            unsigned whenp1=first.decisionR1();
            unsigned whenp2=second.decisionR2();
            map<unsigned,double>::iterator found;
            if(whenp1 <= whenp2)
                found=results.find(whenp1);
            else found=results.find(whenp2);
            found->second +=tmp;
        }
    }

    map<unsigned,double>::iterator start = results.begin();
    map<unsigned,double>::iterator stop = results.end(); // I assume the iteration is always in the same order
    
    of << v1 << "\t" << eps << "\t";
    unsigned numb= (unsigned)results.size();
    unsigned count=0;
    while(start != stop){
        cout << start->first << "\t";
        of << (start->second/double(repeats));
        if(count < (numb-1))
            of << "\t";
        start++;
        count++;
    }
    of << endl;
    cout << endl;
    return true;
}

bool Kernel::beliefDistribution(StrategySpace* strats, double v1, ofstream& of){
    map<string,double> results;
    for(unsigned i=0; i < strats->size();i++){
        double tmp = gsl_vector_get(_stationary, i);
        Strategy* elm = (*strats)[i];
        stringstream ss;
        ss << "(" << elm->beliefR1() << "," << elm->beliefR2() <<")";
        map<string,double>::iterator found =results.find(ss.str());
        if(found != results.end()){
            found->second +=tmp;
        }
        else results[ss.str()] = tmp;
    }

    map<string,double>::iterator start = results.begin();
    map<string,double>::iterator stop = results.end(); // I assume the iteration is always in the same order
    
    of << v1 << "\t";
    unsigned numb= (unsigned)results.size();
    unsigned count=0;
    while(start != stop){
        cout << start->first << "\t";
        of << start->second;
        if(count < (numb-1))
            of << "\t";
        start++;
        count++;
    }
    of << endl;
    cout << endl;
    return true;
}

double Kernel::averageBelief(StrategySpace* strats){
    vector<double> results (_game.length()+1,0.0) ;
    for(unsigned i=0; i < strats->size();i++){
        double tmp = gsl_vector_get(_stationary, i);
        Strategy* elm = (*strats)[i];
        results[elm->beliefR1()] += tmp;
    }
    double total =0;
    for (unsigned i=0; i < results.size(); i++){
        total += (i*results[i]);
    }
    return total;
}

double Kernel::averageLevel(StrategySpace* strats){
    vector<double> results (_game.length()+1,0.0) ;
    for(unsigned i=0; i < strats->size();i++){
        double tmp = gsl_vector_get(_stationary, i);
        Strategy* elm = (*strats)[i];
        results[elm->level()] += tmp;
    }
    double total =0;
    for (unsigned i=0; i < results.size(); i++){
        total += (i*results[i]);
    }
    return total;
}

double Kernel::averageDecision(StrategySpace* strats, unsigned repeats, double eps, RanGen* ran){
    vector<double> results (_game.length()+1,0.0) ;
    for(unsigned i=0; i < strats->size();i++){
        double tmp = gsl_vector_get(_stationary, i);
        Strategy* elm = (*strats)[i];
        for(unsigned j = 0; j < repeats; j++ ){
            Strategy first = *elm;
            Strategy second = *elm;
            first.stochasticInferDecision(&_game, eps, ran);
            second.stochasticInferDecision(&_game, eps, ran);
            unsigned whenp1=first.decisionR1();
            unsigned whenp2=second.decisionR2();
            if(whenp1 <= whenp2)
                results[whenp1] += (tmp);
            else results[whenp2] += (tmp);
        }
    }
    double total =0;
    for (unsigned i=0; i < results.size(); i++){
        double tmp = results[i]/repeats;
        total += (i*tmp);
    }
    return total;
}


bool Kernel::stepsdistribution(StrategySpace* strats, double v1, double eps, unsigned steps, unsigned repeats, RanGen* ran, ofstream& of){
    vector<double> results(steps,0.0);
    
    for(unsigned i=0; i < strats->size();i++){
        double tmp = gsl_vector_get(_stationary, i);
        Strategy* elm = (*strats)[i];
        for(unsigned j=0; j < repeats; j++){
            Strategy first = *elm;
            Strategy second = *elm;
            first.stochasticInferDecision(&_game, eps, ran);
            second.stochasticInferDecision(&_game, eps, ran);
            unsigned whenp1=first.decisionR1();
            if(whenp1%2!=0 && whenp1 != _game.length())
                whenp1+=1;
            unsigned whenp2=second.decisionR2();
            if(whenp2%2==0 && whenp2 != _game.length())
                whenp2+=1;
            if(whenp1 <= whenp2)
                results[whenp1] += tmp;
            else results[whenp2] += tmp;
        }
    }
    of << v1 << "\t" << eps << "\t";
    double sum=0;
    for (unsigned j=0; j < results.size(); j++){
        of << (results[j]/double(repeats));
        sum+=results[j];
        if (j < (results.size()-1))
            of << "\t";
    }
    of << endl;
//    cout << "Sum is " << sum << endl;
    return true;
}


double Kernel::calcMSE(StrategySpace* strats, CpieEntry* elm, double eps, unsigned repeats, RanGen* ran){
    double result = 0.0;
    unsigned len = elm->steps()+1;
//    cout << *elm << endl;
    vector<double> stpvec(len,0.0);
    for(unsigned i=0; i < strats->size();i++){
        double tmp = gsl_vector_get(_stationary, i);
        Strategy* s = (*strats)[i];
        for(unsigned j=0; j < repeats; j++){
            Strategy first = *s;
            Strategy second = *s;
            first.stochasticInferDecision(&_game, eps, ran);
            second.stochasticInferDecision(&_game, eps, ran);
            unsigned whenp1=first.decisionR1();
            if(whenp1%2!=0 && whenp1 != _game.length())
                whenp1+=1;
            unsigned whenp2=second.decisionR2();
            if(whenp2%2==0 && whenp2 != _game.length())
                whenp2+=1;
            if(whenp1 <= whenp2)
                stpvec[whenp1] += tmp;
            else stpvec[whenp2] += tmp;
        }
    }
    for(unsigned i=0 ; i < (len); i++){
        double tmp=(stpvec[i]/double(repeats));
        double diff = (tmp-(*elm)[i]);
        result += pow(diff, 2.0);
    }
    return (result/double(len));
}

double Kernel::averageFitness(StrategySpace* strats){

    double total = 0;
    for(unsigned i=0; i < strats->size();i++){
        double sdstrat1 = gsl_vector_get(_stationary, i);
        for(unsigned j=0; j < strats->size(); j++){
            double sdstrat2 = gsl_vector_get(_stationary, i);
            double payoff = gsl_matrix_get(_payoffs, i, j);
            total += (payoff*sdstrat1*sdstrat2);
        }
    }
    return (total/2.0);
}


double Kernel::calcLevelTransition(StrategySpace* strats, unsigned numbeliefs, unsigned from, unsigned to){
    double result(0);
    vector<double> partial(numbeliefs,0.0);

    for(unsigned i=0; i < strats->size();i++){
        Strategy* first = (*strats)[i];
        if(first->level() == from){
            double sdi = gsl_vector_get(_stationary, i);
            for (unsigned j=0; j < _fprobs->size2; j++){
                Strategy* second = (*strats)[j];
                if(second->level() == to){
                    double tmp = gsl_matrix_get(_fprobs, i, j)*sdi;
                    partial[second->beliefR1()] += tmp; //symmetric beliefs
                }
            }
        }
    }
    for(unsigned i=0; i < numbeliefs; i++){
        result += partial[i];
    }
    
    return result;
}

double Kernel::calcBeliefTransition(StrategySpace* strats, unsigned numlevels, unsigned from, unsigned to){
    double result(0);
    vector<double> partial(numlevels,0.0);

    for(unsigned i=0; i < strats->size();i++){
        Strategy* first = (*strats)[i];
        if(first->beliefR1() == from){
            double sdi = gsl_vector_get(_stationary, i);
            for (unsigned j=0; j < _fprobs->size2; j++){
                Strategy* second = (*strats)[j];
                if(second->beliefR1() == to){
                    double tmp = gsl_matrix_get(_fprobs, i, j)*sdi;
                    partial[second->level()] += tmp; //symmetric beliefs
                }
            }
        }
    }
    for(unsigned i=0; i < numlevels; i++){
        result += partial[i];
    }
    
    return result;
}

void Kernel::printStationary(StrategySpace* strats, double threshold, ofstream& of){
    of << "id,b,k,sd\n";
    for (unsigned i = 0; i < _stationary->size; i++) {
        Strategy* s = (*strats)[i];
        double val = gsl_vector_get (_stationary,i);
        if(val >= threshold)
            of << int(i) <<"," << int(s->beliefR1())<< "," << int(s->level())<< "," << scientific << val << endl;
    }
}

void Kernel::printFixation(StrategySpace* strats, double threshold, ofstream& of){
    of << "from,to,weight\n";
    for (unsigned i = 0; i < _fprobs->size1; i++) {
        for (unsigned j = 0; j < _fprobs->size2; j++) {
            double val = gsl_matrix_get(_fprobs, i, j);
            if(val > threshold)
                of << int(i) << "," << int(j) << "," << scientific << val << endl;
        }
    }
}
