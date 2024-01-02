//
//  cpiegame.cpp
//  signalling
//
//  Created by Tom Lenaerts on 03/02/2020.
//  Copyright © 2020 Tom Lenaerts. All rights reserved.
//

#include "cpiegame.hpp"
#include <cmath>
#include <iomanip>
#define FIXED_FLOAT(x) std::fixed <<std::setprecision(2)<<(x)


CpieGame::CpieGame(const CpieGame& other){
    _p1=other.p1();
    _p2=other.p2();
    _factor=other.factor();
    _length=other.length();
    _payoffsp1.clear();
    _payoffsp2.clear();
    _payoffsp1.push_back(other.p1());
    _payoffsp2.push_back(other.p2());
    for(unsigned i = 1; i < (other.length()+1); i++){
        double part = _payoffsp1[i-1]* (1.0 -other.factor());
        double diff = _payoffsp1[i-1] - part;
        double prev = _payoffsp2[i-1];
            _payoffsp1.push_back(part);
            _payoffsp2.push_back(prev+diff);
    }
    for(unsigned i=1; i < _payoffsp1.size(); i++){
        if(i%2==0){
            double tmp = _payoffsp1[i];
            _payoffsp1[i]=_payoffsp2[i];
            _payoffsp2[i]= tmp;
        }
    }
}

CpieGame& CpieGame::operator=(const CpieGame& other){
    _p1=other.p1();
    _p2=other.p2();
    _factor=other.factor();
    _length=other.length();
    _payoffsp1.clear();
    _payoffsp2.clear();
    _payoffsp1.push_back(other.p1());
    _payoffsp2.push_back(other.p2());
    for(unsigned i = 1; i < (other.length()+1); i++){
        double part = _payoffsp1[i-1]* (1.0 -other.factor());
        double diff = _payoffsp1[i-1] - part;
        double prev = _payoffsp2[i-1];
            _payoffsp1.push_back(part);
            _payoffsp2.push_back(prev+diff);
    }
    for(unsigned i=1; i < _payoffsp1.size(); i++){
        if(i%2==0){
            double tmp = _payoffsp1[i];
            _payoffsp1[i]=_payoffsp2[i];
            _payoffsp2[i]= tmp;
        }
    }
    return *this;
}



std::ostream& CpieGame::display(std::ostream& os) const {
    os << "[[";
    for(unsigned i =0 ; i < _payoffsp1.size(); i++){
        os << _payoffsp1[i];
        if(i< (_payoffsp1.size()-1)){
            os<< ", ";
        }
    }
    os << "],\n[";
    for(unsigned i =0 ; i < _payoffsp2.size(); i++){
        os << _payoffsp2[i];
        if(i< (_payoffsp2.size()-1)){
            os<< ", ";
        }
    }
    os << "]]";
    return os;
}

//Payoff as defined in the Mckelvey and palfrey paper
double CpieGame::payoff(unsigned round,bool which){
    if(round >=0 and round <= _length){
        if(which)
            return _payoffsp1[round];
        return _payoffsp2[round];
    }
    return NAN;
}

double CpieGame::payoff(unsigned rowA, unsigned rowB, unsigned colA, unsigned colB){//A and B define roles
    unsigned minimumA =min(rowA, colB); //row moves first
    double tmp1 = payoff(minimumA, true); // so payoff for row as first player
    unsigned minimumB =min(rowB, colA); // col moves first
    double tmp2 = payoff(minimumB, false); // so payoff for row as second player
    
    //need to adjust payoffs depending on who took the money when
    if(rowA != colB){
        if(rowA< colB && (rowA%2)!=0)
            tmp1 = payoff(minimumA+1, true);
        else if(colB < rowA && (colB%2)==0)
            tmp1 = payoff(minimumA+1, true);
    }
    if(rowB != colA){
        if(colA< rowB && (colA%2)!=0)
            tmp2 = payoff(minimumB+1, false);
        else if(rowB < colA && (rowB%2)==0)
            tmp2 = payoff(minimumB+1, false);
    }
    return (tmp1 + tmp2)/2.0; // return average payoff in both role A and role B
}
