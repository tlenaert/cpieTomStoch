//
//  cpiegame.hpp
//  signalling
//
//  Created by Tom Lenaerts on 03/02/2020.
//  Copyright © 2020 Tom Lenaerts. All rights reserved.
//

#ifndef cpiegame_hpp
#define cpiegame_hpp

#include <string.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <cassert>

using namespace std;

class CpieGame {
public:
    CpieGame(double p1, double p2, double factor, unsigned length):_p1(p1), _p2(p2), _factor(factor),_length(length){
        _payoffsp1.push_back(p1);
        _payoffsp2.push_back(p2);
        for(unsigned i = 1; i < (length+1); i++){
            double part = _payoffsp1[i-1]* (1.0 -factor);
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
    };
    CpieGame(const CpieGame& other);
    
    double p1() const {return _p1;}
    bool setP1(double val) {_p1=val;return (_p1==val);}
    double p2() const {return _p2;}
    bool setP2(double val) {_p2=val;return (_p2==val);}
    double factor() const {return _factor;}
    bool setFactor(double val) {_factor=val;return (_factor==val);}
    unsigned length() const {return _length;}
    bool setLength(unsigned val) {_length=val;return (_length==val);}
    
    virtual double payoff(unsigned rounds, bool which);
    virtual double payoff(unsigned rowA, unsigned rowB, unsigned colA, unsigned colB);

    friend ostream & operator<<(ostream &o, CpieGame& g){return g.display(o);}
    virtual CpieGame& operator=(const CpieGame& other);

protected:
    virtual ostream& display(ostream& os) const ;

    
    bool check();
    double _p1;
    double _p2;
    double _factor;
    unsigned _length;
    
    vector<double> _payoffsp1;
    vector<double> _payoffsp2;
};


#endif /* cpiegame_hpp */
