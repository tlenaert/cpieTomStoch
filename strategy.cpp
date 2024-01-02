//
//  strategy.cpp
//  signal
//
//  Created by Tom Lenaerts on 03/02/2020.
//  Copyright © 2020 Tom Lenaerts. All rights reserved.
//

#include <algorithm>
#include "strategy.hpp"
#include <iomanip>

#define TWO_FIXED_FLOAT(x) std::fixed <<setprecision(2)<<(x)


Strategy::Strategy(const Strategy& other){
    _level= other.level();
    _belief1=other.beliefR1();
    _belief2=other.beliefR2();
    _decision1=other.decisionR1();
    _decision2=other.decisionR2();
}

Strategy& Strategy::operator=(const Strategy& other){
    _level= other.level();
    _belief1=other.beliefR1();
    _belief2=other.beliefR2();
    _decision1=other.decisionR1();
    _decision2=other.decisionR2();
    return *this;
}

ostream& Strategy::display(ostream& os) const {
    os << "[B" << _belief1 << "-L"<< _level << "]";
    return os;
}

bool Strategy::operator==(const Strategy& other) const{ // same when beliefs are the same
    return (_level == other.level() && _belief1 == other.beliefR1() && _belief2==other.beliefR2());
    //the decision given the beliefs may vary as there is stochasticity given the beliefs
}

bool Strategy::operator!=(const Strategy& other) const{
    return !(*this == other);
}

unsigned Strategy::belief(unsigned role) const{
    if(role==1)
        return _belief1;
    return _belief2;
}

unsigned Strategy::decision(unsigned role) const{
    if(role==1)
        return _decision1;
    return _decision2;
}

bool Strategy::setBelief(unsigned role, unsigned value){
    if(role==1){
        _belief1 = value;
        return true;
    }
    else if(role==2){
        _belief2 = value;
        return true;
    }
    return false;
}

bool Strategy::setDecision(unsigned role, unsigned value){
    if(role==1){
        _decision1 = value;
        return true;
    }
    else if(role==2){
        _decision2 = value;
        return true;
    }
    return false;
}

bool Strategy::inferR2Decision(CpieGame* game, unsigned& R1, unsigned& R2){
    double R1pi = game->payoff(R1, false);
    double R1pi_alt = game->payoff(R1-1, false); // what do de when R1 is 0 ? not chnage R2
    if(R1pi_alt > R1pi && R1 > 0)
       R2=R1-1;
    return true;
}

bool Strategy::inferR1Decision(CpieGame* game, unsigned& R1, unsigned& R2){
    double R2pi = game->payoff(R2, true);
    double R2pi_alt = game->payoff(R2-1, true);
    if(R2pi_alt > R2pi && R2 > 0)
       R1=R2-1;
    return true;
}

bool Strategy::inferDecision(CpieGame* game){
    unsigned tmpR1 = ((_belief1%2!=0 && _belief1 < game->length())?_belief1+1:_belief1); //when can one take
    unsigned tmpR2 = ((_belief2 !=game->length() && _belief2%2==0)?_belief2+1:_belief2); //when can two take
    unsigned level=0;
    while(level < _level){
        unsigned testR1(tmpR1);
        inferR1Decision(game, tmpR1, tmpR2); //modifies tempR1
        inferR2Decision(game, testR1, tmpR2); // modifies tempR2
        level+=1;
    }
    _decision2=tmpR2;
    _decision1=tmpR1;
    return true;
}

int Strategy::addNoise(RanGen* ran, CpieGame* game, int base, double epsilon){
    double test = ran->randouble();
    int res=base;
    if(test < epsilon){
        vector<int> vals={-1,1};
        int step=ran->ranval(0, 1);
        res = base + vals[step];
        if(res < 0) res=0;
        else if(res > game->length()) res=game->length();
    }
    return res;
}


bool Strategy::stochasticInferDecision(CpieGame* game, double epsilon, RanGen* ran){
    unsigned tmpR1 = _belief1;
    unsigned tmpR2 = _belief2;  
    unsigned level=0;
    while(level < _level){
        unsigned testR1(tmpR1);
        inferR1Decision(game, tmpR1, tmpR2); //modifies R1
        //add noise to tempR1
        tmpR1 = addNoise(ran, game, tmpR1, epsilon);
        inferR2Decision(game, testR1, tmpR2); // modifies R2
        //add noise to tempR2
        tmpR2 = addNoise(ran, game, tmpR2, epsilon);
        level+=1;
    }
    _decision2=(tmpR2%2==0 && tmpR2 != game->length()?tmpR2+1:tmpR2);
    _decision1=(tmpR1%2!=0 && tmpR1 != game->length()?tmpR1+1:tmpR1);
    return true;
}


///StrategySpace methods
bool StrategySpace::createRandNowakStrategies(CpieGame& game, unsigned levels){
    if(_cleared){
        //all L0 strategies
        for(unsigned i = 0; i <= (levels); i++){
            for(unsigned j = 0; j <= (levels); j++){
                Strategy *elm = new Strategy(0, i, j);
                elm->inferDecision(&game);
                _space.push_back(elm);
            }
        }
        _cleared = false;
    }
    return _cleared;
}

bool StrategySpace::createSubRandNowakStrategies(CpieGame& game, unsigned levels){
    if(_cleared){
        //all L0 strategies
        for(unsigned i = 0; i <= (levels); i+=2){
            for(unsigned j = 0; j <= (levels); j+=2){
                Strategy *elm = new Strategy(0, i, j);
                elm->inferDecision(&game);
                _space.push_back(elm);
            }
        }
        _cleared = false;
    }
    return _cleared;
}

//bool StrategySpace::createAllStrategies(CpieGame& game, unsigned levels){
//    if(_cleared){
//        //all L0 strategies
//        for(unsigned i = 0; i <= (levels); i+=1){
//            for(unsigned i = 0; i <= (levels); i+=1){
//                for(unsigned j = 0; j <= (levels); j+=1){
//                    Strategy *elm = new Strategy(0, i, j);
//                    elm->inferDecision(&game);
//                    _space.push_back(elm);
//                }
//            }
//        }
//        _cleared = false;
//    }
//    return _cleared;
//}

bool StrategySpace::createAltruistStrategies(CpieGame& game, unsigned levels){
    if(_cleared){
        //all L0 strategies
        for(unsigned i = 0; i <= (levels); i++){
            Strategy *elm = new Strategy(i, levels, levels);
            elm->inferDecision(&game);
            _space.push_back(elm);
        }
        _cleared = false;
    }
    return _cleared;
}

bool StrategySpace::createAllSymmetricStrategies(CpieGame& game, unsigned levels){
    if(_cleared){
        for(unsigned i = 0; i <= (levels); i++){
            for(unsigned j = 0; j <= (levels); j++){
                Strategy *elm = new Strategy(j, i, i);
                elm->inferDecision(&game);
                _space.push_back(elm);
            }
        }
        _cleared = false;
    }
    return _cleared;
}

bool StrategySpace::createSymmetricStrategies(CpieGame& game, unsigned beliefs, unsigned level){
    if(_cleared){
        for(unsigned i = 0; i <= (beliefs); i++){
            Strategy *elm = new Strategy(level, i, i);
            elm->inferDecision(&game);
            _space.push_back(elm);
        }
        _cleared = false;
    }
    return _cleared;
}

bool StrategySpace::createSymmetricUTStrategies(CpieGame& game, unsigned beliefs, unsigned level){
    if(_cleared){
        for(unsigned i = 0; i <= (beliefs); i++){
            for(unsigned j=0; j <=level; j++ ){
                Strategy *elm = new Strategy(j, i, i);
                elm->inferDecision(&game);
                _space.push_back(elm);
            }
        }
        _cleared = false;
    }
    return _cleared;
}

bool StrategySpace::createSymmetricReducedStrategies(CpieGame& game, unsigned beliefs, unsigned level){
    if(_cleared){
        for(unsigned i = 0; i <= (beliefs); i++){
            Strategy* prev = NULL;
            for(unsigned j=0; j <=level; j++ ){
                Strategy *elm = new Strategy(j, i, i);
                elm->inferDecision(&game);
                if(prev == NULL){
                    prev = elm;
                    _space.push_back(elm);
                }
                else if (prev->decisionR1() != elm->decisionR1() || prev->decisionR2() != elm->decisionR2()){
                    _space.push_back(elm);
                    prev = elm;
                }
            }
        }
        _cleared = false;
    }
    return _cleared;
}

bool StrategySpace::createNecessaryStrategies(CpieGame& game, unsigned levels){
    if(_cleared){
        //all L0 strategies
        for(unsigned i = 0; i <= (levels); i+=2){
            for(unsigned j = 0; j <= (levels); j+=2){
                for(unsigned l = 0; l <= (levels); l++){
                    Strategy *elm = new Strategy(l, i, j);
                    elm->inferDecision(&game);
                    _space.push_back(elm);
                }
            }
        }
        _cleared = false;
    }
    return _cleared;
}

bool StrategySpace::createRestrictedStrategies(CpieGame& game, unsigned levels){
    if(_cleared){
        //all L0 strategies
        for(unsigned i = 0; i <= (levels); i+=2){
            for(unsigned j = 0; j <= (levels); j+=2){
                Strategy* prev = NULL;
                for(unsigned l = 0; l <= (levels); l++){
                    Strategy *elm = new Strategy(l, i, j);
                    elm->inferDecision(&game);
                    if(prev == NULL){
                        prev = elm;
                        _space.push_back(elm);
                    }
                    else if (prev->decisionR1() != elm->decisionR1() || prev->decisionR2() != elm->decisionR2()){
                        _space.push_back(elm);
                        prev = elm;
                    }
                }
            }
        }
        _cleared = false;
    }
    return _cleared;
}


Strategy* StrategySpace::operator[](unsigned pos) const{
    if (!_cleared && pos >=0 && pos < _space.size())
        return _space[pos];
    return NULL;
}

struct compareStrategies{
    Strategy _tf;
    compareStrategies(Strategy tf): _tf(tf){};
    
    bool operator()(Strategy* elm){
        return (*elm == _tf);
    }
};

int StrategySpace::find(Strategy& s){
    int loc=-1;
    vector<Strategy*>::iterator itr = std::find_if(_space.begin(),_space.end(), compareStrategies(s));
    if(itr != _space.cend())
        loc=(unsigned) std::distance(_space.begin(), itr);
    return loc;
}

struct compareToString{
    unsigned _level;
    unsigned _b1;
    unsigned _b2;
    compareToString(unsigned level,unsigned b1, unsigned b2): _level(level), _b1(b1), _b2(b2){};
    
    bool operator()(Strategy* elm){
        return (elm->level() == _level && elm->beliefR1() == _b1 && elm->beliefR2() == _b2);
    }
};

int StrategySpace::find(unsigned level, unsigned t1, unsigned t2){
    int loc=-1;
    vector<Strategy*>::iterator itr = std::find_if(_space.begin(),_space.end(), compareToString(level, t1, t2));
    if(itr != _space.cend())
        loc=(unsigned) std::distance(_space.begin(), itr);
    return loc;

}


ostream& StrategySpace::display(std::ostream& os) const {
    os << "[";
    for (unsigned i=0 ; i< _space.size(); i++){
        os << *_space[i];
        if(i < (_space.size() -1))
            os <<",";
    }
    os << "]";
    return os;
}


