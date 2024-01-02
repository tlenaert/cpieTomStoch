//
//  cpiedata.cpp
//  ctpbeliefs
//
//  Created by Tom Lenaerts on 11/02/2022.
//

#include "cpiedata.hpp"
#include <cmath>

CpieEntry::CpieEntry(const CpieEntry& other){
    _steps=other.steps();
    _probs.clear();
    for(unsigned i=0; i < other.size(); i++){
        _probs.push_back(other[i]);
    }
}

CpieEntry& CpieEntry::operator=(const CpieEntry& other){
    _steps=other.steps();
    _probs.clear();
    for(unsigned i=0; i < other.size(); i++){
        _probs.push_back(other[i]);
    }
    return *this;
}


double CpieEntry::operator[](unsigned pos) const{
    if(pos>=0 && pos < _probs.size())
        return _probs[pos];
    return NAN;
}

void CpieEntry::set(unsigned pos, double val){
    if(pos>=0 && pos < _probs.size()){
        _probs[pos] = val;
    }
}

ostream& CpieEntry::display(ostream& os) const {
    os << "(" << _steps << "){";
    for(unsigned i=0; i < _probs.size(); i++){
        os << _probs[i];
        if(i < (_probs.size()-1))
            os <<",";
    }
    os << "}";
    return os;
}

bool CpieEntry::operator==(const CpieEntry& other) const{
    bool ret=(_steps == other.steps() && _probs.size() == other.size());
    if(!ret)
        return ret;
    
    for(unsigned i=0; i < _probs.size(); i++){
        ret = ret && (_probs[i] == other[i]);
        if(!ret)
            return ret;
    }
    return ret;
}

bool CpieEntry::operator!=(const CpieEntry& other) const{
    return !(*this == other);
}


void CpieData::add(string key, CpieEntry* value){
    _dict[key]=value;
}

CpieEntry* CpieData::operator[](string key) const{
    map<string,CpieEntry*>::const_iterator elm = _dict.find(key);
    CpieEntry* ret=NULL;
    if(elm != _dict.end())
        ret = elm->second;
    return ret;
}


ostream& CpieData::display(ostream& os) const {
    os << "DICT=["<< endl;
    map<string,CpieEntry*>::const_iterator start = _dict.begin();
    map<string,CpieEntry*>::const_iterator stop = _dict.end();
    while (start != stop) {
        string tmp = start->first;
        CpieEntry* e = start->second;
        os << "\t["<<tmp<<","<<*e<<"]";
        os << endl;
        start++;
    }
    os << "]" << endl;
    return os;
}



