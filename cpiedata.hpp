//
//  cpiedata.hpp
//
//  Created by Tom Lenaerts on 11/02/2022.
//

#ifndef cpiedata_hpp
#define cpiedata_hpp

#include <stdio.h>
#include <vector>
#include <iostream>
#include <map>

using namespace std;

class CpieEntry{
public:
    CpieEntry(){}
    CpieEntry(unsigned steps):_steps(steps){
        _probs.clear();
        for(unsigned i=0;i<=steps;i++){
            _probs.push_back(0.0);
        }
    }
    CpieEntry(unsigned steps, vector<double> data):_steps(steps){
        _probs.clear();
        for(unsigned i=0;i<data.size();i++){
            _probs.push_back(data[i]);
        }
    }

    CpieEntry(const CpieEntry& other);
    CpieEntry& operator=(const CpieEntry& other);

    ~CpieEntry(){
        _probs.clear();
    }

    unsigned steps() const {return _steps;}
    unsigned int size() const {return (unsigned)_probs.size();}
    double operator[](unsigned pos) const;
    void set(unsigned pos, double val);
    
    bool operator==(const CpieEntry& other) const;
    bool operator!=(const CpieEntry& other) const;

    std::vector<double>::iterator begin() {return _probs.begin();}
    std::vector<double>::iterator end() {return _probs.end();}
    std::vector<double>::const_iterator begin() const{return _probs.begin();}
    std::vector<double>::const_iterator end() const {return _probs.end();}

    friend std::ostream & operator<<(std::ostream &o, CpieEntry& s){return s.display(o);}
        
protected:
    virtual ostream& display(ostream& os) const ;

    unsigned _steps;
    vector<double> _probs;
};

class CpieData {
public:
    CpieData(){
        _dict.clear();
    }
    ~CpieData(){
        _dict.clear();
    }
    
    void add(string key, CpieEntry* value);
    CpieEntry* operator[](string key) const;

    friend std::ostream & operator<<(std::ostream &o, CpieData& s){return s.display(o);}

protected:
    virtual ostream& display(ostream& os) const ;
    map<string,CpieEntry*> _dict;
    
};
#endif /* cpiedata_hpp */
