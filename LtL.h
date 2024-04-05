//--------------------------------------
//            CA Algorithms
//--------------------------------------

#pragma once

#include <iostream>
#include <memory>
#include <vector>
#include <deque>
#include <string>
#include <cmath>

#include "random.h"

using std::vector; using std::deque; using std::string;
using std::unique_ptr; using std::shared_ptr;


typedef struct rule
{
    rule (vector<int> &vals, bool &i_m, bool &vN) {
        if (vals.size() != 5) cerr << "incorrect number of parameters\n";
        range = std::move(vals[0]);
        birth_min = std::move(vals[1]);
        birth_max = std::move(vals[2]);
        survive_min = std::move(vals[3]);
        survive_max = std::move(vals[4]);
        vonNeumann = std::move(vN);
    };
    rule (string file) {
        ifstream is("./Rules/"+file+".dat");
        is >> range;
        is >> birth_min; is >> birth_max;
        is >> survive_min; is >> survive_max;
        is >> vonNeumann;
        is.close();
    };
    rule(rule& rl) {
        birth_min = rl.birth_min;
        birth_max = rl.birth_max;
        survive_min = rl.survive_min;
        survive_max = rl.survive_max;
        range = rl.range;
        vonNeumann = rl.vonNeumann;
    };
    rule (void) {return;};
    int birth_min;
    int birth_max;
    int survive_min;
    int survive_max;
    int range;
    bool vonNeumann;
} rule;


class Universe
{
    public:
    Universe(vector<unsigned>&,float,vector<int>,rule&);
    Universe(unsigned,unsigned,rule&);
    Universe(unsigned,unsigned);
    Universe(Universe&);
    Universe(void) {return;};
    ~Universe(void) {rule_set.reset();};
    void share_load_point(vector<int>& lp)
        {pattern_load_point = lp;};
    void SetRandomSeed(long seed) {rs.SetRandomSeed(seed);};
    void set_rule(rule& rl) {rule_set = std::make_shared<rule>(rl);};
    void create_space(unsigned,unsigned);
    void resize_space(int,int,int,int);
    
    // modify
    void update(unsigned n = 1);
    void random_pattern(void);
    void random_pattern(float);
    void random_pattern(float,int,int,int,int);
    void empty(void);
    void rotate(short);
    void flip(void);
    void shift(int,int);
    void shift_with_crop(int,int);
    void set(int x, int y, bool val) {(*space)[y-lb_y][x-lb_x] = val;};
    void load_pattern(string filename)
        {ifstream is("Patterns/"+filename+".dat"); load_pattern(is);};
    void load_pattern_at(string filename,int at_x,int at_y,bool crop = 0)
        {ifstream is("Patterns/"+filename+".dat"); load_pattern_at(is,at_x,at_y,crop);};
    void load_pattern_toroidal(string filename,int at_x,int at_y)
        {ifstream is("Patterns/"+filename+".dat"); load_pattern_toroidal(is,at_x,at_y);};
    void load_pattern_vector(string,int,int,vector<bool>&);
    void load_bug_starter(int,int,int,int);

    // access
    bool operator()(int x, int y) const {
        if (x < lb_x || x >= ub_x || y < lb_y || y >= ub_y) {
            cerr << "space index out of bounds: (" << x << ", " << y << ")\n";
            exit(1);
        } return (*space)[y-lb_y][x-lb_x];
    };
    deque<deque<bool>>& get_space(void) {return *space;};
    void copy_space(Universe &u1) {*space=u1.get_space();};
    int get_xbound(void) {return ub_x;};
    int get_ybound(void) {return ub_y;};
    int get_xLbound(void) {return lb_x;};
    int get_yLbound(void) {return lb_y;};
    shared_ptr<rule>& get_rules(void) {return rule_set;};
    const long get_RandomSeed(void) {return rs.GetRandomSeed();};
    void save_pattern(string filename,int lby, int uby, int lbx, int ubx)
        {ofstream os("Patterns/"+filename+".dat"); save_pattern(os,lby,uby,lbx,ubx);};
    void save_pattern(string filename)
        {ofstream os("Patterns/"+filename+".dat"); save_pattern(os);};

    // analysis
    bool emptyQ(void);
    unsigned find_pattern(string);
    vector<int> find_pattern_loc(string,vector<int>&);
    unsigned find_cycle(string,unsigned);
    unsigned count_neighborhood(int x, int y) {
        if (x < lb_x || x >= ub_x || y < lb_y || y >= ub_y) {
            cerr << "space index out of bounds: (" << x << ", " << y << ")\n";
            exit(1);
        } return neighborhood_sum(x-lb_x,y-lb_y);
    };
    vector<double> density_vector(int,int);
    unsigned density_threshold(int,int);
    
    private:
    void update_private(void);
    unsigned neighborhood_sum(int,int);
    unsigned block_sum(int,int);
    unsigned vonNeumann_sum(int,int);
    void load_pattern(ifstream&);
    void load_pattern_at(ifstream&,int,int,bool);
    void load_pattern_toroidal(ifstream&,int,int);
    void save_pattern(ofstream&);
    void save_pattern(ofstream&,int,int,int,int);
    
    
    // space parameters
    shared_ptr<deque<deque<bool>>> space;
    int lb_x,ub_x,lb_y,ub_y;
    unsigned x_size,y_size;
    shared_ptr<rule> rule_set;

    // pattern parameters
    vector<int> pattern_load_point;
    RandomState rs;
    float random_bias;
};

// overload << operator to write universe to file/output
template<class EltType>
inline std::ostream &operator<<(std::ostream &os, const deque<deque<EltType>> &arr)
{
    // up y, right x
    for (auto it1 = arr.rbegin(); it1 != arr.rend(); ++it1) {
        for (auto it2 = it1->begin(); it2 != --(it1->end()); ++it2)
            os << *it2 << ' ';
        os << *(--(it1->end())) << '\n';
    }
    return os;
}

// overload == operator to compare universes
template<class EltType>
inline bool operator==(deque<deque<bool>> &u1, deque<deque<bool>> &u2)
{
    // check size match
    if (u1.size() != u2.size() || u1[0].size() != u2[0].size())
        return 0;
    for (auto y1=u1.begin(),y2=u2.begin(); y1!=u1.end(); ++y1,++y2) {
        for (auto x1=y1->begin(),x2=y2->begin();x1!=y1->end(); ++x1,++x2) {
            if (*x1!=*x2) return 0;
        }
    }
    return 1;
}

void operator++(vector<bool>& v);
bool last_binQ(vector<bool>&);
std::ostream& operator<<(std::ostream&,vector<bool>&);