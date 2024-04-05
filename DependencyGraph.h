/* ------------------------------------------ */
/*  Methods For Extracting Dependency Graphs  */
/* ------------------------------------------ */

#pragma once

#include <list>
#include <map>
#include <set>
#include <algorithm>
#include <sstream>

#include "LtL.h"

using std::list; using std::pair; using std::map; using std::set;
using std::istringstream;


/* --------------- */
/* Symmetry Method */
/* --------------- */

ostream& print_space_with_dont_cares (ostream&,Universe&);

void generate_symmetries(vector<string>&,rule&);


/* -------------------------- */
/* Process Analysis Functions */
/* -------------------------- */

vector<short> process_from_space (int,int,Universe&);

short type_from_state (bool,vector<short>&);

ofstream& save_non_null_processes (string,string,unsigned&,Universe&);

// takes list of base names and corresponding maximum index
void generate_process_set (vector<pair<string,unsigned>>&,rule&,string);


/* ----------------- */
/* Transition Search */
/* ----------------- */

void load_perturbation (Universe&,unsigned,string&);

void find_transitions(string,vector<pair<string,unsigned>>&,vector<string>&,rule&);


/* ----------------------------- */
/* Dependency Analysis Functions */
/* ----------------------------- */

unsigned index_lookup (vector<short>&,string,unsigned,unsigned,set<unsigned>&);
unsigned index_lookup (vector<short>&,string,unsigned,unsigned,vector<pair<unsigned,unsigned>>&,unsigned);

// returns type-first process from GoL universe
vector<short> cell_to_GoL_process (Universe&,int,int);

bool sort_vec (vector<unsigned>&,vector<unsigned>&);
void sort_links (string);

void generate_links (string,vector<pair<string,unsigned>>&,rule&);

void rotate (vector<short>&);
void flip (vector<short>&);
bool check_link_symmetry (string);


/* Utilities */

template<typename T>
inline ostream& operator<<(ostream& os, vector<T>& v)
{
    for (auto it = v.cbegin(); it != --v.cend(); ++it)
        os << *it << ' ';
    os << *(--v.cend());
    return os;
}
template<typename T>
inline istream& operator>>(istream& is,vector<T>& v)
{
    for (auto v_ = v.begin(); v_ != v.end(); ++v_)
        is >> *v_;
    return is;
}

inline string to_filename (string& base_name,unsigned& idx)
{
    return "Patterns/" + base_name + to_string(idx) + ".dat";
}