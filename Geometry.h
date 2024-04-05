#pragma once

#include <iostream>
#include <vector>
#include <deque>
#include <list>
#include <memory>
#include <map>
#include <string>
#include <initializer_list>

using std::vector; using std::map; using std::deque; using std::list;
using std::shared_ptr; using std::unique_ptr; using std::make_shared;
using std::string; using std::istream;

/* --------- */
/* Utilities */
/* --------- */

const map<short,short> type_to_val = {
    {-2,-1}, {-1,2}, {0,0}, {1,0}, {2,0}, {3,1}, {4,1}, {5,1}
};
const map<short,short> val_to_type = {
    {-1,-2}, {0,0}, {1,3}, {2,-1}
};
const map<short,short> type_to_nextval = {
    {-2,-1}, {-1,-1}, {1,1}, {2,0}, {4,0}, {5,1}
};
const map<short,vector<short>> nine_to_3x3 = {
    {1,{0,0}}, {2,{1,0}}, {3,{2,0}},
    {4,{0,1}}, {5,{1,1}}, {6,{2,1}},
    {7,{0,2}}, {8,{1,2}}, {9,{2,2}}
};


/* declarations */
struct geometry;
struct link;
struct link_grid;
struct process;
class layout;

/* --------------- */
/* Geometry Object */
/* --------------- */
typedef struct geometry {
    geometry(void);
    geometry(unsigned,unsigned);
    geometry(link_grid&,layout*);
    ~geometry(void) {return;};
    inline unsigned& operator[](vector<unsigned> idx) {
        return space[idx[1]][idx[0]];
    };
    inline unsigned size(void) {return x_size*y_size;};
    inline unsigned& flat_index(unsigned i /* 0-index */)
        {return space[i/x_size][i%x_size];};
    void trim(void);
    void bash_display(void);
    deque<deque<unsigned>> space;
    unsigned x_size, y_size;
    // bool done;
} geometry;

/* ---------------- */
/*   Link Objects   */
/* ---------------- */
typedef struct link {
    link(void);
    link& operator=(link& l) {
        if (this == &l) return *this;
        from = l.from;
        to = l.to;
        cell = l.cell;
        type = l.type;
    };
    ~link(void) {
        from = 0; to = 0;
    };
    unsigned from, to;
    short cell, type;
} link;
istream& operator>>(istream&,link&);

inline short invert_cell(short cell) {return 10-cell;}

inline short relative_moore(short g1,short g2) {return g2 + (5-g1);}

bool adjacentQ(short,short);

typedef struct link_grid {
    link_grid(void);
    ~link_grid(void) {return;};
    bool neighborQ(process&,process&,short,short);
    inline link& operator[](short idx) {
        return grid[idx-1];
    };
    vector<link> grid;
} link_grid;

/* -------------- */
/* Process Object */
/* -------------- */
typedef struct process {
    process(void);
    ~process(void) {return;};
    bool linksQ(unsigned,short);
    void find_link_sets(vector<process>*);
    short& operator[](short);
    unsigned index;
    short type, state;
    bool used_links;
    vector<short> moore;
    vector<link> links;
    list<link_grid> link_sets;
} process;
istream& operator>>(istream&,process&);

/* ---------------- */
/* Derivation Class */
/* ---------------- */
class layout {
    public:
    layout(const string);
    ~layout(void);
    void initialize(unsigned);
    bool find(void);
    void reset(void);

    void display_to_file(string name);
    void bash_display(void);
    unsigned support_set_size(void) {return reference_set.size();};
    process& support_set_idx(unsigned i) {return reference_set[i];};
    unsigned num_geometries(void) {return geometries.size();};

    bool usedQ(unsigned);
    void use_process(unsigned);
    
    private:
    void use_links(unsigned);
    void merge_search(void);
    void merge(list<shared_ptr<geometry>>::iterator&,
                list<shared_ptr<geometry>>::iterator&,
                unsigned,unsigned);
    bool check_bounding(shared_ptr<geometry>&);
    bool all_usedQ(void);

    list<shared_ptr<geometry>> geometries;
    list<unsigned> support_set;
    vector<process> reference_set;
};


/* General Utilities */
template<typename EltType>
std::istream& operator>>(std::istream&,vector<EltType>&);
template<typename EltType>
std::ostream& operator<<(std::ostream&,vector<EltType>&);