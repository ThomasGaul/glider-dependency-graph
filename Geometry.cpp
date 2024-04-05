#include "Geometry.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

using std::cout; using std::ofstream;

/* General Utilities */

template<typename EltType>
std::istream& operator>>(std::istream& is, vector<EltType>& v)
{
    for (auto& e_ : v) is >> e_;
    return is;
}

template<typename EltType>
std::ostream& operator<<(std::ostream& os, vector<EltType>& v)
{
    for (auto it = v.cbegin(); it != --v.cend(); ++it)
        os << *it << ' ';
    os << *(--v.cend());
    return os;
}

bool adjacentQ (short g1, short g2)
{
    switch (g1)
    {
    case 1:
        switch (g2)
        {
        case 2: return 1; break;
        case 4: return 1; break;
        case 5: return 1; break;
        default: return 0; break;
        }
        break;
    case 2:
        switch (g2)
        {
        case 1: return 1; break;
        case 3: return 1; break;
        case 4: return 1; break;
        case 5: return 1; break;
        default: return 0; break;
        }
        break;
    case 3:
        switch (g2)
        {
        case 2: return 1; break;
        case 5: return 1; break;
        case 6: return 1; break;
        default: return 0; break;
        }
        break;
    case 4:
        switch (g2)
        {
        case 1: return 1; break;
        case 2: return 1; break;
        case 5: return 1; break;
        case 7: return 1; break;
        case 8: return 1; break;
        default: return 0; break;
        }
        break;
    case 5:
        if (g2 > 0 && g2 < 10) return 1;
        else return 0;
        break;
    case 6:
        switch (g2)
        {
        case 2: return 1; break;
        case 3: return 1; break;
        case 5: return 1; break;
        case 8: return 1; break;
        case 9: return 1; break;
        default: return 0; break;
        }
        break;
    case 7:
        switch (g2)
        {
        case 4: return 1; break;
        case 5: return 1; break;
        case 8: return 1; break;
        default: return 0; break;
        }
        break;
    case 8:
        switch (g2)
        {
        case 4: return 1; break;
        case 5: return 1; break;
        case 6: return 1; break;
        case 7: return 1; break;
        case 9: return 1; break;
        default: return 0; break;
        }
        break;
    case 9:
        switch (g2)
        {
        case 5: return 1; break;
        case 6: return 1; break;
        case 8: return 1; break;
        default: return 0; break;
        }
        break;
    default:
        return 0;
        break;
    }
}

/* ---------------- */
/* Geometry Methods */
/* ---------------- */

// default constructor
geometry::geometry (void)
{
    x_size = 3; y_size = 3;
    deque<unsigned> tmp(x_size,0);
    space.resize(y_size,std::move(tmp));
    return;
}

// size constructor
geometry::geometry (unsigned x, unsigned y)
{
    x_size = x; y_size = y;
    deque<unsigned> tmp(x_size,0);
    space.resize(y_size,std::move(tmp));
    return;
}

// construct from link_grid
geometry::geometry (link_grid& lg,layout* lo)
{
    x_size = 3; y_size = 3;
    deque<unsigned> tmp(x_size,0);
    space.resize(y_size,std::move(tmp));
    for (std::size_t cell = 1; cell <= 9; ++cell) {
        if (lg[cell].from != 0) {
            space[(cell-1)/3][(cell-1)%3] = lg[cell].to;
            if (!lo->usedQ(lg[cell].to))
                lo->use_process(lg[cell].to);
        }
    }
    return;
}

void geometry::trim (void)
{
    // search rows
    deque<unsigned> empty_row(x_size,0);
    bool empty = 1;
    do {
        if (space[0] != empty_row) empty = 0;
        else {space.pop_front(); --y_size;}
    } while (empty);
    empty = 1;
    do {
        if (space[y_size-1] != empty_row) empty = 0;
        else {space.pop_back(); --y_size;}
    } while (empty);
    
    // search columns
    empty = 1;
    do {
        for (std::size_t r = 0; r < y_size; ++r)
            if (space[r][0] != 0) {empty = 0; break;}
        if (empty) {
            for (std::size_t r = 0; r < y_size; ++r)
                space[r].pop_front();
            --x_size;
        }
    } while (empty);
    empty = 1;
    do {
        for (std::size_t r = 0; r < y_size; ++r)
            if (space[r][x_size-1] != 0) {empty = 0; break;}
        if (empty) {
            for (std::size_t r = 0; r < y_size; ++r)
                space[r].pop_back();
            --x_size;
        }
    } while (empty);
    return;
}

// display to terminal
void geometry::bash_display()
{
    for (auto& row : space) {
        for (auto& p_ : row) {
            if (p_ == 0) cout << "_\t";
            else cout << p_ << '\t';
        }
        cout << '\n';
    }
}

/* ------------ */
/* Link Methods */
/* ------------ */

// default constructor
link::link (void)
{
    from = 0;
    to = 0;
    cell = 0;
    type = 0;
}

// link_grid constructor
link_grid::link_grid (void)
{
    grid.resize(9);
    return;
}

// stream operator
istream& operator>> (istream& is, link& l)
{
    is >> l.from;
    is >> l.to;
    is >> l.type;
    is >> l.cell;
    return is;
}

bool indicate13 = 0;

// check whether to processes can be valid neighbors
bool link_grid::neighborQ (process& p1,process& p2,short p1_n,short p2_n)
{
    short m_rel = relative_moore(p1_n,p2_n);
    if (!adjacentQ(p1_n,p2_n)) goto check_links;

    // check states
    if (p1[m_rel] != p2.state || p2[invert_cell(m_rel)] != p1.state)
        return 0;
    // check moore neighbors
    for (short n = 1; n <= 9; ++n) {
        if (!adjacentQ(n,m_rel)) continue;
        else if (p1[n] != p2[relative_moore(m_rel,n)])
            return 0;
    }
    // check for shared link
    check_links:
    for (auto& lnk1 : p1.links) {
        for (auto& lnk2 : p2.links) {
            if (lnk1.to != lnk2.to) continue;
            if (adjacentQ(p1_n,p2_n) && relative_moore(lnk1.cell,lnk2.cell) == m_rel)
                return 1;
            else if (lnk1.cell - lnk2.cell == p1_n - p2_n) return 1;
        }
    }
    return 0;
}

/* --------------- */
/* Process Methods */
/* --------------- */

// default constructor
process::process (void)
{
    index = 0;
    type = -2; state = -1;
    used_links = 0;
    moore.resize(8,-1);
    // loc.resize(2);
    // links.resize(9);
}

// stream operator
istream& operator>> (istream& is, process& p)
{
    is >> p.index;
    is >> p.type;
    p.state = type_to_val.at(p.type);
    is >> p.moore;
    return is;
}

short& process::operator[](short m)
{
    switch (m)
    {
    case 1: return moore[0]; break;
    case 2: return moore[1]; break;
    case 3: return moore[2]; break;
    case 4: return moore[3]; break;
    case 5: return state; break;
    case 6: return moore[4]; break;
    case 7: return moore[5]; break;
    case 8: return moore[6]; break;
    case 9: return moore[7]; break;
    default:
        std::cerr << "invalid moore index: " << m << '\n';
        exit(1);
        break;
    }
}

bool process::linksQ (unsigned index, short cell)
{
    for (auto& l : links)
        if (l.to == index && (cell == 0 || l.cell == cell))
            return 1;
    return 0;
}

// find link sets
void process::find_link_sets (vector<process>* ref_set)
{
    // check that links exist
    if (links.size() == 0) return;
    
    // group links by dependency cell
    list<vector<unsigned>> sets;
    for (unsigned l_idx = 0; l_idx < links.size(); ++l_idx) {
        bool match = 0;
        for (auto& set : sets) {
            if (links[set.front()].cell == links[l_idx].cell) {
                set.push_back(l_idx);
                match = 1;
            }
        }
        if (!match) {
            vector<unsigned> tmp({l_idx});
            sets.push_back(std::move(tmp));
        }
    }

    // find most repeated dependency cell
    short max_cell = 0, max = 0;
    list<vector<unsigned>>::iterator max_ptr;
    for (auto set = sets.begin(); set != sets.end(); ++set) {
        if (set->size() > max) {
            max_cell = invert_cell(links[set->front()].cell);
            max = set->size();
            max_ptr = set;
        }
    }

    // copy links to list (for efficient removal without replacement)
    list<link> links_no_r;
    links_no_r.insert(links_no_r.begin(),links.begin(),links.end());

    
    // initialize dependency grids
    unsigned num_rm = 0;
    // short init_loc = invert_cell(links[max_ptr->back()].cell);
    for (auto& l_idx : *max_ptr) {
        link_grid set;
        link_sets.push_back(std::move(set));
        link_sets.back()[max_cell] = links[l_idx];
        list<link>::iterator rm = links_no_r.begin();
        std::advance(rm,l_idx-num_rm);
        links_no_r.erase(rm);
        ++num_rm;
    }

    // fill grid adjacent to initial process
    for (std::size_t cell = 1; cell <= 9; ++cell) {
        for (auto& grid : link_sets) {
            for (auto lnk = links_no_r.begin(); lnk != links_no_r.end(); ++lnk) {
                // check basic compatability
                if (grid[max_cell].type != lnk->type || invert_cell(lnk->cell) != cell)
                    continue;
                // check with all adjacent neighbors
                bool valid_neighbor = 1;
                for (std::size_t n = 1; n <= 9; ++n) {
                    if (grid[n].to == 0)
                        continue;
                    if (!grid.neighborQ((*ref_set)[grid[n].to-1],
                                             (*ref_set)[lnk->to-1],n,invert_cell(lnk->cell)))
                    {valid_neighbor = 0; break;}
                }

                if (valid_neighbor) {
                    grid[cell] = *lnk;
                    lnk = links_no_r.erase(lnk);
                    --lnk;
                    break;
                }
            }
        }
    }
    return;
}


/* ------------------ */
/* Derivation Methods */
/* ------------------ */

// name constructor
layout::layout (const string object)
{
    // allocate memory based on object
    if (object == "block") {
        reference_set.resize(16);
        support_set.resize(16,0);
    }
    else if (object == "blinker") {
        reference_set.resize(30);
        support_set.resize(30,0);
    }
    else if (object == "glider") {
        reference_set.resize(352);
        support_set.resize(352,0);
    }
    else {
        std::cerr << "invalid object name ("
            << object << ") provided at initialization\n";
        exit(1);
    }

    // load in process definitions and dependency links
    std::ifstream is1("DG_files/"+object+"_process_set.dat"),
                  is2("DG_files/"+object+"_links.dat");
    string line;
    
    for (unsigned idx = 0; idx < reference_set.size(); ++idx) {
        is1 >> reference_set[idx];
    } while (std::getline(is2,line)) {
        std::istringstream iss(line);
        // load link
        link l; iss >> l;
        // assign link to process
        reference_set[l.from-1].links.push_back(std::move(l));
    }

    // find link sets
    vector<process>* ref_set_ptr = &reference_set;
    for (auto& p : reference_set)
        p.find_link_sets(ref_set_ptr);

    return;
}

// destructor
layout::~layout (void)
{
    // release geometry pointers
    return;
}

// initialize geometry algorithm
void layout::initialize (unsigned p_idx /*1-index*/)
{
    if (p_idx < 1) {
        std::cerr << "invalid index (must be 1-indexed)\n";
        exit(1);
    }
    shared_ptr<geometry> g_init = make_shared<geometry>(1,1);
    geometries.push_back(std::move(g_init));
    geometries.front()->space[0][0] = p_idx;

    use_links(p_idx);
    return;
}

// reset algorthim and geometries
void layout::reset (void)
{
    // release geometries
    for (auto& ptr : geometries) ptr.reset();
    geometries.clear();

    // reset support set
    support_set.resize(reference_set.size());
    unsigned index = 0;
    for (auto& p_ : support_set) p_ = ++index;

    // unmark used links
    for (auto& p_ : reference_set)
        p_.used_links = 0;

    return;
}

// remove process from support set
void layout::use_process (unsigned p_idx)
{
    support_set.remove(p_idx);
    return;
}

// check for process accessability
bool layout::usedQ (unsigned index)
{
    return !(std::binary_search(support_set.begin(),support_set.end(),index));
}

void layout::merge_search (void)
{
    unsigned mergers = 0;
    // compare all geometries
    for (auto g1 = geometries.begin(); g1 != geometries.end(); ++g1) {
        for (auto g2 = geometries.begin(); g2 != geometries.end(); ++g2) {
            if (g1 == g2) continue;
            bool merged = 0;
            for (std::size_t p1 = 0; p1 < (*g1)->size(); ++p1) {
                for (std::size_t p2 = 0; p2 < (*g2)->size(); ++p2) {
                    if ((*g1)->flat_index(p1) == (*g2)->flat_index(p2)
                            && (*g1)->flat_index(p1) != 0) {
                        merge(g1,g2,p1,p2);
                        g2 = geometries.erase(g2); --g2;
                        merged = 1;
                        break;
                    }
                }
                if (merged) break;
            }
        }
    }
    return;
}

void layout::merge (list<shared_ptr<geometry>>::iterator& g1,
                    list<shared_ptr<geometry>>::iterator& g2,
                    unsigned p1, unsigned p2)
{
    vector<unsigned> idx1({p1%(*g1)->x_size,p1/(*g1)->x_size}),
                     idx2({p2%(*g2)->x_size,p2/(*g2)->x_size});

    /* resize g1 */
    int front_col = idx2[0] - idx1[0],
        back_col = (*g2)->x_size - idx2[0] - ((*g1)->x_size - idx1[0]),
        front_row = idx2[1] - idx1[1],
        back_row = (*g2)->y_size - idx2[1] - ((*g1)->y_size - idx1[1]);
    // adjust values below zero
    if (front_col < 0) front_col = 0;
    if (back_col < 0) back_col = 0;
    if (front_row < 0) front_row = 0;
    if (back_row < 0) back_row = 0;
    
    // add columns at front
    for (int c = 0; c < front_col; ++c) {
        for (auto& row : (*g1)->space)
            row.push_front(0);
    }
    // add columns at back
    for (int c = 0; c < back_col; ++c) {
        for (auto& row : (*g1)->space)
            row.push_back(0);
    }
    (*g1)->x_size += ((front_col>0)?front_col:0) + ((back_col>0)?back_col:0);
    deque<unsigned> empty_row((*g1)->x_size,0);
    // add rows at front
    for (int r = 0; r < front_row; ++r)
        (*g1)->space.push_front(empty_row);
    // add rows at back
    for (int r = 0; r < back_row; ++r)
        (*g1)->space.push_back(empty_row);
    (*g1)->y_size += ((front_row>0)?front_row:0) + ((back_row>0)?back_row:0);

    // placement of g2 on new g1
    vector<int> overlay_point = {
        (int)idx1[0]-(int)idx2[0]+front_col,
        (int)idx1[1]-(int)idx2[1]+front_row
        };

    // merge
    for (std::size_t p2x = 0; p2x < (*g2)->x_size; ++p2x) {
        for (std::size_t p2y = 0; p2y < (*g2)->y_size; ++p2y) {
            unsigned p1x = overlay_point[0] + p2x,
                     p1y = overlay_point[1] + p2y;            
            if ((*g2)->space[p2y][p2x] == 0) continue;
            (*g1)->space[p1y][p1x] = (*g2)->space[p2y][p2x];
        }
    }

    // clean up
    g2->reset();

    return;
}

void layout::use_links (unsigned p_idx)
{
    for (auto& lg : reference_set[p_idx-1].link_sets) {
        shared_ptr<geometry> dep = make_shared<geometry>(lg,this);
        dep->trim();
        for (short p = 0; p < dep->size(); ++p) 
            use_process(dep->flat_index(p));
        geometries.push_back(std::move(dep));
    }
    reference_set[p_idx-1].used_links = 1;
    return;
}

bool layout::find (void)
{
    for (auto g = geometries.begin(); g != geometries.end(); ++g) {
        for (std::size_t p = 0; p < (*g)->size(); ++p) {
            unsigned p_idx = (*g)->flat_index(p);
            if (p_idx == 0 || reference_set[p_idx-1].used_links)
                continue;
            use_links(p_idx);
        }
        if (all_usedQ()) break;
    }
    merge_search();
    // for (auto& g : geometries)
    //     if (!check_bounding(g)) return 0;
    return 1;
}

// display to command line
void layout::bash_display(void)
{
    for (auto& g : geometries) {
        g->bash_display();
        cout << "\n\n";
    }
}

// check that all geometries are bounded
    // does not work yet
bool layout::check_bounding (shared_ptr<geometry>& g)
{    
    // from top
    unsigned top_cells = 3 * g->x_size;
    unsigned top_cnt = 0, steps = 0;
    for (std::size_t r = 0; r < g->y_size; ++r,++steps) {
        for (auto& p_ : g->space[r]) {
            if (p_ == 0) continue;
            // check top row of moore neighbors
            for (short n = 0; n < 3; ++n)
                if (reference_set[p_-1].moore[n] == 2)
                    ++top_cnt;
            if (top_cnt + steps == top_cells) break;
        }
        if (top_cnt + steps == top_cells) break;
    }

    // from bottom
    unsigned bot_cells = 3 * g->x_size;
    unsigned bot_cnt = 0; steps = 0;
    for (std::size_t r = g->y_size-1; r >= 0 && r < g->y_size; --r,++steps) {
        for (auto& p_ : g->space[r]) {
            if (p_ == 0) continue;
            // check top row of moore neighbors
            for (short n = 5; n < 8; ++n)
                if (reference_set[p_-1].moore[n] == 2)
                    ++bot_cnt;
            if (bot_cnt + steps == bot_cells) break;
        }
        if (bot_cnt + steps == bot_cells) break;
    }

    // from left
    unsigned left_cells = 3 * g->y_size;
    unsigned left_cnt = 0; steps = 0;
    for (std::size_t c = 0; c < g->x_size; ++c,++steps) {
        for (std::size_t r = 0; r < g->y_size; ++r) {
            unsigned p_ = g->space[r][c];
            if (p_ == 0) continue;
            if (reference_set[p_-1].moore[0] == 2)
                ++left_cnt;
            if (reference_set[p_-1].moore[3] == 2)
                ++left_cnt;
            if (reference_set[p_-1].moore[5] == 2)
                ++left_cnt;
            if (left_cnt + steps == left_cells) break;
        }
        if (left_cnt + steps == left_cells) break;
    }

    // from right
    unsigned right_cells = 3 * g->y_size;
    unsigned right_cnt = 0; steps = 0;
    for (std::size_t c = g->x_size-1; c >= 0 && c < g->x_size; --c,++steps) {
        for (std::size_t r = 0; r < g->y_size; ++r) {
            unsigned p_ = g->space[r][c];
            if (p_ == 0) continue;
            if (reference_set[p_-1].moore[2] == 2)
                ++right_cnt;
            if (reference_set[p_-1].moore[4] == 2)
                ++right_cnt;
            if (reference_set[p_-1].moore[7] == 2)
                ++right_cnt;
            if (right_cnt + steps == right_cells) break;
        }
        if (right_cnt + steps == right_cells) break;
    }

    return 1;
}

// check if all links have been used
bool layout::all_usedQ (void)
{
    for (auto& p : reference_set)
        if (!p.used_links) return 0;
    return 1;
}

void layout::display_to_file (string name)
{
    ofstream os("geometry/"+name+".dat");
    for (auto& g : geometries) {
        os << g->x_size << ' ' << g->y_size << '\n';
        for (std::size_t r = 0; r < g->y_size; ++r) {
            for (std::size_t c = 0; c < g->x_size; ++c) {
                unsigned p_ = g->space[r][c];
                if (p_ == 0) continue;
                os << c << ' ' << r << ' ';
                os << reference_set[p_-1].type << ' ';
                os << reference_set[p_-1].moore << '\n';
            }
        }
        os << "#\n";
    }
    os.close();
}