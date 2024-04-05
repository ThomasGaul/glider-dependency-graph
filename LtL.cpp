#include "LtL.h"

// standard constructor
Universe::Universe (vector<unsigned> &bounds, float prob, vector<int> p_load, rule &rl)
{
    // check for errors in passed arguments
    if (bounds.size()!=2) {std::cerr << "incorrect number of bound parameters\n"; exit(1);}
    // initialize universe
    ub_x = bounds[0];
    ub_y = bounds[1];
    lb_x = 0; lb_y = 0;
    x_size = ub_x - lb_x;
    y_size = ub_y - lb_y;
    pattern_load_point = p_load;
    random_bias = prob;
    deque<bool> line(bounds[0],0);
    shared_ptr<deque<deque<bool>>> buffer(new deque<deque<bool>>(y_size,std::move(line)));
    space = std::move(buffer);
    rule_set = std::make_shared<rule>(rl);
    return;
}

// basic constructor
Universe::Universe (unsigned ubx, unsigned uby, rule& rl)
{
    ub_x = ubx;
    ub_y = uby;
    lb_x = 0; lb_y = 0;
    x_size = ub_x - lb_x;
    y_size = ub_y - lb_y;
    deque<bool> line(x_size,0);
    shared_ptr<deque<deque<bool>>> buffer(new deque<deque<bool>>(y_size,std::move(line)));
    space = std::move(buffer);
    rule_set = std::make_shared<rule>(rl);
    random_bias = 0.5;
    pattern_load_point = {0,0};
    return;
}

// space constructor
Universe::Universe (unsigned ubx, unsigned uby)
{
    ub_x = ubx;
    ub_y = uby;
    lb_x = 0; lb_y = 0;
    x_size = ub_x - lb_x;
    y_size = ub_y - lb_y;
    deque<bool> line(ubx,0);
    shared_ptr<deque<deque<bool>>> buffer(new deque<deque<bool>>(y_size,std::move(line)));
    space = std::move(buffer);
    return;
}

// copy constructor
Universe::Universe (Universe& u)
{
    ub_x = u.get_xbound();
    ub_y = u.get_ybound();
    lb_x = u.get_xLbound();
    lb_y = u.get_yLbound();
    x_size = ub_x - lb_x;
    y_size = ub_y - lb_y;
    random_bias = u.random_bias;
    pattern_load_point = u.pattern_load_point;
    if (u.get_rules() != NULL)
        rule_set = u.get_rules();
    {
        deque<bool> line(x_size,0);
        shared_ptr<deque<deque<bool>>> buffer(new deque<deque<bool>>(y_size,std::move(line)));
        space = std::move(buffer);
    } copy_space(u);
    return;
}

// set space parameters (after null construction)
void Universe::create_space (unsigned x,unsigned y)
{
    ub_x = x;
    ub_y = y;
    lb_x = 0; lb_y = 0;
    x_size = ub_x - lb_x;
    y_size = ub_y - lb_y;
    deque<bool> tmp(x,0);
    space->resize(y,std::move(tmp));
    return;
}

// resize space
void Universe::resize_space (int xl,int xu,int yl,int yu)
{
    // remove rows
    for (int i = 0; i < -yu; ++i)
        (*space).pop_back();
    for (int i = 0; i < yl; ++i)
        (*space).pop_front();

    // add and remove columns
    for (int i = 0; i < space->size(); ++i) {
        // lower bound
        for (int j = 0; j < -xl; ++j)
            space->operator[](i).push_front(0);
        for (int j = 0; j < xl; ++j)
            space->operator[](i).pop_front();
        // upper bound
        for (int j = 0; j < xu; ++j)
            space->operator[](i).push_front(0);
        for (int j = 0; j < -xu; ++j)
            space->operator[](i).pop_front();
    }
    
    ub_x += xu; lb_x += xl;
    x_size = ub_x - lb_x;

    // add rows
    deque<bool> tmp(x_size,0);
    for (int i = 0; i < yu; ++i)
        (*space).push_back(tmp);
    for (int i = 0; i < -yl; ++i)
        (*space).push_front(tmp);

    ub_y += yu; lb_y += yl;
    y_size = ub_y - lb_y;
    return;
}


// update universe
unsigned Universe::neighborhood_sum (int x,int y)
{
    unsigned sum;
    if (rule_set->vonNeumann)
        sum = vonNeumann_sum(x,y);
    else
        sum = block_sum(x,y);
    return sum;
}
unsigned Universe::block_sum (int x,int y)
{
    unsigned sum = 0;
    
    for (int dy = -rule_set->range; dy <= rule_set->range; ++dy) {
        for (int dx = -rule_set->range; dx <= rule_set->range; ++dx) {
            sum += (unsigned)(*space)[(y_size+y+dy)%y_size][(x_size+x+dx)%x_size];
        }
    }
    return sum;
}
unsigned Universe::vonNeumann_sum (int x,int y)
{
    unsigned sum = 0;
    
    for (int dy = -rule_set->range; dy <= rule_set->range; ++dy) {
        for (int dx = -rule_set->range; dx <= rule_set->range; ++dx) {
            if (abs(dx)+abs(dy) > rule_set->range) continue;
            sum += (unsigned)(*space)[(y_size+y+dy)%y_size][(x_size+x+dx)%x_size];
        }
    }
    return sum;
}


void Universe::update_private (void)
{
    // construct pointer to buffer space
    deque<bool> tmp(x_size,0);
    unique_ptr<deque<deque<bool>>> buffer(new deque<deque<bool>>(y_size,std::move(tmp)));

    // function pointer to neighborhood type
    unsigned (Universe::*nsum)(int,int) = NULL;
    nsum = (rule_set->vonNeumann) ? &Universe::vonNeumann_sum : &Universe::block_sum;

    // update grid
    for (int y = 0; y < y_size; ++y) {
        for (int x = 0; x < y_size; ++x) {
            unsigned sum = (*this.*nsum)(x,y);
            if ((*space)[y][x] && (sum >= rule_set->survive_min) && (sum <= rule_set->survive_max))
                (*buffer)[y][x] = 1;
            else if (!(*space)[y][x] && (sum >= rule_set->birth_min) && (sum <= rule_set->birth_max))
                (*buffer)[y][x] = 1;
            else (*buffer)[y][x] = 0;
        }
    }
    space = std::move(buffer);
    return;
}
void Universe::update (unsigned n)
{
    for (std::size_t t = 0; t < n; ++t)
        update_private();
    return;
}


// generate random pattern in universe
void Universe::random_pattern (void)
{
    for (unsigned x = lb_x; x < ub_x; ++x)
        for (unsigned y = lb_y; y < ub_y; ++y)
            set(x,y,rs.ProbabilisticChoice(random_bias));
    return;
}
void Universe::random_pattern (float bias)
{
    for (unsigned x = lb_x; x < ub_x; ++x)
        for (unsigned y = lb_y; y < ub_y; ++y)
            set(x,y,rs.ProbabilisticChoice(bias));
    return;
}
// bounded pattern
void Universe::random_pattern (float bias,int x0,int y0,int dx,int dy)
{
    for (int x = x0; x <= x0+dx; ++x) {
        for (int y = y0; y <= y0+dy; ++y)
        set(x,y,rs.ProbabilisticChoice(bias));
    }
    return;
}


// empty universe
void Universe::empty (void)
{
    for (auto ity = space->begin(); ity != space->end(); ++ity)
        for (auto itx = ity->begin(); itx != ity->end(); ++itx)
            *itx = 0;
    return;
}

// rotate grid in 90 degree intervals counter-clockwise
void Universe::rotate(short n)
{
    n = n%4;
    if (n == 0) return;
    
    deque<bool> tmp;
    unique_ptr<deque<deque<bool>>> buffer(new deque<deque<bool>>);
    if (n % 2 != 0) {
        tmp.resize(y_size,0);
        buffer->resize(x_size,std::move(tmp));
    } else {
        tmp.resize(x_size,0);
        buffer->resize(y_size,std::move(tmp));
    }
    
    switch (n)
    {
    case 1:
        for (int y = 0; y < x_size; ++y)
            for (int x = 0; x < y_size; ++x)
                (*buffer)[y][x] = (*space)[y_size-1-x][y];
        break;
    case 2:
        for (int y = 0; y < y_size; ++y)
            for (int x = 0; x < x_size; ++x)
                (*buffer)[y][x] = (*space)[y_size-1-y][x_size-1-x];
        break;
    case 3:
        for (int y = 0; y < x_size; ++y)
            for (int x = 0; x < y_size; ++x)
                (*buffer)[y][x] = (*space)[x][x_size-1-y];
        break;
    default:
        cerr << "invalid rotation number: " << n << '\n';
        exit(1);
        break;
    }
    
    if (n % 2 != 0) {
        x_size = buffer->front().size();
        y_size = buffer->size();
        int tmp = lb_x; lb_x = lb_y;
        lb_y = tmp;
        tmp = ub_x; ub_x = ub_y;
        ub_y = tmp;
    }
    space = std::move(buffer);
    return;
}

// flip grid over x-axis
void Universe::flip (void)
{
    deque<bool> tmp(x_size,0);
    for (int y = 0; y < y_size / 2; ++y) {
        tmp = (*space)[ub_y-1-y];
        (*space)[ub_y-1-y] = (*space)[y];
        (*space)[y] = tmp;
    }
    return;
}


// load pattern from stream
void Universe::load_pattern (ifstream &is)
{
    unsigned x_bound, y_bound;
    is >> x_bound;
    is >> y_bound;
    unsigned start_y = (pattern_load_point)[1] + y_bound;  // point furthest up
    unsigned end_x = (pattern_load_point)[0] + x_bound;    // point furthest right
    if (end_x > ub_x || start_y > ub_y)
        {cerr << "loaded pattern is too large or does not fit into space at given point\n"; exit(1);}
    for (unsigned y = start_y-1; y >= (pattern_load_point)[1]; --y) {
        for (unsigned x = (pattern_load_point)[0]; x < end_x; ++x) {
            short val; is >> val;
            if (val != 2) set(x,y,val);
        }
    }
    is.close();
    return;
}

// load pattern from stream at given point
void Universe::load_pattern_at (ifstream &is, int at_x, int at_y,bool crop = 0)
{
    int x_bound, y_bound;
    is >> x_bound;
    is >> y_bound;
    int start_y = at_y + y_bound;   // point furthest up
    int end_x = at_x + x_bound;     // point furthest right
    if (!crop && (end_x > ub_x || start_y > ub_y || at_x < lb_x || at_y < lb_y)) {
        cerr << "loaded pattern is too large or does not fit into space at given point: ("
            << at_x << ',' << at_y << ")\n";
        exit(1);
    }
    for (int y = start_y-1; y >= at_y; --y) {
        for (int x = at_x; x < end_x; ++x) {
            short val; is >> val;
            if (crop && (y>ub_y || y<lb_y || x>ub_x || x<lb_x)) continue;
            if (val != 2) set(x,y,val);
        }
    }
    is.close();
    return;
}

// load pattern from stream
void Universe::load_pattern_toroidal (ifstream& is,int at_x,int at_y)
{
    if (at_x < lb_x || at_x > ub_x || at_y < lb_y || at_y > ub_y) {
        cerr << "given coordinates out of bounds: (" << at_x << ',' << at_y << ")\n";
        exit(1);
    }
    at_x = at_x-lb_x;
    at_y = at_y-lb_y;
    int x_bound, y_bound;
    is >> x_bound; is >> y_bound;
    for (int y_n = y_bound-1; y_n >= 0; --y_n) {
        for (int x_n = 0; x_n < x_bound; ++x_n) {
            int x = (at_x + x_n)%x_size,
                y = (at_y + y_n)%y_size;
            short val; is >> val;
            if (val != 2) set(x,y,val);
        }
    }
    is.close();
    return;
}

// save pattern to file
deque<bool> get_vertical (const deque<deque<bool>> &u, int &lb, int &ub, int &x)
{
    // lb,ub inclusive
    deque<bool> line(ub-lb+1,0);
    for (int i = lb; i <= ub; ++i)
        line[i-lb] = u[i][x];
    return line;
}
void Universe::save_pattern (ofstream &os)
{
    // find subset of non-empty space
    int lby,uby,lbx,ubx; // inclusive
    // vertical bounds
    deque<bool> empty_xline(x_size,0);
    for (int i = 0; i < y_size; ++i) {
        if ((*space)[i] == empty_xline) lby = i+1;
        else break;
    }
    for (int i = y_size-1; i >= 0; --i) {
        if ((*space)[i] == empty_xline) uby = i-1;
        else break;
    }
    // horizontal bounds
    deque<bool> empty_yline(uby-lby+1,0);
    for (int i = 0; i < x_size; ++i) {
        if (empty_yline == get_vertical(*space,lby,uby,i)) lbx = i+1;
        else break;
    }
    for (int i = ub_x-1; i >= 0; --i) {
        if (empty_yline == get_vertical(*space,lby,uby,i)) ubx = i-1;
        else break;
    }
    // write dimensions to file
    os << ubx-lbx+1 << ' ' << uby-lby+1 << '\n';
    // write pattern to file
    for (int y = uby; y > lby; --y) {
        for (int x = lbx; x < ubx; ++x)
            os << (*space)[y][x] << ' ';
        os << (*space)[y][ubx] << '\n';
    }
    for (int x = lbx; x < ubx; ++x)
        os << (*space)[lby][x] << ' ';
    os << (*space)[lby][ubx];
    os.close();
    return;
}

void Universe::save_pattern(ofstream &os, int lby, int uby, int lbx, int ubx)
{
    // inclusive lb,ub
    // write dimensions to file
    os << ubx-lbx+1 << ' ' << uby-lby+1 << '\n';
    // write pattern to file
    for (int y = uby; y > lby; --y) {
        for (int x = lbx; x < ubx; ++x)
            os << (*space)[y][x] << ' ';
        os << (*space)[y][ubx] << '\n';
    }
    for (int x = lbx; x < ubx; ++x)
        os << (*space)[lby][x] << ' ';
    os << (*space)[lby][ubx];
    os.close();
    return;
}


// determine if space is empty
bool Universe::emptyQ(void)
{
    for (int y = lb_y; y < ub_y; ++y) {
        for (int x = lb_x; x < ub_x; ++x)
            if (this->operator()(x,y)) return 0;
    }
    return 1;
}

// search space for given pattern
unsigned Universe::find_pattern(string filename)
{
    ifstream is(filename);
    int x_lim, y_lim;
    is >> x_lim;
    is >> y_lim;
    is.close();
    int x_bound = (int)x_size - x_lim;
    int y_bound = (int)y_size - y_lim;
    unsigned instances = 0;
    Universe u(x_size,y_size);
    for (int x = 0; x < x_bound; ++x) {
        for (int y = 0; y < y_bound; ++y) {
            u.copy_space(*this);
            u.load_pattern_at(filename,x,y);
            if (u.get_space() == *space) ++instances;
        }
    }
    return instances;
}

// only finds first instance
vector<int> Universe::find_pattern_loc(string filename,vector<int>& null_lp)
{
    ifstream is("Patterns/"+filename+".dat");
    int x_lim, y_lim;
    is >> x_lim;
    is >> y_lim;
    is.close();
    int x_bound = (int)x_size - x_lim + 1;
    int y_bound = (int)y_size - y_lim + 1;
    Universe u(x_size,y_size);
    vector<int> idx(2,0);
    { // check null first
        u.copy_space(*this);
        u.load_pattern_at(filename,null_lp[0],null_lp[1]);
        if (u.get_space() == *space) return null_lp;
    }
    for (int x = 0; x < x_bound; ++x) {
        for (int y = 0; y < y_bound; ++y) {
            u.copy_space(*this);
            u.load_pattern_at(filename,x,y);
            if (u.get_space() == *space) {
                idx = {x+lb_x,y+lb_y};
                return idx;
            }
        }
    }
    idx = {-1,-1};
    return idx;
}

unsigned Universe::find_cycle(string filename, unsigned search_lim)
{
    unsigned steps = 0;
    unique_ptr<deque<deque<bool>>> buffer(new deque<deque<bool>>);
    *buffer = *space;
    empty();
    load_pattern(filename);
    for (unsigned i = 0; i < search_lim; ++i, ++steps) {
        update();
        if (find_pattern(filename)) break;
    }
    space = std::move(buffer);
    return (++steps)*(!(steps==search_lim));
}


void Universe::load_pattern_vector(string filename, int at_x, int at_y, vector<bool>& v)
{
    std::ifstream is(filename);
    int x_bound, y_bound;
    is >> x_bound;
    is >> y_bound;
    int start_y = at_y + y_bound;   // point furthest up
    int end_x = at_x + x_bound;     // point furthest right
    vector<bool>::iterator it = v.begin();
    if (end_x > ub_x || start_y > ub_y)
        {cerr << "loaded pattern is too large or does not fit into space at given point\n"; exit(1);}
    for (int y = start_y-1; y >= at_y; --y) {
        for (int x = at_x; x < end_x; ++x) {
            short val; is >> val;
            if (val != 2) set(x,y,*it);
            ++it;
        }
    }
    is.close();
    return;
}

// load an initial configuration to shape into a bug
double pNorm (double p,int x,int y) { return pow(pow(abs(x),p) + pow(abs(y),p),1/p); }
void Universe::load_bug_starter (int r1,int r2,int x,int y)
{
    unsigned center_gap = rule_set->range / 4;
    for (int cx = -r2; cx <= r2; ++cx) {
        for (int cy = -r2; cy <= r2; ++cy) {
            double dist1 = pNorm(2,cx,cy),
                   dist2 = pNorm(2,cx,cy+center_gap);
            if (dist1 > r2) continue;
            else if (dist2 <= r1) set(x+cx,y+cy,0);
            else set(x+cx,y+cy,1);
        }
    }
    return;
}

// calculate average position of on-cell in neighborhood
vector<double> Universe::density_vector (int x,int y)
{
    vector<double> n({0,0});
    unsigned N = 0;
    for (int dx = -rule_set->range; dx <= rule_set->range; ++dx) {
        for (int dy = -rule_set->range; dy <= rule_set->range; ++dy) {
            if (this->operator()(x,y) == 0 || (!dx && !dy)) continue;
            n[0] += dx; n[1] += dy;
            ++N;
        }
    }
    if (N == 0) return n;
    n[0] /= N; n[1] /= N;
    return n;
}

// calculate perturbation magnitude to change process product
unsigned Universe::density_threshold (int x,int y)
{
    unsigned N = count_neighborhood(x,y);
    if (this->operator()(x,y)) {
        if (N < rule_set->survive_min)
            return rule_set->survive_min - N;
        else if (N <= rule_set->survive_max)
            return rule_set->survive_max - N;
    } else {
        if (N < rule_set->birth_min)
            return rule_set->birth_min - N;
        else if (N <= rule_set->birth_max)
            return rule_set->birth_max - N;
    } return 0;
}


void operator++(vector<bool>& v)
{
    for (auto rit = v.rbegin(); rit != v.rend(); ++rit) {
        if (!*rit) {*rit = 1; break;}
        else *rit = 0;
    }
    return;
}
bool last_binQ(vector<bool>& v)
{
    for (auto bit : v) if (!bit) return 0;
    return 1;
}
std::ostream& operator<<(std::ostream& os,vector<bool>& v)
{
    for (auto it = v.cbegin(); it != --v.cend(); ++it)
        os << *it << ' ';
    os << *(--v.cend());
    return os;
}