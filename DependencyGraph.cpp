#include "DependencyGraph.h"

/* ----------------------------- */
/* Symmetry Generation Functions */
/* ----------------------------- */

ostream& print_space_with_dont_cares (ostream& os,Universe& u)
{
    for (int y = u.get_ybound()-1; y >= 0; --y) {
        for (int x = 0; x < u.get_xbound()-1; ++x) {
            if (!u.count_neighborhood(x,y))
                os << 2 << ' ';
            else os << u(x,y) << ' ';
        }
        if (!u.count_neighborhood(u.get_xbound()-1,y))
            os << 2 << '\n';
        else os << u(u.get_xbound()-1,y) << '\n';
    }
    return os;
}

void generate_symmetries(vector<string>& patterns,rule& rl)
{
    list<list<Universe>> geometries;
    unsigned x_size, y_size;

    for (string name : patterns) {
        ifstream is("Patterns/"+name+".dat");
        is >> x_size; is >> y_size;
        is.close();
        list<Universe> tmp_l;
        tmp_l.emplace_back(x_size,y_size,rl);
        tmp_l.back().load_pattern_at(name,0,0);
        geometries.push_back(std::move(tmp_l));
    }

    // generate all symmetries
    for (list<Universe>& u_l : geometries) {
        u_l.emplace_back(u_l.front());
        u_l.back().flip();
        for (short n = 1; n <= 3; ++n) {
            u_l.emplace_back(u_l.front());
            u_l.back().rotate(n);
            u_l.emplace_back(u_l.back());
            u_l.back().flip();
        }
    }
    
    // remove redundant symmetries
    for (auto& sym : geometries) {
        for (auto g1 = sym.begin(); g1 != sym.end(); ++g1) {
            for (auto g2 = sym.begin(); g2 != sym.end(); ++g2) {
                if (g1 == g2) continue;
                else if (g1->get_space() == g2->get_space()) {
                    g2 = sym.erase(g2);
                    --g2;
                }
            }
        }
    }

    // save geometries to file
    list<list<Universe>>::iterator it = geometries.begin();
    for (string name : patterns) {
        unsigned index = 1;
        string path = "Patterns/" + name;
        for (auto& sym : *it) {
            ofstream os(path+to_string(index)+".dat");
            os << sym.get_xbound() << ' ';
            os << sym.get_ybound() << '\n';
            print_space_with_dont_cares(os,sym);
            // os << sym.get_space();
            os.close();
            ++index;
        }
        ++it;
    }

    return;
}


/* -------------------------- */
/* Process Analysis Functions */
/* -------------------------- */

vector<short> process_from_space (int x,int y,Universe& u)
{
    vector<short> neighbors(8,0); // assume GoL rule
    short n = 0;
    for (int dy = u.get_rules()->range; dy >= -u.get_rules()->range; --dy) {
        for (int dx = -u.get_rules()->range; dx <= u.get_rules()->range; ++dx) {
            if (!dx && !dy) continue;
            int real_x = (u.get_xbound()+x+dx)%u.get_xbound(),
                real_y = (u.get_ybound()+y+dy)%u.get_ybound();
            if (u.count_neighborhood(real_x,real_y) == 0)
                neighbors[n] = 2;
            else neighbors[n] = u(real_x,real_y);
            ++n;
        }
    }
    return neighbors;
}

short type_from_state (bool state,vector<short>& neighbors)
{
    // assumes GoL rule
    unsigned sum = 0;
    if (!state) {
        for (short n = 0; n < 8; ++n) {
            if (neighbors[n] == 2) return 0;
            sum += neighbors[n];
        } return (sum == 3)?1:2;
    } else {
        for (short n = 0; n < 8; ++n)
            sum += neighbors[n];
        return (sum==2 || sum==3)?5:4;
    }

    cerr << "invalid state (" << state << ") provided\n";
    exit(1);
}

ofstream& save_non_null_processes (string filename,string group,unsigned& index,Universe& u)
{
    ofstream os(filename,ofstream::app);
    ofstream os1("process_sets/"+group+".dat");
    for (int y = u.get_ybound()-1; y >= 0; --y) {
        for (int x = 0; x < u.get_xbound(); ++x) {
            if (u.count_neighborhood(x,y) == 0) continue;
            vector<short> neighbors = process_from_space(x,y,u);
            os << ++index << ' ' << type_from_state(u(x,y),neighbors) << ' '
                << neighbors << '\n';
            os1 << index << ' ';
        }
    }
    os.close(); os1.close();
}

// takes list of base names and corresponding indicies
void generate_process_set (vector<pair<string,unsigned>>& names,rule& rl,string set_name)
{
    ofstream os("process_sets/glider_process_set.dat"); os.close();
    unsigned index = 0;
    for (auto& g : names) {
        for (unsigned idx = 1; idx <= g.second; ++idx) {
            ifstream is(to_filename(g.first,idx));
            unsigned x_sz, y_sz;
            is >> x_sz; is >> y_sz;
            is.close();
            Universe u(x_sz+2,y_sz+2,rl);
            u.load_pattern_at(g.first+to_string(idx),1,1);
            string filename = "process_sets/"+set_name+"_process_set.dat";
            save_non_null_processes(filename,g.first+to_string(idx),index,u);
        }
    }
    return;
}


/* ----------------- */
/* Transition Search */
/* ----------------- */

void load_perturbation (Universe& u,unsigned orientation_num,string& p_class)
{
    u.load_pattern_at(p_class,0,0);
    
    u.rotate((orientation_num-1)/2);
    if (orientation_num%2 == 0) u.flip();
    
    return;
}

void find_transitions(string name,vector<pair<string,unsigned>>& geometries,vector<string>& classes,rule& rl)
{
    unsigned total_geometries = 0;
    for (auto& set : geometries)
        total_geometries += set.second;
    vector<string> names(total_geometries); {
        unsigned idx = 0;
        for (auto& set : geometries) {
            for (int n = 1; n <= set.second; ++n) {
                names[idx] = set.first + to_string(n);
                ++idx;
    }}}
    Universe u1(7,7,rl);

    ofstream os("DG_files/"+name+"_transitions.dat");
    
    for (auto g1 = names.begin(); g1 != names.end(); ++g1) {
        for (auto g2 = names.begin(); g2 != names.end(); ++g2) {
            for (auto& p_class : classes) {
                
                /* glider-only condition */
                if ((*g1)[1] == 'r' && (p_class != "grey" && p_class != "green"))
                    continue;
                else if ((*g1)[1] == 'w' && (p_class == "grey" || p_class == "green"))
                    continue;

                u1.empty();
                load_perturbation(u1,(unsigned)g1->back(),p_class);
                u1.load_pattern_at(*g1,1,1);
                u1.update();
                vector<int> loc({1,1});
                loc = u1.find_pattern_loc(*g2,loc);
                if (*g1 == "gw5" && *g2 == "gr6 ") cout << p_class << ' ' << loc << '\n';
                if (loc[0] == -1) continue;
                os << *g1 << ' ' << *g2 << ' '
                    << (int)loc[0]-1 << ' ' << (int)loc[1]-1 << '\n';
            }
        }
    }
    os.close();
    return;
}


/* ----------------------------- */
/* Dependency Analysis Functions */
/* ----------------------------- */

unsigned index_lookup (vector<short>& process /* type-first form */,string set_name,
                        unsigned idx_lb,unsigned idx_ub,set<unsigned>& used)
{
    ifstream is("process_sets/"+set_name+"_process_set.dat");
    string line;
    vector<short> holder(9,0); unsigned index;
    while (getline(is,line)) {
        istringstream iss(line);
        iss >> index;
            if (index < idx_lb) continue;
            else if (index > idx_ub) {index = 0; break;}
        iss >> holder;
        if (process != holder) continue;
        else if (std::binary_search(used.begin(),used.end(),index))
            continue;
        break;
    }
    is.close();
    if (index == 0) {
        cerr << "could not find process in given bounds: [" << idx_lb << ',' << idx_ub << "]\n";
        exit(1);
    }
    used.insert(index);
    return index;
}
unsigned index_lookup (vector<short>& process /* type-first form */,string set_name,
                        unsigned idx_lb,unsigned idx_ub,
                        vector<vector<unsigned>>& used,unsigned n,unsigned p1)
{
    ifstream is("process_sets/"+set_name+"_process_set.dat");
    string line;
    vector<short> holder(9,0); unsigned index;
    while (getline(is,line)) {
        istringstream iss(line);
        iss >> index;
            if (index < idx_lb) continue;
            else if (index > idx_ub) {index = 0; break;}
        iss >> holder;
        if (process != holder) continue;
        bool cont = 0;
        for (auto& p : used)
            if (index == p[0] && (p1 == p[1] || n == p[2]))
                {cont = 1; break;}
        if (cont) continue;
        break;
    }
    is.close();
    if (index == 0) {
        cerr << "could not find process in given bounds: [" << idx_lb << ',' << idx_ub << "]\n";
        exit(1);
    }
    vector<unsigned> tmp({index,p1,n});
    used.push_back(std::move(tmp));
    return index;
}

// returns type-first process from GoL universe
vector<short> cell_to_GoL_process (Universe& u,int x,int y)
{
    vector<short> proc(9,0);
    // get neighbor cells
    vector<short>::iterator cell = ++proc.begin();
    unsigned neighborhood_sum = 0;
    for (int dy = 1; dy >= -1; --dy) {
        for (int dx = -1; dx <= 1;  ++dx) {
            int real_x = (u.get_xbound()+x+dx)%u.get_xbound();
            int real_y = (u.get_ybound()+y+dy)%u.get_ybound();
            if (dx == 0 && dy == 0) {proc[0] = u(x,y); continue;}
            else if (!u.count_neighborhood(real_x,real_y)) *cell = 2;
            else {
                *cell = u(real_x,real_y);
                neighborhood_sum += u(real_x,real_y);
            }
            ++cell;
        }
    }
    // determine process type
    for (int i = 1;  i < 9; ++i)
        if (proc[i] == 2) return proc;
    if (!proc[0]) {
        if (neighborhood_sum == 3) {proc[0] = 1; return proc;}
        else {proc[0] = 2; return proc;}
    } else {
        if (neighborhood_sum < 2 || neighborhood_sum > 3) {proc[0] = 4; return proc;}
        else {proc[0] = 5; return proc;}
    } return proc;
}

bool sort_vec (vector<unsigned>& v1,vector<unsigned>& v2)
{
    for (std::size_t i = 0; i < v1.size();  ++i) {
        if (v1[i] < v2[i]) return 1;
        else if (v1[i] > v2[i]) return 0;
    }
    return 1;
}
void sort_links (string file)
{
    list<vector<unsigned>> links;
    vector<unsigned> holder(4);
    ifstream is("DG_files/"+file+"_links.dat"); string line;
    while (getline(is,line)) {
        istringstream iss(line);
        iss >> holder;
        links.push_back(holder);
    }
    is.close();
    links.sort(sort_vec);
    ofstream os("DG_files/"+file+"_links.dat");
    for (auto& l : links)
        os << l << '\n';
    return;
}

void generate_links (string name,vector<pair<string,unsigned>>& sets,rule& rl)
{
    /* copy process set to DG_files folder */
    ifstream is_cp("process_sets/"+name+"_process_set.dat");
    ofstream os_cp("DG_files/"+name+"_process_set.dat");
    string line;
    while (getline(is_cp,line)) {
        os_cp << line << '\n';
    } is_cp.close(); os_cp.close();
    
    /* load geometry grouped sets */
    map<string,vector<unsigned>> process_indicies; {
        for (auto& config : sets) {
            for (unsigned idx = 1; idx <= config.second; ++idx) {
                ifstream is("process_sets/"+config.first+to_string(idx)+".dat");
                vector<unsigned> tmp(2,0);
                unsigned val; is >> val;
                tmp[0] = val;
                while (!is.eof()) is >> val;
                tmp[1] = val;
                is.close();
                process_indicies.emplace(config.first+to_string(idx),std::move(tmp));
            }
        }
    };

    unsigned x_size, y_size;
    ifstream is("DG_files/"+name+"_transitions.dat"); {
        /* get universe space size */
        string tmp; is >> tmp;
        ifstream is1("Patterns/"+tmp+".dat");
        is1 >> x_size; is1 >> y_size;
        x_size += 2; y_size += 2;
        is.close(); is1.close();
    };

    string g1,g2; int relx,rely;
    vector<short> proc_holder(9,0);
    unsigned p1, p2;
    ofstream os("DG_files/"+name+"_links.dat");
    is.open("DG_files/"+name+"_transitions.dat");
    while (getline(is,line)) {
        istringstream iss(line);
        // set up spaces
        Universe u1(x_size,y_size,rl),u2(x_size,y_size,rl);
        iss >> g1; iss >> g2;
        iss >> relx; iss >> rely;
        u1.load_pattern_at(g1,1,1);
        u2.load_pattern_at(g2,1+relx,1+rely);

        // find dependencies
        set<unsigned> used_g1;
        vector<vector<unsigned>> used_g2;
        for (int y = y_size-2; y > 0; --y) {
            for (int x = 1; x < x_size-1; ++x) {
                if (u1.count_neighborhood(x,y) == 0) continue;
                proc_holder = cell_to_GoL_process(u1,x,y);
                p1 = index_lookup(proc_holder,name,
                    process_indicies[g1][0],process_indicies[g1][1],
                    used_g1);
                unsigned n = 0;
                for (int dy = 1; dy >= -1; --dy) {
                    for (int dx = -1; dx <= 1; ++dx) { ++n;
                        if (u2.count_neighborhood(x+dx,y+dy) == 0) continue;
                        proc_holder = cell_to_GoL_process(u2,x+dx,y+dy);
                        p2 = index_lookup(proc_holder,name,
                            process_indicies[g2][0],process_indicies[g2][1],
                            used_g2,n,p1);
                        os << p1 << ' ' << p2 << ' ';
                        if (u2.count_neighborhood(x,y) == 0)
                            os << 2 << ' ';
                        else os << u2(x,y) << ' ';
                        os << 10-n << '\n';
                    }
                }
            }
        }
    }
    is.close(); os.close();

    sort_links(name);
    return;
}


void rotate (vector<short>& p)
{
    short tmp = p[3];
    p[3] = p[1];
    p[1] = p[6];
    p[6] = p[8];
    p[8] = tmp;
    //
    tmp = p[5];
    p[5] = p[2];
    p[2] = p[4];
    p[4] = p[7];
    p[7] = tmp;
    return;
}
void flip (vector<short>& p)
{
    short tmp;
    for (short i = 1; i < 4; ++i) {
        tmp = p[i];
        p[i] = p[i+5];
        p[i+5] = tmp;
    }
    return;
}
bool check_link_symmetry (string obj)
{
    vector<short> temp_proc(9);
    vector<vector<short>> process_set(352,temp_proc);
    std::ifstream is("DG_files/"+obj+"_process_set.dat");
    string line;
    unsigned idx = 0;
    while (getline(is,line)) {
        std::istringstream iss(line);
        is >> temp_proc;
        process_set[idx] = temp_proc;
        ++idx;
    } is.close();

    idx = 0;
    
}