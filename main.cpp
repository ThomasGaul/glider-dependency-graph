#include <iostream>
#include <vector>
#include <memory>
#include <string>

#include "LtL.h"
#include "DependencyGraph.h"
#include "Geometry.h"


int main (int argc,char* argv[])
{
    /* -------------- */
    /* Get Symmetries */
    /* -------------- */

    rule GoL((string)"GoL");
    Universe u(7,7,GoL);
    vector<string> base({"gw","gr"});
    generate_symmetries(base,GoL);


    /* --------------------- */
    /* Get Interaction Graph */
    /* --------------------- */

    vector<pair<string,unsigned>> names;
        names.emplace_back("gw",8);
        names.emplace_back("gr",8);
    vector<string> classes({"grey","black","orange","brown","blue","green"});
    find_transitions("glider",names,classes,GoL);


    /* ------------------------ */
    /* Extract Dependency Graph */
    /* ------------------------ */

    generate_process_set(names,GoL,"glider");
    generate_links("glider",names,GoL);
    sort_links("glider");


    // optional
    /* ----------------------- */
    /* Reduce Dependency Graph */
    /* ----------------------- */

    // load in graph
    list<vector<unsigned>> reduced_graph;
    vector<unsigned> hold(4);
    vector<unsigned> counts(352,0);
    
    string line;
    ifstream is("DG_files/glider_links.dat");
    while (getline(is,line)) {
        istringstream iss(line);
        iss >> hold;
        if (hold[2] == 2) continue;
        ++counts[hold[0]-1];
        reduced_graph.push_back(hold);
    }
    is.close();
    for (auto it = reduced_graph.begin(); it != reduced_graph.end(); ++it) {
        if (counts[it->operator[](1)-1] == 0) {
            it = reduced_graph.erase(it);
            --it;
        }
    }

    // overwrite with reduced graph
    ofstream os("DG_files/glider_links.dat");
    for (auto l : reduced_graph)
        os << l << '\n';
    os.close();


    /* -------------------- */
    /* Reconstruct Geometry */
    /* -------------------- */

    layout lo("glider");
    lo.initialize(2 /* can by any process in graph */);
    lo.find();
    lo.bash_display();

    return 0;
}