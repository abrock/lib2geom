#include <2geom/d2.h>
#include <2geom/intersection-graph.h>
#include <2geom/path.h>
#include <2geom/sbasis.h>
#include <2geom/svg-path-parser.h>
#include <2geom/transforms.h>
#include <2geom/sbasis-to-bezier.h>
#include <2geom/svg-path-writer.h>
#include <2geom/curve.h>

//#include <toys/path-cairo.h>
//#include <toys/toy-framework-2.h>

//#include <cairo/cairo.h>
//#include <gtk-2.0/gdk/gdkevents.h>

#include <cstdlib>
#include <boost/algorithm/minmax_element.hpp>

#include <iostream>
#include <fstream>
#include <random>
#include <boost/filesystem.hpp>

#include "helper/geom-pathstroke.h"

#include "RunningStats.h"

using namespace Geom;


#include "hair.h"

namespace fs = boost::filesystem;

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " outlinefile curvefile" << std::endl;
        return 0;
    }
    if (argc > 1) {
        std::cerr << "Argument 1: " << argv[1] << std::endl;
    }

    fs::path full_path( fs::initial_path<fs::path>() );
    std::cerr << "Current path: " << fs::current_path() << std::endl;
    full_path = fs::system_complete( fs::path( argv[0] ) );
    std::cerr << full_path << std::endl;
    //Without file name
    std::cerr << full_path.stem() << std::endl;
    fs::current_path(full_path.remove_filename());


    const std::string outline_file(argv[1]);
    const std::string curve_file(argv[2]);
    std::cout << "Outline file: " << outline_file << std::endl
              << "Curve file: " << curve_file << std::endl;
    PathVector boundary = read_svgd(outline_file.c_str());
    PathVector curve = read_svgd(curve_file.c_str());

/*
    Hair hair(boundary[0], curve[0]);


    hair.run();
    hair.write(outline_file + ".svg");
    hair.writeCurves(outline_file + "-curves.svg");
    hair.writeStitches(hair.greedySolution, outline_file + "-greedy.txt");

    hair.printStats();

    */

    std::vector<std::string> all_stats;

    std::vector<std::vector<Geom::Point> > total_vector;
    for (size_t ii = 0; ii < boundary.size() && ii < curve.size(); ++ii) {
        std::cout << "Running curve #" << ii+1 << " out of " << std::min(boundary.size(), curve.size()) << std::endl;
        Hair hair;
        hair.setDensity(3.05);

        hair.setOutline(boundary[ii]);
        hair.setCurve(curve[ii]);

        hair.run();
        all_stats.push_back(hair.getStats());
        total_vector.push_back(hair.greedySolution);
    }
    Hair hair;
    hair.writeStitches(total_vector, outline_file + "-all-greedy.txt");
    std::cout << "All stitch length stats:" << std::endl;
    for (auto const & it : all_stats) {
        std::cout << it << std::endl;
    }

    /*
    hair.writeAreas("tail-1-areas.svg");

    hair.writeForwardAreas("tail-1-forward-areas.svg");
    hair.writeReverseAreas("tail-1-reverse-areas.svg");
    */


    return 0;
}
