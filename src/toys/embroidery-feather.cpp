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

std::vector<Geom::Point> makeCircleN(double const radius, size_t const n) {
    std::vector<Geom::Point> result;
    for (size_t ii = 0; ii <= n; ++ii) {
        double const angle = 2. * M_PI * static_cast<double>(ii) / n;
        double const x = radius * std::cos(angle);
        double const y = radius * std::sin(angle);
        result.push_back(Geom::Point(x,y));
    }
    return result;
}

std::vector<Geom::Point> makeCircleLength(double const radius, double const length) {
    size_t const n = static_cast<size_t>(std::floor(radius * 2. * M_PI / length));
    return makeCircleN(radius, n);
}

std::vector<Geom::Point> makeArchimedesSpiral(
        double const radius,
        double const length,
        double const rotation_distance) {
    // r = alpha * phi
    double const alpha = rotation_distance / (2. * M_PI);
    double phi = 0;
    std::vector<Geom::Point> result;
    result.push_back(Geom::Point(0,0));
    while (alpha * phi < radius) {
        Geom::Point next_point (alpha * phi * cos(phi), alpha * phi * sin(phi));
        while (Geom::distance(next_point, result.back()) < length) {
            phi += 1e-6;
            next_point = Geom::Point(alpha * phi * cos(phi), alpha * phi * sin(phi));
        }
        result.push_back(next_point);
    }
    return result;
}

std::vector<Geom::Point> makeArchimedesSpiralIncreasingLength(
        double const radius,
        double const max_length,
        double const rotation_distance) {
    // r = alpha * phi
    double const alpha = rotation_distance / (2. * M_PI);
    double const circumference = radius * 2. * M_PI;
    double const d_phi = max_length / radius;
    double phi = 0;
    std::vector<Geom::Point> result;
    result.push_back(Geom::Point(0,0));
    while (alpha * phi < radius) {
        phi += d_phi;
        Geom::Point next_point (alpha * phi * cos(phi), alpha * phi * sin(phi));
        result.push_back(next_point);
    }
    return result;
}

namespace fs = boost::filesystem;

int main(int argc, char **argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " outlinefile curvefile templatefile" << std::endl;
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
    const std::string template_file(argv[3]);
    std::cout << "Outline file: " << outline_file << std::endl
              << "Curve file: " << curve_file << std::endl;
    PathVector boundary = read_svgd(outline_file.c_str());
    PathVector curve = read_svgd(curve_file.c_str());
    PathVector template_path = read_svgd(template_file.c_str());

    auto spiral2 = makeArchimedesSpiralIncreasingLength(90, 6, 10);
    Hair::addStartStop(spiral2);
    Hair::writeStitches(spiral2, "spiral_inc.txt");

    auto spiral = makeArchimedesSpiral(90, 3.6, 10);

    Hair::addStartStop(spiral);
    Hair::writeStitches(spiral, "spiral.txt");

    Hair hair;

    hair.setOutline(boundary[0]);
    hair.setFeatherCurve(curve[0]);
    hair.setFeatherTemplate(template_path[0]);

    hair.runFeather();

    hair.write(outline_file + ".svg");
    hair.writeCurves(outline_file + "-curves.svg");
    hair.writeStitches(hair.greedySolution, outline_file + "-greedy.txt");

    hair.printStats();

    std::cout << "Boundary: " << boundary.size() << std::endl;

    std::vector<std::vector<Geom::Point> > total_vector;
    for (size_t ii = 0; ii < boundary.size() && ii < curve.size() && ii < template_path.size(); ++ii) {
        Hair hair;
        hair.setDensity(3.0);

        hair.setOutline(boundary[ii]);
        hair.setFeatherCurve(curve[ii]);
        hair.setFeatherTemplate(template_path[ii]);

        hair.runFeather();
        total_vector.push_back(hair.greedySolution);
    }
    hair.writeStitches(total_vector, outline_file + "-all-greedy.txt");



    /*
    hair.writeAreas("tail-1-areas.svg");

    hair.writeForwardAreas("tail-1-forward-areas.svg");
    hair.writeReverseAreas("tail-1-reverse-areas.svg");
    */


    return 0;
}
