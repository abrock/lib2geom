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

#include "geom-pathstroke.cpp"

#include "RunningStats.h"

using namespace Geom;

int main(int argc, char **argv) {
    if (argc > 1) {
        std::cout << "Argument 1: " << argv[1] << std::endl;
        PathVector path = read_svgd(argv[1]);
        for (Path curve1 : path) {
            for (Path curve2 : path) {
                auto intersections = curve1.intersect(curve2);
            }
        }
    }

    return 0;
}
