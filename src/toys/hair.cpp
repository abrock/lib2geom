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

#include "geom-pathstroke.h"

#include "RunningStats.h"

using namespace Geom;


#include "hair.h"

int main(int argc, char **argv) {
    if (argc > 1) {
        std::cout << "Argument 1: " << argv[1] << std::endl;
    }
    PathVector boundary = read_svgd("svgd/tail-outline-1.svgd");
    PathVector curve = read_svgd("svgd/tail-curve.svgd");

    Hair hair(boundary[0], curve[0]);


    hair.run2();
    hair.write("tail-1-results.svg");
    hair.writeStitches("tail-1-results.txt");

    hair.writeAreas("tail-1-areas.svg");

    return 0;
}
