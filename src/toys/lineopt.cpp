#include <iostream>
#include <2geom/path.h>
#include <2geom/pathvector.h>
#include <2geom/svg-path-parser.h>

#include "lineoptimization.h"

int main(int argc, char ** argv) {
    srand(time(NULL));

    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <line data svgd file>" << std::endl;
        return 0;
    }

    std::string const lines_file = argv[1];

    Geom::PathVector lines = Geom::read_svgd(lines_file.c_str());

    std::cout << "file has " << lines.size() << " lines" << std::endl;

    for (auto const& curve : lines) {
        std::cout << curve.size() << " ";
    }
    std::cout << std::endl;

    LineOptimization opt;
    opt.addPaths(lines);
    opt.splitAtForks();

    LineOptimization::Solution sol = opt.greedySolution();
    LineOptimization::Solution sol2 = opt.randomGreedySolution();

    return 0;
}
