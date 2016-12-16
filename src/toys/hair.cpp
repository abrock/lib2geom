#include "hair.h"
#include <cmath>

#include "embroideryoptimization.h"

void Hair::setOutline(Path const& p) {
    _outline = p;
}

void Hair::setFeatherCurve(Path const& p) {
    _feather_curve = p;
}

void Hair::setFeatherTemplate(Path const& p) {
    _feather_template = p;
}


PathTime plus(PathTime t, double dt) {
    t.t += dt;
    while (t.t >= 1.0) {
        t.t -= 1;
        t.curve_index += 1;
    }
    return t;
}

PathTime plus(PathTime t, double dt, Path const& curve) {
    t = plus(t, dt);
    if (t.curve_index >= curve.size()) {
        t.curve_index--;
        t.t = 1.0;
    }
    return t;
}


/**
 * @brief writeStitches Write the stitches into a text file suitable for libembroidery-convert
 * @param filename
 */
void Hair::writeStitches(const char* filename) {
    std::ofstream out(filename);
    writeStitches(out);
}

void Hair::writeStitches(const std::string filename) {
    writeStitches(filename.c_str());
}

void writeSVGHead(std::ostream& out) {
    out << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?><svg" << std::endl
        << "   xmlns:dc=\"http://purl.org/dc/elements/1.1/\""
        << "  xmlns:cc=\"http://creativecommons.org/ns#\""
        << "  xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\""
        << "  xmlns:svg=\"http://www.w3.org/2000/svg\""
        << "  xmlns=\"http://www.w3.org/2000/svg\""
        << "  xmlns:sodipodi=\"http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd\""
        << "  xmlns:inkscape=\"http://www.inkscape.org/namespaces/inkscape\""
        << "  width=\"400mm\""
        << "  height=\"350mm\""
        << "  viewBox=\"0 0 400 350\""
        << "  id=\"svg12765\""
        << "  version=\"1.1\">";
}

std::string Hair::outputStitch(Point p) {
    std::stringstream out;
    out << p.x() << "," << -p.y();
    return out.str();
}

void Hair::writeStitches(std::ostream& out) {
    if (countStitches() < 1) {
        return;
    }
    const Point offset = stitches.front().front();
    out <<  countStitches() + 1 + stitches.size() << std::endl;

    out << "0.0,0.0 color:0 flags:1" << std::endl;
    bool firstPatch = true;
    for (const auto& patch : stitches) {
        if (!firstPatch) {
            out << outputStitch(patch.front()-offset) << " color:0 flags:2" << std::endl;
        }
        firstPatch = false;
        out << outputStitch(patch.front()-offset) << " color:0 flags:1" << std::endl;
        for (const auto& stitch : patch) {
            out << outputStitch(stitch-offset) << " color:0 flags:0" << std::endl;
        }
    }
    out << outputStitch(stitches.back().back()-offset) << " color:0 flags:16";
}

void Hair::writeStitches(const std::vector<Point>& stitches, std::ostream& out) {
    if (stitches.empty()) {
        return;
    }
    const Point offset = getCenter(stitches);

    out << "0.0,0.0 color:0 flags:1" << std::endl
        << "0.0,0.0 color:0 flags:1" << std::endl
        << "0.0,0.0 color:0 flags:1" << std::endl
        << "0.0,0.0 color:0 flags:1" << std::endl;

    out << outputStitch(stitches.front() - offset) << " color:0 flags:1" << std::endl;
    out << outputStitch(stitches.front() - offset) << " color:0 flags:1" << std::endl;
    for (const Point& stitch : stitches) {
        out << outputStitch(stitch - offset) << " color:0 flags:0" << std::endl;
    }
    out << outputStitch(stitches.back() - offset) << " color:0 flags:16";
}

void Hair::writeStitches(const std::vector<Point>& stitches, const char* filename) {
    std::ofstream out(filename);
    writeStitches(stitches, out);
}

void Hair::writeStitches(const std::vector<Point>& stitches, const std::string filename) {
    writeStitches(stitches, filename.c_str());
}

double Hair::getStitchLengthes(const std::vector<Point>& stitches, std::vector<double>& lengths) {
    lengths.resize(stitches.size());
    if (stitches.size() < 2) {
        return 0;
    }
    lengths[0] = (stitches[0] - stitches[1]).length();
    double sum = lengths[0];
    for (size_t ii = 1; ii < stitches.size(); ++ii) {
        lengths[ii] = (stitches[ii] - stitches[ii-1]).length();
        sum += lengths[ii];
    }
    return sum;
}

Hair::Hair(Path _outline, Path _curve) : _outline(_outline), _curve(_curve) {}

void Hair::run() {
    clock_t start = clock();
    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;

    //makeAreaLarger(outline, areaGrow);

    clock_t start2 = clock();
    getCurves();
    clock_t stop2 = clock();
    std::cerr << "getCurves took " << (static_cast<double>(stop2-start2)) / CLOCKS_PER_SEC << std::endl << std::endl;
    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;

    start2 = clock();
    getStitches();
    stop2 = clock();
    std::cerr << "getStitches took " << (static_cast<double>(stop2-start2)) / CLOCKS_PER_SEC << std::endl << std::endl;
    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;

    start2 = clock();
    purgeOutside();
    stop2 = clock();
    std::cerr << "purgeOutside took " << (static_cast<double>(stop2-start2)) / CLOCKS_PER_SEC << std::endl << std::endl;
    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;

    start2 = clock();
    assignOutlineIntersections();
    stop2 = clock();
    std::cerr << "assignOutlineIntersections took " << (static_cast<double>(stop2-start2)) / CLOCKS_PER_SEC << std::endl << std::endl;
    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;

    start2 = clock();
    getOutlineIntermediateStitches();
    stop2 = clock();
    std::cerr << "getOutlineIntersections took " << (static_cast<double>(stop2-start2)) / CLOCKS_PER_SEC << std::endl << std::endl;
    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;

    start2 = clock();
    assembleAreas();
    stop2 = clock();
    std::cerr << "assemblePatches took " << (static_cast<double>(stop2-start2)) / CLOCKS_PER_SEC << std::endl << std::endl;
    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;

    start2 = clock();
    assembleGreedySolution();
    stop2 = clock();
    std::cerr << "assembleGreedySolution took " << (static_cast<double>(stop2-start2)) / CLOCKS_PER_SEC << std::endl << std::endl;
    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;

    clock_t stop = clock();
    std::cerr << "Calculation took " << (static_cast<double>(stop-start)) / CLOCKS_PER_SEC << std::endl << std::endl;
}

size_t getCorner(Path const& p) {
    size_t result = 0;
    double corner_angle = 0;
    for (size_t ii = 1; ii < p.size(); ++ii) {
        Point d1 = Geom::unit_vector(p[ii-1].pointAndDerivatives(1, 1)[1]);
        Point d2 = Geom::unit_vector(p[ii]  .pointAndDerivatives(0, 1)[1]);
        double current_angle = std::abs(Geom::angle_between(d1,d2));
        if (current_angle > corner_angle) {
            corner_angle = current_angle;
            result = ii;
        }
    }
    return result;
}

void Hair::runFeather() {

    _feather_corner = getCorner(_feather_template);

    std::cout << "Feather corner: " << _feather_corner << std::endl;

    getFeatherStitches();

    return;


    clock_t start = clock();
    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;

    //makeAreaLarger(outline, areaGrow);

    clock_t start2 = clock();
    getCurves();
    clock_t stop2 = clock();
    std::cerr << "getCurves took " << (static_cast<double>(stop2-start2)) / CLOCKS_PER_SEC << std::endl << std::endl;
    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;

    start2 = clock();
    getStitches();
    stop2 = clock();
    std::cerr << "getStitches took " << (static_cast<double>(stop2-start2)) / CLOCKS_PER_SEC << std::endl << std::endl;
    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;

    start2 = clock();
    purgeOutside();
    stop2 = clock();
    std::cerr << "purgeOutside took " << (static_cast<double>(stop2-start2)) / CLOCKS_PER_SEC << std::endl << std::endl;
    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;

    start2 = clock();
    assignOutlineIntersections();
    stop2 = clock();
    std::cerr << "assignOutlineIntersections took " << (static_cast<double>(stop2-start2)) / CLOCKS_PER_SEC << std::endl << std::endl;
    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;

    start2 = clock();
    getOutlineIntermediateStitches();
    stop2 = clock();
    std::cerr << "getOutlineIntersections took " << (static_cast<double>(stop2-start2)) / CLOCKS_PER_SEC << std::endl << std::endl;
    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;

    start2 = clock();
    assembleAreas();
    stop2 = clock();
    std::cerr << "assemblePatches took " << (static_cast<double>(stop2-start2)) / CLOCKS_PER_SEC << std::endl << std::endl;
    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;

    start2 = clock();
    assembleGreedySolution();
    stop2 = clock();
    std::cerr << "assembleGreedySolution took " << (static_cast<double>(stop2-start2)) / CLOCKS_PER_SEC << std::endl << std::endl;
    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;

    clock_t stop = clock();
    std::cerr << "Calculation took " << (static_cast<double>(stop-start)) / CLOCKS_PER_SEC << std::endl << std::endl;
}

void Hair::assembleGreedySolution() {

    greedySolution.clear();

    greedySolution.insert(greedySolution.end(), _areas[0].forward_stitches.begin(), _areas[0].forward_stitches.end());
    OutlineIntersection last_point = _areas[0].stop;

    for (size_t ii = 1; ii < _areas.size(); ++ii) {
        const EmbroideryArea& current_area = _areas[ii];
        std::vector<Point> forward_bridge = getShortestConnection(last_point, current_area.start);
        std::vector<Point> reverse_bridge = getShortestConnection(last_point, current_area.r_start);
        if (forward_bridge.size() < reverse_bridge.size()) {
            greedySolution.insert(greedySolution.end(), forward_bridge.begin(), forward_bridge.end());
            greedySolution.insert(greedySolution.end(),
                                  current_area.forward_stitches.begin(),
                                  current_area.forward_stitches.end());
            last_point = current_area.stop;
        } else {
            greedySolution.insert(greedySolution.end(), reverse_bridge.begin(), reverse_bridge.end());
            greedySolution.insert(greedySolution.end(),
                                  current_area.reverse_stitches.begin(),
                                  current_area.reverse_stitches.end());
            last_point = current_area.r_stop;
        }
    }

    addStartStop(greedySolution);

    std::vector<Point> greedy_solution_2;

    greedy_solution_2.insert(greedy_solution_2.end(), _areas[0].reverse_stitches.begin(), _areas[0].reverse_stitches.end());
    last_point = _areas[0].r_stop;

    for (size_t ii = 1; ii < _areas.size(); ++ii) {
        const EmbroideryArea& current_area = _areas[ii];
        std::vector<Point> forward_bridge = getShortestConnection(last_point, current_area.start);
        std::vector<Point> reverse_bridge = getShortestConnection(last_point, current_area.r_start);
        if (forward_bridge.size() < reverse_bridge.size()) {
            greedy_solution_2.insert(greedy_solution_2.end(), forward_bridge.begin(), forward_bridge.end());
            greedy_solution_2.insert(greedy_solution_2.end(),
                                     current_area.forward_stitches.begin(),
                                     current_area.forward_stitches.end());
            last_point = current_area.stop;
        } else {
            greedy_solution_2.insert(greedy_solution_2.end(), reverse_bridge.begin(), reverse_bridge.end());
            greedy_solution_2.insert(greedy_solution_2.end(),
                                     current_area.reverse_stitches.begin(),
                                     current_area.reverse_stitches.end());
            last_point = current_area.r_stop;
        }
    }

    addStartStop(greedy_solution_2);

    std::cout << "Greedy 1: " << greedySolution.size() << ", greedy 2: " << greedy_solution_2.size() << std::endl;

    if (greedy_solution_2.size() < greedySolution.size()) {
        greedySolution = greedy_solution_2;
    }
}

double Hair::curveLength(const Geom::Path& curve,
                         const Geom::PathTime start,
                         const Geom::PathTime stop,
                         const double dt) {
    if (stop < start) {
        return curveLength(curve, stop, start, dt);
    }
    double length = 0;
    Geom::Point lastPoint = curve.pointAt(start);
    Geom::PathTime currentTime = start;
    currentTime += dt;
    for (size_t ii = 0; ii < 100*1000*1000 && currentTime < stop; ++ii, currentTime += dt) {
        Geom::Point currentPoint = curve.pointAt(currentTime);
        length += (lastPoint - currentPoint).length();
        lastPoint = currentPoint;
    }
    return length + (lastPoint - curve.pointAt(stop)).length();
}

void Hair::getOutlineIntermediateStitches() {

    //* // Variant with simple outline
    // Conservative estimation for the step size we need:
    // d = 1000mm: maximum distance between two points of a embroidery design
    // h = 0.05mm: Smallest step we might want in the line length integration
    // dt = h/d = 1.0 / 20.000: time-step

    const double dt = 1.0 / (20 * 1000);

    const PathTime endTime = _intersections.back().time;
    PathTime lastTime = _intersections.back().time;
    double currentLength = 0;

    PathTime currentTime(lastTime);
    while(currentTime < PathTime(_outline.size(), 0.0)) {
        currentLength += (_outline.pointAt(lastTime) - _outline.pointAt(currentTime)).length();
        if (currentLength >= _stitch_length) {
            currentLength = 0;
            _outline_stitches.push_back(OutlineIntersection(
                                            _outline.pointAt(currentTime),
                                            currentTime,
                                            -1,
                                            -1,
                                            false));
        }
        lastTime = currentTime;
        currentTime += dt;
    }

    lastTime = PathTime();
    currentTime = PathTime();
    while(currentTime < PathTime(_outline.size(), 0.0)) {
        currentLength += (_outline.pointAt(lastTime) - _outline.pointAt(currentTime)).length();
        if (currentLength >= _stitch_length) {
            currentLength = 0;
            _outline_stitches.push_back(OutlineIntersection(
                                            _outline.pointAt(currentTime),
                                            currentTime,
                                            -1,
                                            -1,
                                            false));
        }
        lastTime = currentTime;
        currentTime += dt;
    }
    // */

    // Add all corners of the outline.
    for (size_t ii = 1; ii < _outline.size(); ++ii) {
        const Point a = -_outline[ii-1].pointAndDerivatives(.9999, 1)[1];
        const Point b = _outline[ii].pointAndDerivatives(.0001, 1)[1];
        const double angle = 1.0 + dot(a,b) / std::sqrt(dot(a,a) * dot(b,b));
        if (angle > 0.05) {
            currentTime = PathTime(ii, 0.0);
            _outline_stitches.push_back(OutlineIntersection(
                                            _outline.pointAt(currentTime),
                                            currentTime,
                                            -1,
                                            -1,
                                            false));
        }
    }
    {
        const Point a = -_outline.back().pointAndDerivatives(.9999, 1)[1];
        const Point b = _outline.front().pointAndDerivatives(.0001, 1)[1];
        const double angle = 1.0 + dot(a,b) / std::sqrt(dot(a,a) * dot(b,b));
        if (angle > 0.05) {
            currentTime = PathTime(0, 0.0);
            _outline_stitches.push_back(OutlineIntersection(
                                            _outline.pointAt(currentTime),
                                            currentTime,
                                            -1,
                                            -1,
                                            false));
        }
    }


    /* // Variant with prettier outline, not finished.
    for (size_t ii = 1; ii < intersections.size(); ++ii) {
        const OutlineIntersection& a = intersections[ii-1];
        const OutlineIntersection& b = intersections[ii];
        if (curveLength(outline, a.time, b.time) > stitchLength) {

        }
    }
    // */

    // Insert the outlineIntersections and make them optional stitches.
    for (const auto & p : _intersections) {
        _outline_stitches.push_back(OutlineIntersection(
                                        _outline.pointAt(p.time),
                                        p.time,
                                        true));
    }
    std::sort(_outline_stitches.begin(), _outline_stitches.end());

    // Now we need every point in the intersection vector to know
    // its corresponding point in the OutlineIntersectiones vector.
    size_t correspondingIndex = 0;
    for (size_t ii = 0; ii < _intersections.size(); ++ii) {
        OutlineIntersection& currentIntersection = _intersections[ii];
        bool failed = false;
        while (currentIntersection.time != _outline_stitches[correspondingIndex].time) {

            correspondingIndex++;
            if (correspondingIndex >= _outline_stitches.size()) {
                failed = true;
                break;
            }
        }
        if (failed) {
            std::cerr << "Failed to find corresponding point in OutlineIntersections in file "
                      << __FILE__ << ", line " << __LINE__ << std::endl;
            assert(false);
            break;
        }
        currentIntersection.outline_stitch_index = correspondingIndex;
    }
    for (EmbroideryLineLevel & level : _levels) {
        for (EmbroideryLine & line : level) {
            line.startInter.outline_stitch_index = _intersections[line.startInter.index].outline_stitch_index;
            line.endInter.outline_stitch_index   = _intersections[line.endInter.index].outline_stitch_index;
        }
    }
}

void Hair::writeArea(std::ofstream& out, const EmbroideryArea& area) {
    const std::string color = getColor(0,0);
    out << "<path style=\"display:inline;fill:none;fill-opacity:1;stroke:#"
        << color
        << ";stroke-width:" << _line_width << ";stroke-miterlimit:6;stroke-dasharray:none;stroke-opacity:1;enable-background:new"
        << "\" d=\"";
    for (const EmbroideryLine& line : area) {
        out << write_svg_path(getPath(line)) << " ";
    }
    out << "\"/>" << std::endl;
    out << "<g>";

    for (const EmbroideryLine& line: area) {
        for (const Point& p : line) {
            writeCircle(out, p, _stitch_point_radius, color);
        }
    }

    out << "</g>";
}

void Hair::writeAreas(const char* filename) {
    std::ofstream out(filename);
    writeSVGHead(out);

    for (const EmbroideryArea area: _areas) {
        writeArea(out, area);
    }

    const std::string color = getColor(0,0);
    out << "<path style=\"display:inline;fill:none;fill-opacity:1;stroke:#"
        << color
        << ";stroke-width:" << _line_width << ";stroke-miterlimit:6;stroke-dasharray:none;stroke-opacity:1;enable-background:new"
        << "\" d=\"";
    out << write_svg_path(getPath(_outline_stitches)) << " ";
    out << "\"/>" << std::endl;

    out << "<g>";

    for (const OutlineIntersection & s : _outline_stitches) {
        writeCircle(out, Point(s.x(), s.y()), _stitch_point_radius, color);
    }

    out << "</g>";

    out << "</svg>";
}

void Hair::writeForwardAreas(const char* filename) {
    std::ofstream out(filename);
    writeSVGHead(out);

    for (const EmbroideryArea area: _areas) {
        const std::string color = getColor(0,0);
        out << "<path style=\"display:inline;fill:none;fill-opacity:1;stroke:#"
            << color
            << ";stroke-width:" << _line_width << ";stroke-miterlimit:6;stroke-dasharray:none;stroke-opacity:1;enable-background:new"
            << "\" d=\"";
        out << write_svg_path(getPath(area.forward_stitches)) << " ";
        out << "\"/>" << std::endl;
        writeCircles(out, area.forward_stitches, _stitch_point_radius, color);

    }

    out << "</svg>";
}

void Hair::writeReverseAreas(const char* filename){
    std::ofstream out(filename);
    writeSVGHead(out);

    for (const EmbroideryArea area: _areas) {
        const std::string color = getColor(0,0);
        out << "<path style=\"display:inline;fill:none;fill-opacity:1;stroke:#"
            << color
            << ";stroke-width:" << _line_width << ";stroke-miterlimit:6;stroke-dasharray:none;stroke-opacity:1;enable-background:new"
            << "\" d=\"";
        out << write_svg_path(getPath(area.reverse_stitches)) << " ";
        out << "\"/>" << std::endl;
        writeCircles(out, area.reverse_stitches, _stitch_point_radius, color);
    }

    out << "</svg>";
}

bool Hair::isAdjacentLine(const EmbroideryLine& a, const EmbroideryLine& b) {
    if (   1 == moduloDist(a.startInter.index, b.startInter.index, _intersections.size())
           && 1 == moduloDist(a.endInter.index,   b.endInter.index,   _intersections.size())) {
        return true;
    }
    if (   1 == moduloDist(a.startInter.index, b.endInter.index,   _intersections.size())
           && 1 == moduloDist(a.endInter.index,   b.startInter.index, _intersections.size())) {
        return true;
    }
    return false;
}

void Hair::assembleAreas() {
    std::vector<std::vector<BOOL> > visited(_levels.size());

    bool foundFirstUnvisited = false;
    size_t firstUnvisitedLevel = 0;
    size_t totalLines = 0;
    for (size_t level = 0; level < _levels.size(); ++level) {
        std::cout << "We have " << _levels[level].size() << " lines on level " << level << ": ";
        for (size_t ii = 0; ii < _levels[level].size(); ++ii) {
            totalLines++;
            std::cout << _levels[level][ii].size() << ", ";
        }
        if (!foundFirstUnvisited && !_levels[level].empty()) {
            foundFirstUnvisited = true;
            firstUnvisitedLevel = level;
        }
        std::cout << std::endl;
        visited[level] = std::vector<BOOL>(_levels[level].size(), false);
    }
    EmbroideryLine& firstUnvisited = _levels[firstUnvisitedLevel][0];
    while (true) {
        bool foundUnvisited = false;

        for (size_t level = 0; level < _levels.size(); ++level) {
            for (size_t lineIndex = 0; lineIndex < _levels[level].size(); ++lineIndex) {
                if (!visited[level][lineIndex]) {
                    foundUnvisited = true;
                    firstUnvisited = _levels[level][lineIndex];
                    firstUnvisitedLevel = level;
                    visited[level][lineIndex] = true;
                    break;
                }
            }
            if (foundUnvisited) {
                break;
            }
        }
        if (!foundUnvisited) {
            break;
        }

        EmbroideryArea currentArea(firstUnvisited, firstUnvisitedLevel);
        EmbroideryLine& lastLine = firstUnvisited;
        for (size_t level = firstUnvisitedLevel + 1; level < _levels.size(); ++level) {
            bool foundAdjacent = false;
            for (size_t lineIndex = 0; lineIndex < _levels[level].size(); ++lineIndex) {
                if (!visited[level][lineIndex] && isAdjacentLine(lastLine, _levels[level][lineIndex])) {
                    lastLine = _levels[level][lineIndex];
                    currentArea.push_back(_levels[level][lineIndex]);
                    foundAdjacent = true;
                    visited[level][lineIndex] = true;
                    break;
                }
            }
            if (!foundAdjacent) {
                //std::cout << "no more adjacent lines at level " << level << std::endl;
                break;
            }
        }
        _areas.push_back(currentArea);
    }

    size_t visitedLines = 0;
    std::cout << "got " << _areas.size() << "areas with sizes:" << std::endl;
    for (const EmbroideryArea & area : _areas) {
        std::cout << area.size() << ", ";
        visitedLines += area.size();
    }
    std::cout << std::endl;
    std::cout << "min_levels: ";
    for (const EmbroideryArea & area : _areas) {
        std::cout << area.min_level << ", ";
    }
    std::cout << std::endl;
    std::cout << "visited " << visitedLines << " lines out of " << totalLines << std::endl;


    for (EmbroideryArea& area : _areas) {
        area.finish((*this));
    }

    std::sort(_areas.begin(), _areas.end());
    std::cout << "Areas:" << std::endl;
    std::cout << "min_level, forward_stitches.size(), size()" << std::endl;
    for (const EmbroideryArea& area : _areas) {
        std::cout << area.min_level << "\t" << area.forward_stitches.size() << "\t" << area.size() << std::endl;
    }
}

std::vector<Point> Hair::getShortestConnection(OutlineIntersection a, OutlineIntersection b) const {
    std::vector<Point> result_up;
    if (a == b) {
        return result_up;
    }
    std::vector<Point> result_down;
    if (b.outline_stitch_index < a.outline_stitch_index) {
        std::vector<Point> result = getShortestConnection(b, a);
        std::reverse(result.begin(), result.end());
        return result;
    }
    for (size_t index = a.outline_stitch_index; index < b.outline_stitch_index; ++index) {
        OutlineIntersection current = _outline_stitches[index];
        if (!current.optional) {
            result_up.push_back(Point(current.x(), current.y()));
        }
    }
    for (size_t index = b.outline_stitch_index; index < _outline_stitches.size(); ++index) {
        OutlineIntersection current = _outline_stitches[index];
        if (!current.optional) {
            result_down.push_back(Point(current.x(), current.y()));
        }
    }
    for (size_t index = 0; index < a.outline_stitch_index; ++index) {
        OutlineIntersection current = _outline_stitches[index];
        if (!current.optional) {
            result_down.push_back(Point(current.x(), current.y()));
        }
    }

    if (result_up.size() > result_down.size()) {
        std::reverse(result_down.begin(), result_down.end());
        return result_down;
    }
    return result_up;
}

size_t Hair::countStitches() {
    return countStitches(stitches);
}

template<class Inner>
size_t Hair::countStitches(const std::vector<Inner>& elements) {
    size_t result = 0;
    for (auto it : elements) {
        result += it.size();
    }
    return result;
}

std::vector<Point> Hair::purgeSmallStitches(const std::vector<Point>& orig, const double minLength) {
    size_t lastSize = orig.size();
    std::vector<Point> result = orig;
    while (true) {
        result = purgeSmallStitchesSub(result, minLength);
        if (lastSize == result.size()) {
            return result;
        }
        lastSize = result.size();
    }
    return result;
}

std::vector<Point> Hair::purgeSmallStitchesSub(const std::vector<Point>& orig, const double minLength) {
    std::vector<Point> result;
    result.reserve(orig.size());
    result.push_back(orig.front());
    size_t ii = 0;
    if (minLength > (orig[0] - orig[1]).length()) {
        result.push_back(orig[2]);
        ++ii;
    }
    for (; ii+3 < orig.size(); ++ii) {
        const Point& a = orig[ii];
        const Point& b = orig[ii+1];
        const Point& c = orig[ii+2];
        const Point& d = orig[ii+3];
        if ((b-c).length() > minLength) {
            result.push_back(b);
            continue;
        }
        if ((a-b).length() + (b-d).length() > (a-c).length() + (c-d).length()) {
            result.push_back(b);
            ii++;
        }
        else {
            result.push_back(c);
            ii++;
        }

    }
    if (minLength < (orig[orig.size() - 2] - orig[orig.size() - 1]).length()) {
        result.push_back(orig[orig.size() - 2]);
    }
    result.push_back(orig[orig.size() - 1]);
    return result;
}

std::vector<Point> Hair::subDivideLargeStitches(const std::vector<Point>& orig, const double maxLength) {
    if (orig.size() < 2) {
        return orig;
    }
    std::vector<Point> result;
    result.reserve(orig.size());
    result.push_back(orig.front());
    for (size_t ii = 0; ii+1 < orig.size(); ++ii) {
        const Point& a = orig[ii];
        const Point& b = orig[ii+1];
        if (maxLength < (a-b).length()) {

            const size_t intervals = static_cast<size_t>(std::ceil((a-b).length() / maxLength));
            for (size_t ii = 1; ii < intervals; ++ii) {
                result.push_back(a + ((b-a) * ii) / intervals);
            }
        }
        result.push_back(b);
    }
    return result;
}

/*
void Hair::getBoundaryDiscretization() {
    std::cerr << "########### getBoundaryDiscretization ###########" << std::endl;
    std::vector<Point> tmp;
    for (auto it = outline.begin(); it != outline.end(); ++it) {
        Path curve;
        curve.append(*it);
        std::vector<Point> tmpStitches;
        getStitches(tmpStitches, curve, discreteOutlineParam);
        if (tmp.empty()) {
            tmp = tmpStitches;
        }
        else if (tmpStitches.size() > 0) {
            if ((tmp.back() - tmpStitches.back()).length() < (tmp.back() - tmpStitches.front()).length()) {
                std::reverse(tmpStitches.begin(), tmpStitches.end());
                std::cerr << "Reversing tmp ordering" << std::endl;
            }
            tmp.insert(tmp.end(), tmpStitches.begin(), tmpStitches.end());
        }
        //std::cerr << "#it " << ii++ << std::endl;
    }
    std::cerr << "Outline size: " << tmp.size() << std::endl;
    lengthStats(tmp);
    //tmp = purgeSmallStitches(tmp, discreteOutlineParam*0.1);
    std::cerr << "Outline size after purging: " << tmp.size() << std::endl;
    lengthStats(tmp);
    //tmp = subDivideLargeStitches(tmp, discreteOutlineParam*1.5);
    std::cerr << "Outline size after subdividing: " << tmp.size() << std::endl;
    lengthStats(tmp);

    discreteOutline = tmp;
}
// */

void Hair::getBoundaryDiscretization() {
    std::cerr << "########### getBoundaryDiscretization ###########" << std::endl;
    std::vector<Point> tmp;
    for (auto it = _outline.begin(); it != _outline.end(); ++it) {
        Path curve;
        curve.append(*it);
        std::vector<Point> tmpStitches;
        getStitches(tmpStitches, curve, _discrete_outline_param);
        if (tmp.empty()) {
            tmp = tmpStitches;
        }
        else if (tmpStitches.size() > 0) {
            if ((tmp.back() - tmpStitches.back()).length() < (tmp.back() - tmpStitches.front()).length()) {
                std::reverse(tmpStitches.begin(), tmpStitches.end());
                std::cerr << "Reversing tmp ordering" << std::endl;
            }
            tmp.insert(tmp.end(), tmpStitches.begin(), tmpStitches.end());
        }
        //std::cerr << "#it " << ii++ << std::endl;
    }
    std::cerr << "Outline size: " << tmp.size() << std::endl;
    lengthStats(tmp);
    //tmp = purgeSmallStitches(tmp, discreteOutlineParam*0.1);
    std::cerr << "Outline size after purging: " << tmp.size() << std::endl;
    lengthStats(tmp);
    //tmp = subDivideLargeStitches(tmp, discreteOutlineParam*1.5);
    std::cerr << "Outline size after subdividing: " << tmp.size() << std::endl;
    lengthStats(tmp);

    _discrete_outline = tmp;
}

void Hair::lengthStats(const std::vector<Point>& points) {
    RunningStats stat;
    for (size_t ii = 0; ii+1 < points.size(); ++ii) {
        stat.push((points[ii] - points[ii+1]).length());
    }
    std::cerr << "Length stats: " << stat.print() << std::endl;
}

std::vector<std::pair<size_t, size_t> > Hair::endpointOutlineMatches(const std::vector<Point>& outline, const std::vector<std::vector<Point> >& lines) {
    std::vector<std::pair<size_t, size_t> > result(lines.size());
    for (size_t ii = 0; ii < lines.size(); ++ii) {

        result[ii] = std::pair<size_t, size_t>(bestMatch(outline, lines[ii].front()), bestMatch(outline, lines[ii].back()));
    }
    return result;
}

size_t Hair::bestMatch(const std::vector<Point>& points, const Point& target) {
    double bestDist = std::numeric_limits<double>::max();
    size_t match = 0;
    for (size_t jj = 0; jj < points.size(); ++jj) {
        const Point& o = points[jj];
        if (bestDist > (o-target).length()) {
            bestDist = (o-target).length();
            match = jj;
        }
    }
    return match;
}

size_t Hair::moduloDist(const size_t _a, const size_t _b, const size_t _total) {
    const int a = static_cast<int>(_a);
    const int b = static_cast<int>(_b);
    const int total = static_cast<int>(_total);
    size_t result = static_cast<size_t>(std::min(std::abs(a-b), std::min(std::abs(a- (b-total)), std::abs(b- (a-total)))));
    return result;
}

/**
 * @brief moduloDistDirection For given numbers _a and _b find the direction of the fastest
 * way from _a to _b modulo _total
 * @param _a
 * @param _b
 * @param _total
 * @return
 */
int Hair::moduloDistDirection(const size_t _a, const size_t _b, const size_t _total) {
    const int a = static_cast<int>(_a);
    const int b = static_cast<int>(_b);
    if (a < b) {
        if (static_cast<size_t>(b-a) <= moduloDist(_a,_b,_total)) {
            return 1;
        }
        return -1;
    }
    if (static_cast<size_t>(a-b) <= moduloDist(_a,_b,_total)) {
        return -1;
    }
    return 1;
}

void Hair::assembleStitches() {
    const size_t origCount = countStitches(stitches);
    std::vector<Point> bestStitches;
    std::vector<std::vector<Point> > bestBridges;
    size_t bestCount = std::numeric_limits<size_t>::max();

    typedef char BOOL;
    std::vector<BOOL> shouldCheck(stitches.size(), false);
    std::vector<BOOL> hasChecked(stitches.size(), false);
    shouldCheck.front() = true;
    shouldCheck.back() = true;

    while (true) {
        bool checkInLastIteration = false;

        for (size_t ii = 0; ii < stitches.size(); ++ii) {
            if (!shouldCheck[ii] || hasChecked[ii]) {
                continue;
            }
            hasChecked[ii] = true;
            checkInLastIteration = true;

            {
                std::vector<std::vector<Point> > currentBridges;
                const std::vector<Point> currentStitches = assembleStitches(stitches, currentBridges, ii, false, shouldCheck);
                const size_t currentCount = currentStitches.size();
                if (currentCount < bestCount) {
                    bestCount = currentCount;
                    bestStitches = currentStitches;
                    bestBridges = currentBridges;
                    std::cerr << "New best stitches in iteration " << ii << ", non-reversed" << std::endl;
                    std::cerr << "Stitches added by assembleStitches: " << bestCount - origCount
                              << " (" << 100*static_cast<double>(bestCount - origCount)/origCount << "%)" << std::endl;
                    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;
                }
            }
            {
                std::vector<std::vector<Point> > currentBridges;
                const std::vector<Point>currentStitches = assembleStitches(stitches, currentBridges, ii, true, shouldCheck);
                const size_t currentCount = currentStitches.size();
                if (currentCount < bestCount) {
                    bestCount = currentCount;
                    bestStitches = currentStitches;
                    bestBridges = currentBridges;
                    std::cerr << "New best stitches in iteration " << ii << ", non-reversed" << std::endl;
                    std::cerr << "Stitches added by assembleStitches: " << bestCount - origCount
                              << " (" << 100*static_cast<double>(bestCount - origCount)/origCount << "%)" << std::endl;
                    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;
                }
            }
        }
        if (!checkInLastIteration) {
            break;
        }
    }
    RunningStats checkedLines;
    for (auto val : hasChecked) {
        checkedLines.push(val ? 1 : 0);
    }
    std::cout << "Checked lines: " << checkedLines.print() << std::endl;

    stitches.clear();
    stitches.push_back(bestStitches);
    _bridges = bestBridges;
}


std::vector<Point> Hair::assembleStitches(const std::vector<std::vector<Point> >& stitches, std::vector<std::vector<Point> >& bridges, const size_t initialLine, const bool reverseInitial, std::vector<BOOL>& shouldCheck) {
    // vector<bool> is slow in this case, the memory saving doesn't seem
    // to outweigh the bit manipulation.
    typedef char BOOL;
    std::vector<BOOL> visited(stitches.size(), false);
    visited[initialLine] = true;
    std::vector<Point> patch = stitches[initialLine];
    patch.reserve((stitches.size() * 4) / 3);
    bridges.clear();
    if (reverseInitial) {
        std::reverse(patch.begin(), patch.end());
    }
    const size_t discreteSize = _discrete_outline.size();

    std::vector<std::vector<Point> > patches;
    std::vector<Point> smallPatch;

    int lastLine = -1;
    std::vector<Point> bridge;
    bridge.reserve(128);
    while (true) {
        const Point& lastPoint = patch.back();
        int currentLine = getNearestLine(lastPoint, visited, stitches);
        //int currentLine = getNextLine(lastLine, lastPoint, visited, stitches);
        lastLine = currentLine;
        if (currentLine < 0) {
            break;
        }
        size_t line = static_cast<size_t>(currentLine);
        visited[line] = true;
        const std::vector<Point> &linePoints = stitches[line];
        bridge.clear();
        bridge.reserve(128);
        bool reverseCurrentLine = false;
        //if (maxStitchLength < (lastPoint - linePoints.front()).length()) {
        {
            // Now we need to find a path along the outline from the lastPoint to the current line.
            // First we find out if the start or the end of the line are within smaller range.

            const size_t lastPointOutline = bestMatch(_discrete_outline, lastPoint);
            const size_t frontOption = bestMatch(_discrete_outline, linePoints.front());
            const size_t backOption = bestMatch(_discrete_outline, linePoints.back());
            size_t option = frontOption;
            if (moduloDist(lastPointOutline, frontOption, discreteSize)
                    > moduloDist(lastPointOutline, backOption, discreteSize)) {
                option = backOption;
                reverseCurrentLine = true;
            }
            // We add the stitches along the way to the chosen option.
            int direction = moduloDistDirection(lastPointOutline, option, discreteSize);

            const std::size_t maxBridgeSize = discreteSize;
            int ii = static_cast<int>(lastPointOutline);
            while (ii != static_cast<int>(option)) {
                ii += direction;
                if (ii < 0) {
                    ii = static_cast<int>(discreteSize) + ii;
                }
                if (ii >= static_cast<int>(discreteSize)) {
                    ii = 0;
                }
                bridge.push_back(_discrete_outline[static_cast<size_t>(ii)]);
                if (bridge.size() > maxBridgeSize) {
                    std::cerr << "Bridge has grown too large, maximum size is " << maxBridgeSize << ", lasPointOutline is " << lastPointOutline << ", option is " << option << ", direction is " << direction << std::endl;
                    break;
                }
            }
            bridge.push_back(_discrete_outline[option]);
            // We remove unnecessary stitches from the bridge.
            const Point& nextLineStart = reverseCurrentLine ? linePoints.back() : linePoints.front();
            if (bridge.size() > 1) {
                const Point& lastBridge = bridge.back();
                const Point& lastlastBridge = bridge[bridge.size()-2];
                if (dot(lastlastBridge - lastBridge, nextLineStart - lastBridge) > 0
                        && (lastlastBridge - nextLineStart).length() < _stitch_length/2) {
                    bridge.pop_back();
                }
            }
            if (bridge.size() > 1) {
                const Point& firstBridge = bridge.front();
                const Point& secondBridge = bridge[1];
                if (dot(secondBridge - firstBridge, lastPoint - firstBridge) > 0
                        && (secondBridge - lastPoint).length() < _stitch_length/2) {
                    bridge.erase(bridge.begin());
                }
            }
            if (bridge.size() == 1) {
                const Point& bridgeElement = bridge.front();
                if (dot(lastPoint - bridgeElement, nextLineStart - bridgeElement) > 0
                        || (lastPoint - nextLineStart).length() < _stitch_length/3) {
                    bridge.clear();
                }
            }
            if (!bridge.empty()) {
                bridges.push_back(bridge);
                patch.insert(patch.end(), bridge.begin(), bridge.end());
                patches.push_back(smallPatch);
                smallPatch.clear();
                //shouldCheck[line] = true;
                if (static_cast<size_t>(currentLine) >= shouldCheck.size()) {
                    std::cerr << "Error: shouldCheck.size: " << shouldCheck.size() << ", currentLine: " << currentLine << std::endl;
                }
                else {
                    shouldCheck[static_cast<size_t>(currentLine)] = true;
                }
            }
        }
        //else {
        //std::cerr << "added bridge stitches: " << bridge.size() << std::endl;
        //patches.push_back(patch);
        //patch.clear();
        //}
        if (reverseCurrentLine) {
            patch.insert(patch.end(), linePoints.rbegin(), linePoints.rend());
            smallPatch.insert(smallPatch.end(), linePoints.rbegin(), linePoints.rend());
        }
        else {
            patch.insert(patch.end(), linePoints.begin(), linePoints.end());
            smallPatch.insert(smallPatch.end(), linePoints.begin(), linePoints.end());
        }
    }
    return patch;
}

int Hair::getNearestLine(const Point& p, const std::vector<BOOL>& visited, const std::vector<std::vector<Point> >& lines) {
    int result = -1;
    double bestDistance = std::numeric_limits<double>::max();
    for (size_t ii = 0; ii < visited.size(); ++ii) {
        if (visited[ii]) {
            continue;
        }
        if (bestDistance > (p - lines[ii].front()).length()) {
            bestDistance = (p - lines[ii].front()).length();
            result = static_cast<int>(ii);
        }
        if (bestDistance > (p - lines[ii].back()).length()) {
            bestDistance = (p - lines[ii].back()).length();
            result = static_cast<int>(ii);
        }
    }
    return result;
}

int Hair::getNextLine(const int lastLine, const Point& p, const std::vector<BOOL>& visited, const std::vector<std::vector<Point> >& lines) {
    if (lastLine < 0) {
        return 0;
    }
    if (lastLine + 1 < lines.size()) {
        if (_stitch_length >= (p - lines[lastLine+1].front()).length()
                || _stitch_length >=(p - lines[lastLine+1].back()).length()) {
            if (!visited[lastLine+1]) {
                return lastLine+1;
            }
        }
    }
    int result = -1;
    double bestDistance = std::numeric_limits<double>::max();
    for (size_t ii = 0; ii < visited.size(); ++ii) {
        if (visited[ii]) {
            continue;
        }
        return static_cast<int>(ii);
        if (bestDistance > (p - lines[ii].front()).length()) {
            bestDistance = (p - lines[ii].front()).length();
            result = static_cast<int>(ii);
        }
        if (bestDistance > (p - lines[ii].back()).length()) {
            bestDistance = (p - lines[ii].back()).length();
            result = static_cast<int>(ii);
        }
    }
    return result;
}

/**
 * @brief purgeOutside Removes all stitches outside the given boundary
 * and replaces stitches crossing the boundary by stitches ending / starting at the boundary.
 */
void Hair::purgeOutside() {
    _stitch_length_stats.clear();
    std::vector<std::vector<Point> > newStitches;
    _levels = std::vector<EmbroideryLineLevel>(stitches.size());
    //#pragma omp parallel for
    for (size_t level = 0; level < stitches.size(); ++level) {
        const std::vector<Point>& line = stitches[level];
        _levels[level].reserve(line.size());
        if (line.size() < 2) {
            continue;
        }
        std::vector<Point> newLine;
        std::vector<BOOL> inside(line.size());
        for (size_t ii = 0; ii < line.size(); ++ii) {
            inside[ii] = isInsidePoint(line[ii]);
        }
        if (inside[0]) {
            std::cerr << "First point must not be inside the area in file "
                      << __FILE__ << ", line " << __LINE__ << std::endl
                      << ", level " << level << ", point " << line[0] << std::endl;
            return;
            newLine.push_back(line[0]);
        }
        OutlineIntersection startInter, endInter;
        startInter.height = level;
        endInter.height = level;
        for (size_t ii = 0; ii + 1 < line.size(); ++ii) {
            if (inside[ii] && inside[ii+1]) {
                newLine.push_back(line[ii+1]);
                _stitch_length_stats.push((line[ii] - line[ii+1]).length());
                continue;
            }
            BezierCurve _straightLine = BezierCurveN<1>(line[ii], line[ii+1]);
            Path straightLine;
            straightLine.append(_straightLine);
            std::vector<PathIntersection> intersection = straightLine.intersect(_outline);
            if (inside[ii] && !inside[ii+1]) {
                //straightLine.stitchTo(line[ii+1]);
                //std::cerr << "straightLine: " << write_svg_path(straightLine) << std::endl;
                Point intersectPoint = straightLine.pointAt(intersection[0].first);
                endInter.setPoint(intersectPoint);
                endInter.time = intersection[0].second;
                newLine.push_back(intersectPoint);
                if (lineLength(newLine) > _min_stitch_length) {
#pragma omp critical
                    {
                        newStitches.push_back(newLine);
                        _levels[level].push_back(EmbroideryLine(newLine, level, startInter, endInter));
                        startInter.line = &_levels[level].back();
                        endInter.line = &_levels[level].back();
                        _intersections.push_back(startInter);
                        _intersections.push_back(endInter);
                    }
                }
                newLine.clear();
                continue;
            }
            if (!inside[ii] && inside[ii+1]) {
                Path straightLine;
                straightLine.append(_straightLine);
                //straightLine.stitchTo(line[ii+1]);
                //std::cerr << "straightLine: " << write_svg_path(straightLine) << std::endl;
                if (intersection.empty()) {
                    std::cout << "Intersection was unexpectedly empty in file "
                              << __FILE__ << ", line " << __LINE__ << std::endl;
                }
                Point intersectPoint = straightLine.pointAt(intersection[0].first);
                newLine.push_back(intersectPoint);
                newLine.push_back(line[ii+1]);
                startInter.time = intersection[0].second;
                startInter.setPoint(intersectPoint);
                continue;
            }
            if (!inside[ii] && !inside[ii+1]) {
                if (2 <= intersection.size()) {
                    const Point start_point = straightLine.pointAt(intersection[0].first);
                    const Point end_point = straightLine.pointAt(intersection[1].first);
                    newLine.push_back(start_point);
                    newLine.push_back(end_point);
                    startInter.time = intersection[0].second;
                    startInter.setPoint(start_point);
                    endInter.time = intersection[1].second;
                    endInter.setPoint(end_point);
                    if (lineLength(newLine) > _min_stitch_length) {
#pragma omp critical
                        {
                            newStitches.push_back(newLine);
                            _levels[level].push_back(EmbroideryLine(newLine, level, startInter, endInter));
                            startInter.line = &_levels[level].back();
                            endInter.line = &_levels[level].back();
                            _intersections.push_back(startInter);
                            _intersections.push_back(endInter);
                        }
                    }
                    newLine.clear();
                }
                continue;
            }
        }
    }
    stitches = newStitches;
    // std::cout << "Number of intersections: " << intersections.size() << std::endl;
    std::cout << "Stitch stats: " << _stitch_length_stats.print() << std::endl;

}

double Hair::lineLength(const std::vector<Point>& points) {
    double result = 0;
    for (size_t ii = 1; ii < points.size(); ++ii) {
        result += (points[ii-1] - points[ii]).length();
    }
    return result;
}

void Hair::printStats() {
    std::cout << "Stitch length stats: " << _stitch_length_stats.print() << std::endl;
}

Path moveAlongPath(
        Path const& curve,
        Path const& template_path,
        PathTime target
        ) {
    Path result = template_path;
    if (curve.size() < 1) {
        return result;
    }
    if (target.curve_index >= curve.size()) {
        target.curve_index = curve.size() - 1;
        target.t = 1;
    }

    Geom::Affine movement, tmp_mov;
    movement.setIdentity();
    movement.setTranslation(-curve.initialPoint());

    Point const targetAngle = Geom::unit_vector(curve[target.curve_index].pointAndDerivatives(target.t, 1)[1]);
    Point const sourceAngle = Geom::unit_vector(curve[0].pointAndDerivatives(0, 1)[1]);
    tmp_mov.setIdentity();
    tmp_mov.setRotation(Geom::angle_between(targetAngle, sourceAngle));
    movement *= tmp_mov;

    tmp_mov.setIdentity();
    tmp_mov.setTranslation(curve.pointAt(target));

    movement *= tmp_mov;

    result *= movement;

    return result;
}

double meanDist(Path const& template_path, Path const& offsetted) {
    double const stepsize = 0.1;

    double last_distance = 0;
    double sum = 0;
    double weight_sum = 0;
    Point current_point = offsetted.initialPoint();
    for (size_t ii = 0; ii < offsetted.size(); ++ii) {
        for (double t = 0; t <= 1; t += stepsize) {
            Point const next_point = offsetted[ii].pointAt(t);

            double const next_distance = Geom::distance(next_point, current_point);
            double const weight = next_distance + last_distance;
            Point const template_point = template_path.pointAt(template_path.nearestTime(current_point));

            double const dist = Geom::distance(current_point, template_point);

            sum += dist * weight;
            weight_sum += weight;

            current_point = next_point;
            last_distance = next_distance;
        }
    }
    Point const template_point = template_path.pointAt(template_path.nearestTime(current_point));

    double const dist = Geom::distance(current_point, template_point);
    double const weight = last_distance;
    sum += dist * weight;
    weight_sum += weight;

    return sum / weight_sum;
}

bool offsetTemplateOnCurve(
        Path const& template_path,
        Path const& last_path,
        Path const& curve,
        PathTime const last_time,
        PathTime & new_time,
        double const target_distance,
        double &offset) {

    if (offset <= 0) {
        offset = .1;
    }
    double min_offset = 0;
    double max_offset = std::numeric_limits<double>::max();

    Path offsetted;

    for (size_t ii = 0; ii < 100; ++ii) {
        new_time = plus(last_time, offset);
        if (new_time.curve_index >= curve.size()) {
            max_offset = offset;
            std::cout << "overflow max_offset: " << max_offset << std::endl;
            break;
        }
        offsetted = moveAlongPath(curve, template_path, new_time);
        double const current_dist = meanDist(last_path, offsetted);
        if (current_dist > target_distance) {
            max_offset = offset;
            std::cout << "max_offset: " << max_offset << std::endl;
            break;
        }
        offset *= 2;
    }
    if (!std::isfinite(max_offset)) {
        return false;
    }

    for (size_t ii = 0; ii < 100; ++ii) {
        offset = (max_offset + min_offset) / 2.;
        new_time = plus(last_time, offset);
        offsetted = moveAlongPath(curve, template_path, new_time);
        double const current_dist = meanDist(last_path, offsetted);
        if (current_dist > target_distance) {
            max_offset = offset;
        }
        else {
            min_offset = offset;
        }
        if (std::abs(max_offset - min_offset) < 1e-6) {
            break;
        }
    }

    offset = (max_offset + min_offset) / 2.;
    std::cout << "final offset: " << offset << std::endl;

    new_time = plus(last_time, offset);

    return true;
}

void Hair::getFeatherStitches() {
    std::ofstream debug("feather-debug.svg");
    writeSVGHead(debug);
    write (debug, _outline, "ff0000");
    write (debug, _feather_curve, "00ff00");
    write (debug, _feather_template, "0000ff");

    _feather_corner = getCorner(_feather_template);
    _feather_corner_point = _feather_template[_feather_corner].initialPoint();
    PathTime closest = _feather_curve.nearestTime(_feather_corner_point);

    Path feather_template = _feather_template;
    Geom::Affine movement;
    movement.setIdentity();
    movement.setTranslation(- _feather_corner_point);

    feather_template *= movement;

    Point sourceAngle = Geom::unit_vector(_feather_curve[closest.curve_index].pointAndDerivatives(closest.t, 1)[1]);
    Point targetAngle = Geom::unit_vector(_feather_curve[0].pointAndDerivatives(0, 1)[1]);
    movement.setIdentity();
    double const angle = Geom::angle_between(targetAngle, sourceAngle);
    movement.setRotation(angle);

    feather_template *= movement;

    movement.setIdentity();
    movement.setTranslation(_feather_curve.initialPoint());
    feather_template *= movement;

    write(debug, feather_template, "000000");

    PathTime last_time(0,0);
    Path last_path = feather_template;
    double offset = 0;

    for (size_t ii = 0; ii < 500; ++ii) {
        PathTime new_time = last_time;
        bool const success = offsetTemplateOnCurve(
                    feather_template,
                    last_path,
                    _feather_curve,
                    last_time,
                    new_time,
                    _line_spacing,
                    offset);
        if (!success) {
            break;
        }
        last_path = moveAlongPath(
                    _feather_curve,
                    feather_template,
                    new_time
                    );
        last_time = new_time;

        if (!_outline.intersect(last_path).empty()) {
            write(debug, last_path, "000000");
        }

    }


    debug << "</svg>" << std::endl;
}

void Hair::getStitches() {

    const size_t initialIndex = _curves.size()/2;
    Path initialCurve = _curves[initialIndex];

    stitches.assign(_curves.size(), std::vector<Point>());
    getStitches(stitches[initialIndex], initialCurve);

    //enlargeCurves();

    _initial_stitch_lengths_sum = getStitchLengthes(stitches[initialIndex], _initial_stitch_lengths);

    for (size_t ii = initialIndex+1; ii < _curves.size(); ++ii) {
        auto tmp = projection(stitches[ii-1], _curves[ii]);
        stitches[ii] = tmp;
    }
    for (size_t ii = initialIndex; ii+1 > 0 && ii + 1 < stitches.size(); --ii) {
        auto tmp = projection(stitches[ii+1], _curves[ii]);
        stitches[ii] = tmp;
    }
    //        for (size_t ii = 0; ii < curves.size(); ii += 2) {
    //            stitches[ii] = vectorOffset(stitches[ii], offset);
    //        }

    double currentOffset = 0;
    for (size_t ii = 0; ii < _curves.size(); ++ii) {
        stitches[ii] = vectorOffsetOnCurve(stitches[ii], _curves[ii], currentOffset);
        currentOffset = fmod(currentOffset + _offset, 1.0);
    }


}

void Hair::assignOutlineIntersections() {
    std::sort(_intersections.begin(), _intersections.end());
    bool failed = false;
    for (EmbroideryLineLevel & level : _levels) {
        for (EmbroideryLine & line : level) {
            auto it = std::find(_intersections.begin(), _intersections.end(), line.startInter);
            if (_intersections.end() == it) {
                std::cerr << "Couldn't find the intersection in the intersections vector in file "
                          << __FILE__ << ", line " << __LINE__ << std::endl
                          << line << std::endl;
                failed = true;
            }
            else {
                it->index = it - _intersections.begin();
                line.startInter.index = it->index;
                it->line = &line;
            }

            it = std::find(_intersections.begin(), _intersections.end(), line.endInter);
            if (_intersections.end() == it) {
                std::cerr << "Couldn't find the intersection in the intersections vector in file "
                          << __FILE__ << ", line " << __LINE__ << std::endl
                          << line << std::endl;
                failed = true;
            }
            else {
                it->index = it - _intersections.begin();
                line.endInter.index = it->index;
                it->line = &line;
            }
        }
    }
    for (const EmbroideryLineLevel & level : _levels) {
        for (const EmbroideryLine & line : level) {
            if (line.startInter != _intersections[line.startInter.index]) {
                std::cerr << "Outline intersection assignment incorrect in file "
                          << __FILE__ << ", line " << __LINE__ << std::endl;
                failed = true;
            }
        }
    }
    for (const OutlineIntersection & inter : _intersections) {
        if (inter != inter.line->endInter && inter != inter.line->startInter) {
            failed = true;
            std::cout << "Outline intersection assignment incorrent in file "
                      << __FILE__ << ", line " << __LINE__ << std::endl;
            std::cout << "inter: " << inter << std::endl
                      << "startInter: " << inter.line->startInter << std::endl
                      << "endInter: " << inter.line->endInter << std::endl;
        }
    }
    if (failed) {
        std::cout << "Failed assigning outline intersections" << std::endl;
    }
    else {
        std::cout << "Succeeded in assigning outline intersection" << std::endl;
    }
}

void Hair::enlargeCurves() {
    for (auto &curve : _curves) {
        enlargeCurve(curve);
    }
}
void Hair::enlargeCurve(Path& curve) {
    enlargeCurveEnd(curve);
    curve = curve.reversed();
    enlargeCurveEnd(curve);
    curve = curve.reversed();
}

void Hair::enlargeCurveEnd(Path& curve) {
    PathTime t = curve.nearestTime(curve.finalPoint());
    std::vector<Point> data = curve.at(t.curve_index).pointAndDerivatives(t.t, 1);
    const double length = (curve.initialPoint() - curve.finalPoint()).length();
    Point a = data[0];
    Point b = a + data[1] * length / data[1].length();
    BezierCurve line = BezierCurveN<1>(a, b);
    curve.append(line);
}

/**
 * @brief atLeastOneInside
 * @param points
 * @param boundary
 * @return
 */
bool Hair::atLeastOneInside(const std::vector<Point>& points, const Path& boundary) {
    for (auto point : points) {
        if (0 != boundary.winding(point)) {
            return true;
        }
    }
    return false;
}

/**
 * @brief vectorOffset moves every point in a vector in the direction given by the difference between itself and the previous point.
 * @param points The original points
 * @param offset The relative offset, usually between -1 and 1
 * @return The moved points.
 */
std::vector<Point> Hair::vectorOffsetOnCurve(
        std::vector<Point> const& points,
        Path const& curve,
        double const offset) {
    std::vector<Point> result = points;
    if (points.size() < 2) {
        return result;
    }
    Point last_point = points[1];
    for (size_t ii = 0; ii < points.size(); ++ii) {
        PathTime t = curve.nearestTime(points[ii]);

        //double const error = Geom::distance(points[ii], curve.pointAt(t));
        //std::cout << error << std::endl; // This is ok, something between 0 and 1e-12

        const double distance = Geom::distance(last_point, points[ii]);
        result[ii] = getOffsettedPointOnCurve(curve, t, offset * distance);
        last_point = points[ii];
    }

    return result;
}

bool Hair::isInsidePoint(const Point p) {
    return 0 != _outline.winding(p);
}

std::vector<Point> Hair::projection(const std::vector<Point> & src, Path& curve) {
    std::vector<Point> result;
    result.reserve(src.size());
    if (src.size() < 2) {
        return result;
    }

    Point lastPoint = src[1];
    for (auto p : src) {
        PathTime t = curve.nearestTime(p);
        result.push_back(curve.pointAt(t));
        lastPoint = p;
    }
    return result;
}

template<class Points>
void Hair::getStitches(Points& initialCurvePoints, const Path& initialCurve) {
    getStitches(initialCurvePoints, initialCurve, _stitch_length);
}

template<class Points>
void Hair::getStitches(Points& curvePoints, const Path& curve, const double stitchLength) {
    PathTime t;
    Point currentPoint = curve.pointAt(t);
    curvePoints.push_back(currentPoint);
    Point lastPoint = currentPoint;
    RunningStats stat;
    while (true) {
        currentPoint = getOffsettedPointOnCurve(curve, t, stitchLength);
        curvePoints.push_back(currentPoint);
        stat.push((currentPoint - lastPoint).length());
        //std::cerr << "Size: " << curvePoints.size() << std::endl;
        //std::cerr << currentPoint << std::endl;

        // Test if we have reached the end of the curve
        if (0.01 > (currentPoint - curve.finalPoint()).length()
                || stitchLength/2 > (currentPoint - lastPoint).length()
                || stitchLength*2 < (currentPoint - lastPoint).length()) {
            break;
        }
        bool finished = false;
        for (size_t ii = 0; ii+1 < curvePoints.size(); ++ii) {
            if (stitchLength/2 > (currentPoint - curvePoints[ii]).length()) {
                finished = true;
                break;
            }
        }
        if (finished) {
            break;
        }
        lastPoint = currentPoint;
    }
    //std::cerr << "getStitches lengths: " << stat.print() << std::endl;
}


Point Hair::getOffsettedPointOnCurve(const Path& curve, PathTime& t, const double targetOffset) {
    /*
    std::vector<Point> const data = curve.at(t.curve_index).pointAndDerivatives(t.t, 1);

    Point initialGuess = data[0] + data[1] * (targetOffset / data[1].length());

    PathTime projectionTime = curve.nearestTime(initialGuess);

    Point projection = curve.pointAt(projectionTime);

    t = projectionTime;
    //std::cerr << "Relative error after 1 projection: " << std::abs(1.0 - (data[0] - projection).length() / stitchLength) << std::endl;

    double offset = targetOffset - (data[0] - projection).length();

    for (size_t ii = 0; ii < 100 && offset > 1e-16; ++ii) {

        offset = targetOffset - (data[0] - projection).length();

        std::vector<Point> data2 = curve.at(projectionTime.curve_index).pointAndDerivatives(projectionTime.t, 1);

        Point secondGuess = data2[0] + data2[1] * offset / data2[1].length();

        projectionTime = curve.nearestTime(secondGuess);

        projection = curve.pointAt(projectionTime);
        //std::cerr << "Relative error after " << 2 + ii << " projections: " << std::abs(1.0 - (data[0] - projection).length() / targetOffset) << std::endl;
    }
    t = projectionTime;
    return projection;

    // */

    //*
    Geom::Point const orig_point = curve.pointAt(t);
    double dt = .1;
    double max_dt = 0;
    Geom::Point new_point;
    double new_distance = 0;
    for (size_t ii = 0; ii < 100; ++ii) {
        new_point = curve.pointAt(plus(t, dt));
        new_distance = Geom::distance(new_point, orig_point);
        if (new_distance > targetOffset) {
            max_dt = dt;
            break;
        }
        dt *= 2;
    }
    double min_dt = 0;
    for (size_t ii = 0; ii < 100; ++ii) {
        new_point = curve.pointAt(plus(t, dt));
        new_distance = Geom::distance(new_point, orig_point);
        if (new_distance > targetOffset) {
            max_dt = dt;
        }
        else {
            min_dt = dt;
        }
        dt = .5 * (max_dt + min_dt);
        if (std::abs(new_distance - targetOffset) < 1e-6) {
            break;
        }
    }
    new_point = curve.pointAt(plus(t, dt));
    new_distance = Geom::distance(new_point, orig_point);
    if (std::abs(new_distance - targetOffset) > 1e-2) {
        std::cout << "Expected " << targetOffset << ", got " << new_distance << std::endl;
    }
    t = curve.nearestTime(new_point);

    return curve.pointAt(t);
    // */
}

void Hair::makeAreaLarger(Path& curve, const double offset) {
    curve = Inkscape::half_outline(curve, offset, _miter);
}

void Hair::getCurves() {

    PathVector rights, lefts;

    Path left = _curve;
    Path right = _curve.reversed();
    //#pragma omp parallel sections
    // These sections seem independent but somehow they are not
    // and the program crashes with a segfault some times when run in parallel
    {
        //#pragma omp section
        {
            for (size_t ii = 0; ii < _max_iter; ++ii) {
                try {
                    left = Inkscape::half_outline(left, _line_spacing, _miter, _join_type);
                }
                catch(std::runtime_error const& e) {
                    std::cout << "Cauch error in left half-outline part: " << e.what() << std::endl;
                    break;
                }
                if (left.intersect(_outline).empty()) {
                    break;
                }
                lefts.push_back(left);
                std::cout << "l " << ii << " " << left.size() << std::endl;
            }
        }
        //#pragma omp section
        {
            for (size_t ii = 0; ii < _max_iter; ++ii) {
                try {
                    right = Inkscape::half_outline(right, _line_spacing, _miter, _join_type);
                }
                catch(std::runtime_error const& e) {
                    std::cout << "Cauch error in right half-outline part: " << e.what() << std::endl;
                    break;
                }
                if (right.intersect(_outline).empty()) {
                    break;
                }
                rights.push_back(right);
                std::cout << "r " << ii << " " << right.size() << std::endl;
            }
            rights.reverse(true);
            _curves.insert(_curves.begin(), rights.begin(), rights.end());
            _curves.push_back(_curve);
        }
    }

    _curves.insert(_curves.end(), lefts.begin(), lefts.end());
    std::cerr << "We now have " << _curves.size() << " curves" << std::endl;
}

void Hair::write(const char* filename) {
    std::ofstream out(filename);
    write(out);
}

void Hair::write(const std::string filename) {
    write(filename.c_str());
}




void Hair::write(std::ostream& out) {
    writeSVGHead(out);
    out << "<g>";
    for (auto c : _curves) {
        write(out, c, "aaaaaa");
    }
    out << "</g>";
    write(out, _curve, "ff0000");
    write(out, _outline, "0000ff");

    //writeStitchPoints(out, "00ff00");
    out << "<g>";
    std::string color = getColor(0, greedySolution.size());
    write(out, getPath(greedySolution), color);
    writeCircles(out, greedySolution, _stitch_point_radius, color);
    out << "</g>";

    out << "<g>";
    write(out, getPath(_discrete_outline), "ff0000");
    writeCircles(out, _discrete_outline, _stitch_point_radius, "ff0000");
    out << "</g>";

    out << "<g>";
    writePatches(out, _bridges);
    out << "</g>";
    out << "<g>";
    writeOutlineIntersections(out, _intersections);
    out << "</g>";
    out << "</svg>" << std::endl;
}


void Hair::writeCurves(const char* filename) {
    std::ofstream out(filename);
    writeCurves(out);
}

void Hair::writeCurves(const std::string filename) {
    writeCurves(filename.c_str());
}




void Hair::writeCurves(std::ostream& out) {
    writeSVGHead(out);
    out << "<g>";
    for (auto c : _curves) {
        write(out, c, "aaaaaa");
    }
    out << "</g>";
    write(out, _curve, "ff0000");
    out << "</svg>" << std::endl;
}

void Hair::writeOutlineIntersections(std::ostream& out, std::vector<OutlineIntersection> &intersections) {
    std::vector<Point> points;
    for (auto inter: intersections) {
        points.push_back(_outline.pointAt(inter.time));
    }
    const auto path = getPath(points);
    write(out, getPath(points), "ff0000");

}

void Hair::write(std::ostream& out, Path path, std::string color) {
    out << "<path style=\"display:inline;fill:none;fill-opacity:1;stroke:#"
        << color
        << ";stroke-width:" << _line_width << ";stroke-miterlimit:6;stroke-dasharray:none;stroke-opacity:1;enable-background:new"
        << "\" d=\"" << write_svg_path(path) << "\"/>" << std::endl;
}

void Hair::writeStitchPoints(std::ostream& out, std::string color) {
    for (auto s : stitches) {
        for (auto p : s) {
            writeCircle(out, p, _stitch_point_radius, color);
        }
    }
}

std::string Hair::getColor(size_t pos, size_t numColors) {
    const std::string chars("0123456789abcdef");
    std::string result;
    result.reserve(6);

    std::uniform_int_distribution<size_t> distribution(0, chars.size()-1);
    for (size_t ii = 0; ii < 6; ++ii) {
        result.push_back(chars[distribution(generator)]);
    }
    return result;
}

void Hair::addStartStop() {
    for (auto & patch : stitches) {
        addStartStop(patch);
    }
}

void Hair::addStartStop(std::vector<Point>& patch) {
    if (patch.size() < 2) {
        return;
    }
    std::vector<Point> start(2);
    start[0] = patch[0];
    Point direction = patch[1] - patch[0];
    start[1] = patch[0] + 0.3*direction / direction.length();
    patch.insert(patch.begin(), start.begin(), start.end());

    const Point& last = patch.back();
    direction = patch[patch.size()-2] - last;
    patch.push_back(last + 0.3*direction / direction.length());
    patch.push_back(last);
}

template<class Line>
void Hair::writePatches(std::ostream& out, std::vector<Line>& patches) {
    for (size_t ii = 0; ii < patches.size(); ++ii) {
        std::vector<Point> & patch = patches[ii];
        if (patch.size() < 2) {
            continue;
        }
        std::string color = getColor(ii, patches.size());
        write(out, getPath(patch), color);
        writeCircles(out, patch, _stitch_point_radius, color);
    }
}

template<class C>
Path Hair::getPath(const std::vector<C>& points) {
    Path result;
    if (points.size() < 2) {
        return result;
    }
    for (size_t ii = 0; ii+1 < points.size(); ++ii) {
        BezierCurve straightLine = BezierCurveN<1>(points[ii], points[ii+1]);
        result.append(straightLine);
    }
    return result;
}

Point Hair::getCenter(const std::vector<Point>& stitches) {
    if (stitches.empty()) {
        return Point(0,0);
    }
    Point max = stitches.front();
    Point min = max;
    for (const auto p : stitches) {
        if (p.x() < min.x()) {
            min.x() = p.x();
        }
        if (p.y() < min.y()) {
            min.y() = p.y();
        }
        if (p.x() > max.x()) {
            max.x() = p.x();
        }
        if (p.y() > max.y()) {
            max.y() = p.y();
        }
    }
    return (min+max)/2;
}

void Hair::writeCircles(std::ostream& out, const std::vector<Point> & centers, double radius, std::string color) {
    for (auto center : centers) {
        writeCircle(out, center, radius, color);
    }
}

void Hair::writeCircle(std::ostream& out, Point center, double radius, std::string color) {
    out << "<circle"
        << " style=\"opacity:1;fill:#" << color << ";fill-opacity:1;stroke:none;\" "
        << " cx=\"" << center.x() << "\""
        << " cy=\"" << center.y() << "\""
        << " r=\"" << radius << "\" />" << std::endl;
}
