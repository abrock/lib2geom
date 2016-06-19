#include "hair.h"

/**
 * @brief writeStitches Write the stitches into a text file suitable for libembroidery-convert
 * @param filename
 */
void Hair::writeStitches(const char* filename) {
    std::ofstream out(filename);
    writeStitches(out);
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

Hair::Hair(Path _outline, Path _curve) : outline(_outline), curve(_curve) {}

void Hair::run() {
    clock_t start = clock();
    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;

    //makeAreaLarger(outline, areaGrow);

    clock_t start2 = clock();
    getBoundaryDiscretization();
    clock_t stop2 = clock();
    std::cerr << "getBoundaryDiscretization took " << (static_cast<double>(stop2-start2)) / CLOCKS_PER_SEC << std::endl << std::endl;
    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;

    start2 = clock();
    getCurves();
    stop2 = clock();
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
    assembleStitches();
    stop2 = clock();
    std::cerr << "assembleStitches took " << (static_cast<double>(stop2-start2)) / CLOCKS_PER_SEC << std::endl << std::endl;
    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;

    std::cerr << "#stitches: " << countStitches() << std::endl;
    start2 = clock();
    for (auto & patch: stitches) {
        patch = purgeSmallStitches(patch, 0.1);
    }
    stop2 = clock();
    std::cerr << "purgeSmallStitches took " << (static_cast<double>(stop2-start2)) / CLOCKS_PER_SEC << std::endl << std::endl;
    std::cerr << "using " << memoryConsumptionKB()/1024 << "MB of memory" << std::endl;

    addStartStop();
    std::cerr << "#stitches after purge: " << countStitches() << std::endl;


    clock_t stop = clock();
    std::cerr << "Calculation took " << (static_cast<double>(stop-start)) / CLOCKS_PER_SEC << std::endl << std::endl;
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
    bridges = bestBridges;
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
    const size_t discreteSize = discreteOutline.size();

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

            const size_t lastPointOutline = bestMatch(discreteOutline, lastPoint);
            const size_t frontOption = bestMatch(discreteOutline, linePoints.front());
            const size_t backOption = bestMatch(discreteOutline, linePoints.back());
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
                bridge.push_back(discreteOutline[static_cast<size_t>(ii)]);
                if (bridge.size() > maxBridgeSize) {
                    std::cerr << "Bridge has grown too large, maximum size is " << maxBridgeSize << ", lasPointOutline is " << lastPointOutline << ", option is " << option << ", direction is " << direction << std::endl;
                    break;
                }
            }
            bridge.push_back(discreteOutline[option]);
            // We remove unnecessary stitches from the bridge.
            const Point& nextLineStart = reverseCurrentLine ? linePoints.back() : linePoints.front();
            if (bridge.size() > 1) {
                const Point& lastBridge = bridge.back();
                const Point& lastlastBridge = bridge[bridge.size()-2];
                if (dot(lastlastBridge - lastBridge, nextLineStart - lastBridge) > 0
                        && (lastlastBridge - nextLineStart).length() < stitchLength/2) {
                    bridge.pop_back();
                }
            }
            if (bridge.size() > 1) {
                const Point& firstBridge = bridge.front();
                const Point& secondBridge = bridge[1];
                if (dot(secondBridge - firstBridge, lastPoint - firstBridge) > 0
                        && (secondBridge - lastPoint).length() < stitchLength/2) {
                    bridge.erase(bridge.begin());
                }
            }
            if (bridge.size() == 1) {
                const Point& bridgeElement = bridge.front();
                if (dot(lastPoint - bridgeElement, nextLineStart - bridgeElement) > 0
                        || (lastPoint - nextLineStart).length() < stitchLength/3) {
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
        if (stitchLength >= (p - lines[lastLine+1].front()).length()
                || stitchLength >=(p - lines[lastLine+1].back()).length()) {
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
    RunningStats stitchLength;
    std::vector<std::vector<Point> > newStitches;
    levels = std::vector<EmbroideryLineLevel>(stitches.size());
#pragma omp parallel for
    for (size_t level = 0; level < stitches.size(); ++level) {
        const std::vector<Point>& line = stitches[level];
        if (line.size() < 2) {
            continue;
        }
        std::vector<Point> newLine;
        std::vector<BOOL> inside(line.size());
        for (size_t ii = 0; ii < line.size(); ++ii) {
            inside[ii] = isInsidePoint(line[ii]);
        }
        if (inside[0]) {
            newLine.push_back(line[0]);
        }
        OutlineIntersection startInter(outline), endInter(outline);
        for (size_t ii = 0; ii + 1 < line.size(); ++ii) {
            if (inside[ii] && inside[ii+1]) {
                newLine.push_back(line[ii+1]);
                stitchLength.push((line[ii] - line[ii+1]).length());
                continue;
            }
            BezierCurve _straightLine = BezierCurveN<1>(line[ii], line[ii+1]);
            Path straightLine;
            straightLine.append(_straightLine);
            std::vector<PathIntersection> intersection = straightLine.intersect(outline);
            if (inside[ii] && !inside[ii+1]) {
                //straightLine.stitchTo(line[ii+1]);
                //std::cerr << "straightLine: " << write_svg_path(straightLine) << std::endl;
                Point intersectPoint = straightLine.pointAt(intersection[0].first);
                endInter = OutlineIntersection(outline, intersection[0].second, level);
                newLine.push_back(intersectPoint);
                if (lineLength(newLine) > minStitchLength) {
#pragma omp critical
                    {
                        newStitches.push_back(newLine);
                        levels[level].push_back(EmbroideryLine(newLine, level, startInter, endInter));
                        intersections.push_back(startInter);
                        intersections.push_back(endInter);
                    }
                }
                newLine.clear();
            }
            if (!inside[ii] && inside[ii+1]) {
                Path straightLine;
                straightLine.append(_straightLine);
                //straightLine.stitchTo(line[ii+1]);
                //std::cerr << "straightLine: " << write_svg_path(straightLine) << std::endl;
                Point intersectPoint = straightLine.pointAt(intersection[0].first);
                newLine.push_back(intersectPoint);
                newLine.push_back(line[ii+1]);
                startInter = OutlineIntersection(outline, intersection[0].second, level);
            }
            if (!inside[ii] && !inside[ii+1]) {
                if (2 <= intersection.size()) {
                    const Point startPoint = straightLine.pointAt(intersection[0].first);
                    const Point endPoint = straightLine.pointAt(intersection[1].first);
                    newLine.push_back(startPoint);
                    newLine.push_back(endPoint);
                    startInter = OutlineIntersection(outline, intersection[0].second, level);
                    endInter = OutlineIntersection(outline, intersection[1].second, level);
                    if (lineLength(newLine) > minStitchLength) {
#pragma omp critical
                        {
                            newStitches.push_back(newLine);
                            levels[level].push_back(EmbroideryLine(newLine, level, startInter, endInter));
                            intersections.push_back(startInter);
                            intersections.push_back(endInter);
                        }
                    }
                    newLine.clear();
                }
            }
        }
    }
    stitches = newStitches;
    std::cout << "Number of intersections: " << intersections.size() << std::endl;
    std::sort(intersections.begin(), intersections.end());
    std::cerr << "Stitch lengh stats: " << stitchLength.print() << std::endl;
}

double Hair::lineLength(const std::vector<Point>& points) {
    double result = 0;
    for (size_t ii = 1; ii < points.size(); ++ii) {
        result += (points[ii-1] - points[ii]).length();
    }
    return result;
}


void Hair::getStitches() {

    const size_t initialIndex = curves.size()/2;
    Path initialCurve = curves[initialIndex];

    stitches.assign(curves.size(), std::vector<Point>());
    getStitches(stitches[initialIndex], initialCurve);

    //enlargeCurves();

    initialStitchLengthsSum = getStitchLengthes(stitches[initialIndex], initialStichLengths);

    for (size_t ii = initialIndex+1; ii < curves.size(); ++ii) {
        auto tmp = projection(stitches[ii-1], curves[ii]);
        if (!atLeastOneInside(tmp, outline)) {
            break;
        }
        stitches[ii] = tmp;
    }
    for (size_t ii = initialIndex; ii+1 > 0; --ii) {
        auto tmp = projection(stitches[ii+1], curves[ii]);
        if (!atLeastOneInside(tmp, outline)) {
            break;
        }
        stitches[ii] = tmp;
    }
    //        for (size_t ii = 0; ii < curves.size(); ii += 2) {
    //            stitches[ii] = vectorOffset(stitches[ii], offset);
    //        }
    double currentOffset = 0;
    for (size_t ii = 0; ii < curves.size(); ++ii) {
        stitches[ii] = vectorOffsetOnCurve(stitches[ii], curves[ii], currentOffset);
        currentOffset = fmod(currentOffset + offset, 1.0);
    }

}

void Hair::enlargeCurves() {
    for (auto &curve : curves) {
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
std::vector<Point> Hair::vectorOffset(
        const std::vector<Point>& points,
        const double offset) {
    std::vector<Point> result = points;
    if (points.size() < 2) {
        return result;
    }
    result[0] += (points[1] - points[0]) * offset;
    for (size_t ii = 1; ii < points.size(); ++ii) {
        result[ii] += (points[ii] - points[ii-1])*offset;
    }
    return result;
}

/**
 * @brief vectorOffset moves every point in a vector in the direction given by the difference between itself and the previous point.
 * @param points The original points
 * @param offset The relative offset, usually between -1 and 1
 * @return The moved points.
 */
std::vector<Point> Hair::vectorOffsetOnCurve(
        const std::vector<Point>& points,
        const Path& curve,
        const double offset) {
    std::vector<Point> result = points;
    if (points.size() < 2) {
        return result;
    }
    Point lastPoint = points[1];
    for (size_t ii = 0; ii < points.size(); ++ii) {
        PathTime t = curve.nearestTime(result[ii]);
        const double distance = (lastPoint - points[ii]).length();
        result[ii] = getOffsettedPointOnCurve(curve, t, offset * distance);

        lastPoint = points[ii];
    }

    return result;
}

bool Hair::isInsidePoint(const Point p) {
    return 0 != outline.winding(p);
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
    getStitches(initialCurvePoints, initialCurve, stitchLength);
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
    std::vector<Point> data = curve.at(t.curve_index).pointAndDerivatives(t.t, 1);

    Point initialGuess = data[0] + data[1] * targetOffset / data[1].length();

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
}

void Hair::makeAreaLarger(Path& curve, const double offset) {
    curve = Inkscape::half_outline(curve, offset, miter);
}

void Hair::getCurves() {

    Path left = curve;
    Path right = curve.reversed();
    PathVector rights, lefts;

#pragma omp parallel sections
    {
        for (size_t ii = 0; ii < maxIter; ++ii) {
            left = Inkscape::half_outline(left, lineSpacing, miter);
            if (left.intersect(outline).empty()) {
                break;
            }
            lefts.push_back(left);
        }
#pragma omp section
        {
            for (size_t ii = 0; ii < maxIter; ++ii) {
                right = Inkscape::half_outline(right, lineSpacing, miter);
                if (right.intersect(outline).empty()) {
                    break;
                }
                rights.push_back(right);
            }
            rights.reverse(true);
            curves.insert(curves.begin(), rights.begin(), rights.end());
            curves.push_back(curve);
        }
    }

    curves.insert(curves.end(), lefts.begin(), lefts.end());
    std::cerr << "We now have " << curves.size() << " curves" << std::endl;
}

void Hair::write(const char* filename) {
    std::ofstream out(filename);
    write(out);
}

void Hair::write(std::ostream& out) {
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
    out << "<g>";
    for (auto c : curves) {
        write(out, c, "aaaaaa");
    }
    out << "</g>";
    write(out, curve, "ff0000");
    write(out, outline, "0000ff");

    //writeStitchPoints(out, "00ff00");
    out << "<g>";
    writePatches(out, stitches);
    out << "</g>";

    out << "<g>";
    write(out, getPath(discreteOutline), "ff0000");
    writeCircles(out, discreteOutline, stitchPointRadius, "ff0000");
    out << "</g>";

    out << "<g>";
    writePatches(out, bridges);
    out << "</g>";
    out << "<g>";
    writeOutlineIntersections(out, intersections);
    out << "</g>";
    out << "</svg>" << std::endl;
}


void Hair::writeOutlineIntersections(std::ostream& out, std::vector<OutlineIntersection> &intersections) {
    std::vector<Point> points;
    for (auto inter: intersections) {
        points.push_back(inter.outline.pointAt(inter.time));
    }
    const auto path = getPath(points);
    write(out, getPath(points), "ff0000");

}

void Hair::write(std::ostream& out, Path path, std::string color) {
    out << "<path style=\"display:inline;fill:none;fill-opacity:1;stroke:#"
        << color
        << ";stroke-width:" << lineWidth << ";stroke-miterlimit:6;stroke-dasharray:none;stroke-opacity:1;enable-background:new"
        << "\" d=\"" << write_svg_path(path) << "\"/>" << std::endl;
}

void Hair::writeStitchPoints(std::ostream& out, std::string color) {
    for (auto s : stitches) {
        for (auto p : s) {
            writeCircle(out, p, stitchPointRadius, color);
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
        writeCircles(out, patch, stitchPointRadius, color);
    }
}

Path Hair::getPath(std::vector<Point>& points) {
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

void Hair::writeCircles(std::ostream& out, std::vector<Point> & centers, double radius, std::string color) {
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
