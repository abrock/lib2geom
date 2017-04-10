#include "lineoptimization.h"

#include <list>

using namespace Geom;

LineOptimization::LineOptimization() {}

LineOptimization::Solution LineOptimization::greedySolution(const size_t start, const bool start_reversed) {
    if (start >= lines.size()) {
        throw std::runtime_error("Starting index for the greedy solution is out of bounds");
    }
    calcNeighbours();
    Solution current_solution(*this);
    current_solution.addStep(Step(start, start_reversed, false), Geom::Point());
    Step last_step = current_solution.steps.back();
    Point last_stop = start_reversed ? lines[last_step.index].initialPoint() : lines[last_step.index].finalPoint();

    size_t counter = 0;
    while (!current_solution.isComplete()) {
        if (++counter > lines.size() * (1 + maxDraws)) {
            throw std::runtime_error("You're doing it wrong! Too many steps, something went south, aborting.");
        }
        std::vector<Step> options;
        options.reserve(size());

        size_t last_index = last_step.index;

        //continue_options.reserve(lines.size());
        std::vector<EndpointIndex> const& neighbours = last_step.reverse ?
                    lines[last_index].startNeighbours : lines[last_index].endNeighbours;

        for (EndpointIndex const& e : neighbours) {
            if (current_solution.mayDraw(e.index) && (changeDirection || e.start)) {
                options.push_back(Step(e.index, !e.start, false));
            }
        }
        if (options.empty()) {
            for (EndpointIndex const& e : neighbours) {
                if (current_solution.mayRedraw(e.index) && (changeDirection || e.start)) {
                    options.push_back(Step(e.index, !e.start, true));
                }
            }
        }

        if (options.empty()) {
            for (size_t ii = 0; ii < lines.size(); ++ii) {
                if (current_solution.mayDraw(ii)) {
                    options.push_back(Step(ii, false, false, true));
                    if (changeDirection) {
                        options.push_back(Step(ii, true, false, true));
                    }
                }
            }
        }

        if (options.empty()) {
            for (size_t ii = 0; ii < lines.size(); ++ii) {
                if (current_solution.mayRedraw(ii)) {
                    options.push_back(Step(ii, false, true, true));
                    if (changeDirection) {
                        options.push_back(Step(ii, true, true, true));
                    }
                }
            }
        }

        if (options.empty()) {
            throw std::runtime_error("No options for continuing found despite the solution not being complete, aborting.");
        }
        size_t dice = rand() % options.empty();
        current_solution.addStep(options[dice], last_stop);
        last_step = current_solution.steps.back();
        last_stop = last_step.reverse ? lines[last_step.index].initialPoint() : lines[last_step.index].finalPoint();
    }

    return current_solution;
}

LineOptimization::Solution LineOptimization::randomGreedySolution(size_t num_repeat) {
    Solution bestCost(*this), bestStops(*this), bestPareto(*this);
    bestCost = bestStops = bestPareto = greedySolution(rand() % size(), rand() % 2);
    std::cout << "best stops: " << bestStops.stopCount << std::endl;

    for (size_t ii = 0; ii < num_repeat; ++ii) {
        Solution current_sol = greedySolution(rand() % size(), rand() % 2);
        bool printed_id = false;
        if (current_sol.stopCount < bestStops.stopCount) {
            std::cout << "it# " << ii << std::endl;
            printed_id = true;
            std::cout << "new best stops: " << current_sol.stopCount << " (was: " << bestStops.stopCount << ")" << std::endl;
            bestStops = current_sol;
        }
        if (current_sol.costSum < bestCost.costSum) {
            if (!printed_id) {
                std::cout << "it# " << ii << std::endl;
                printed_id = true;
            }
            std::cout << "new best cost: " << current_sol.costSum << " (was: " << bestCost.costSum << ")" << std::endl;
            bestCost = current_sol;
        }
        if (current_sol.paretoBetterThan(bestPareto)) {
            if (!printed_id) {
                std::cout << "it# " << ii << std::endl;
                printed_id = true;
            }
            std::cout << "new best pareto: " << current_sol.costSum << "/" << current_sol.stopCount << " (was: "
                      << bestPareto.costSum << "/" << bestPareto.stopCount << std::endl;
            bestPareto = current_sol;
        }
    }

    return bestPareto;
}

void LineOptimization::addPaths(Geom::PathVector const& in) {
    for (Geom::Path const& c : in) {
        lines.push_back(c);
    }
}

void LineOptimization::splitAtForkEnds(std::vector<OptLine> &in) {

    bool global_did = true;
    size_t global_loops = 0;
    while (global_did && global_loops < 100) {
        global_loops++;
        std::cout << "it# " << global_loops << ": " << in.size() << std::endl;
        global_did = false;
        std::vector<OptLine> replacement;
        for (size_t ii = 0; ii < in.size(); ++ii) {
            Path const& A = in[ii];
            bool loop_1_did = false;
            for (size_t jj = 0; jj < in.size(); ++jj) {
                if (jj == ii) {
                    continue;
                }
                Path const& B = in[jj];
                bool did = false;
                if (!did) {
                    double distance = std::numeric_limits<double>::max();
                    PathTime const nearest = A.nearestTime(B.initialPoint(), &distance);
                    if (distance <= forkThresholdDistance) {
                        Point nearest_point = A.pointAt(nearest);
                        if (
                                Geom::distance(nearest_point, A.finalPoint()) > minimumSegmentLength
                                && Geom::distance(nearest_point, A.initialPoint()) > minimumSegmentLength) {
                            std::vector<PathTime> times;
                            times.push_back(nearest);
                            std::vector<Path> parts = A.subdivide(times);
                            for (Path const& p : parts) {
                                replacement.push_back(p);
                            }
                            if (parts.size() > 1) {
                                did = true;
                            }
                        }
                    }
                }
                if (!did) {
                    double distance = std::numeric_limits<double>::max();
                    PathTime const nearest = A.nearestTime(B.finalPoint(), &distance);
                    if (distance <= forkThresholdDistance) {
                        Point nearest_point = A.pointAt(nearest);
                        if (
                                Geom::distance(nearest_point, A.finalPoint()) > minimumSegmentLength
                                && Geom::distance(nearest_point, A.initialPoint()) > minimumSegmentLength) {
                            std::vector<PathTime> times;
                            times.push_back(nearest);
                            std::vector<Path> parts = A.subdivide(times);
                            for (Path const& p : parts) {
                                replacement.push_back(p);
                            }
                            if (parts.size() > 1) {
                                did = true;
                            }
                        }
                    }
                }
                if (did) {
                    global_did = true;
                    loop_1_did = true;
                    cacheValid = false;
                    break;
                }
            }
            if (!loop_1_did) {
                replacement.push_back(A);
            }
        }
        in = replacement;
    }
}

double calcLength(Path const& p, double const step = 0.1) {
    double sum = 0;
    for (Curve const& c : p) {
        Point last_point = c.initialPoint();

        for (double t = step; t < 1.0; t += step) {
            Geom::Point current_point = c.pointAt(t);
            sum += Geom::distance(last_point, current_point);
            last_point = current_point;
        }
        sum += Geom::distance(last_point, c.finalPoint());
    }
    return sum;
}

void LineOptimization::calcNeighbours()
{
    if (cacheValid) {
        return;
    }
    cacheValid = true;
    for (OptLine& l : lines) {
        l.startNeighbours.clear();
        l.endNeighbours.clear();
    }
    for (size_t a_ind = 0; a_ind < lines.size(); ++a_ind) {
        OptLine& A = lines[a_ind];
        for (size_t b_ind = a_ind+1; b_ind < lines.size(); ++b_ind) {
            OptLine& B = lines[b_ind];
            if (isNeighbour(A.finalPoint(), B.initialPoint())) {
                A.endNeighbours.push_back(EndpointIndex(b_ind, true));
                if (changeDirection) {
                    B.startNeighbours.push_back(EndpointIndex(a_ind, false));
                }
            }
            if (isNeighbour(A.initialPoint(), B.finalPoint())) {
                B.endNeighbours.push_back(EndpointIndex(a_ind, true));
                if (changeDirection) {
                    A.startNeighbours.push_back(EndpointIndex(b_ind, false));
                }
            }
            if (changeDirection) {
                if (isNeighbour(A.initialPoint(), B.initialPoint())) {
                    A.startNeighbours.push_back(EndpointIndex(b_ind, true));
                    B.startNeighbours.push_back(EndpointIndex(a_ind, true));
                }
                if (isNeighbour(A.finalPoint(), B.finalPoint())) {
                    A.endNeighbours.push_back(EndpointIndex(b_ind, false));
                    B.endNeighbours.push_back(EndpointIndex(a_ind, false));
                }
            }
        }
    }
}

bool LineOptimization::isNeighbour(const Point &a, const Point &b) {
    return Geom::distanceSq(a, b) <= forkThresholdDistance * forkThresholdDistance;
}

void LineOptimization::splitAtForks() {
    size_t old_count = lines.size();
    splitAtForkEnds(lines);
    size_t new_count = lines.size();
    std::cout << "Lines count change: " << old_count << " -> " << new_count << std::endl;
}

void LineOptimization::Solution::addStep(const LineOptimization::Step &step, Point const& prev) {
    steps.push_back(step);
    if (step.index >= drawCounts.size()) {
        drawCounts.resize(step.index+1, 0);
    }
    drawCounts[step.index] += step.redraw ? 1 : parent.drawFactor;
    costSum += (step.redraw ? 1 : parent.drawFactor) * parent.lines[step.index].getLength();
    if (step.jump) {
        stopCount++;
        costSum += parent.stopCost;
        costSum += parent.g0Cost * Geom::distance(prev, step.reverse ?
                                                      parent.lines[step.index].finalPoint() :
                                                  parent.lines[step.index].initialPoint());
    }
}

bool LineOptimization::Solution::isDrawn(const size_t index) const {
    return drawCounts[index] >= parent.minDraws;
}

bool LineOptimization::Solution::isComplete() const {
    for (size_t ii = 0; ii < parent.size(); ++ii) {
        if (!isDrawn(ii)) {
            return false;
        }
    }
    return true;
}

bool LineOptimization::Solution::mayDraw(size_t index) const {
    return drawCounts[index] + parent.drawFactor <= parent.maxDraws;
}

bool LineOptimization::Solution::mayRedraw(size_t index) const {
    return drawCounts[index] + 1 <= parent.maxDraws;
}

void LineOptimization::Solution::printStats() const {
    size_t mayDrawCount = 0;
    size_t isDrawnCount = 0;
    size_t mayRedrawCount = 0;
    for (size_t ii = 0; ii < parent.size(); ++ii) {
        if (mayDraw(ii)) mayDrawCount++;
        if (mayRedraw(ii)) mayRedrawCount++;
        if (isDrawn(ii)) isDrawnCount++;
    }
    std::cout << "Sol stats: "
              << "size: " << parent.size() << ", "
              << steps.size() << " steps, "
              << mayDrawCount << " mayDraw, "
              << mayRedrawCount << " mayRedraw, "
              << isDrawnCount << " are done " << std::endl;
}

bool LineOptimization::Solution::paretoBetterThan(const LineOptimization::Solution &other) const {
    return (costSum <  other.costSum && stopCount <= other.stopCount) ||
           (costSum <= other.costSum && stopCount <  other.stopCount);
}

double LineOptimization::OptLine::getLength() const {
    if (cacheValid) {
        return length;
    }
    cacheValid = true;
    length = calcLength(*this);

    return length;
}
