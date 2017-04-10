#ifndef LINEOPTIMIZATION_H
#define LINEOPTIMIZATION_H

#include <2geom/path.h>
#include <2geom/pathvector.h>

using namespace Geom;

class LineOptimization
{
public:
    /**
     * @brief The EndpointIndex class stores an index corresponding to a Path in the lines vector of the LineOptimization class
     * and a boolean (start) which stores the information if the EnpointIndex corresponds to the initial point (true) or endpoint (false).
     */
    class EndpointIndex {
    public:
        size_t index;
        bool start;
        EndpointIndex(size_t _index = 0, bool _start = true) : index(_index), start(_start) {}
    };

    /**
     * @brief The OptLine class stores a Geom::Path together with information about the Path and its neightbours in the
     * lines vector of the LineOptimization class.
     */
    class OptLine : public Path {
    private:
        mutable bool cacheValid = false;
        mutable double length = 0;
    public:
        OptLine(Path const p) : Path(p) {}
        /**
         * @brief startNeighbours List of Paths which have an endpoint close to the current Path's initial point.
         */
        std::vector<EndpointIndex> startNeighbours;
        /**
         * @brief startNeighbours List of Paths which have an endpoint close to the current Path's end point.
         */
        std::vector<EndpointIndex> endNeighbours;

        double getLength() const;
    };

    /**
     * @brief The Step class stores information about one step in an solution to the optimization problem.
     */
    class Step {
    public:
        /**
         * @brief index Index of the Path which is drawn in the current step.
         */
        size_t index;
        /**
         * @brief reverse True if the Path is traversed in the reverse direction.
         */
        bool reverse;
        /**
         * @brief redraw True if the Path is "redrawn" which might or might not differ from drawing it.
         * For example in embroidery one could "draw" with a triple stitch and "redraw" with a normal stitch.
         */
        bool redraw;

        /**
         * @brief jump True if the step requires stopping the machine, for example turning of the laser
         * or trimming the thread in embroidery.
         */
        bool jump;
        Step(size_t _index = 0, bool _reverse = false, bool _redraw = false, bool _jump = false) :
            index(_index), reverse(_reverse), redraw(_redraw), jump(_jump) {}

        Step(Step const& other) :
            index(other.index),
            reverse(other.reverse),
            redraw(other.redraw),
            jump(other.jump)
        {}
    };

    class Solution {
    public:
        typedef char BOOL;

        LineOptimization & parent;
        std::vector<size_t> drawCounts;
        std::vector<Step> steps;

        double costSum = 0;
        size_t stopCount = 0;

        Solution(LineOptimization & _parent) :
            parent(_parent),
            drawCounts(_parent.size(), 0) {
            steps.reserve(parent.size());
        }

        Solution(Solution const& other) :
            parent(other.parent),
            drawCounts(other.drawCounts),
            steps(other.steps),
            costSum(other.costSum),
            stopCount(other.stopCount)
        {
            steps.reserve(other.parent.size());
        }

        Solution operator = (Solution const& other) {
            parent = other.parent;
            drawCounts = other.drawCounts;
            steps = other.steps;
            costSum = other.costSum;
            stopCount = other.stopCount;
            return *this;
        }

        void addStep(size_t index, bool reversed = false, bool redrawn = false, bool jump = false);
        void addStep(Step const& step, Geom::Point const& prev);

        bool isDrawn(size_t const index) const;
        bool isComplete() const;

        bool mayDraw(size_t index) const;
        bool mayRedraw(size_t index) const;

        void printStats() const;

        bool paretoBetterThan (Solution const& other) const;
    };

private:
    double bestCosts = std::numeric_limits<double>::max();

    /**
     * @brief drawingCost Estimated cost for re-drawing 1 length unit.
     */
    double redrawingCost = 1.0/(10*2);

    /**
     * @brief g0Cost Cost for moving the machine 1 length unit without working.
     */
    double g0Cost = 1.0/10;

    /**
     * @brief stopCost Estimated fixed cost for stopping at a point where no direct neightbour exists for
     * continuing. Applies to embroidery where the thread has to be trimmed but usually not to engraving/cutting.
     */
    double stopCost = 2;

    /**
     * @brief connectByRedraw
     */
    bool connectByRedraw = true;

    /**
     * @brief changeDirection Set to "true" for allowing the optimizer to reverse the direction of the line.
     */
    bool changeDirection = true;

    /**
     * @brief lines Set of lines for which the order will be optimized.
     */
    std::vector<OptLine> lines;

    /**
     * @brief forkThresholdDistance Maximum distance between two curves s.t. we consider them to have an intersection.
     */
    double forkThresholdDistance = 1;

    /**
     * @brief minimumSegmentLength Minimum resulting segment length so we consider splitting a Path into sub-paths.
     */
    double minimumSegmentLength = 1;



    void splitAtForkEnds(std::vector<OptLine>& in);

    bool cacheValid = false;

    void calcNeighbours();

    bool isNeighbour(Geom::Point const& a, Geom::Point const& b);

    /**
     * @brief maxDraws This is the maximum number of times a line may be drawn.
     * The optimizer won't produce solutions where the number of (re-)drawings of a line exceeds
     * this value. Redrawings are counted as 1, drawings as @see drawFactor.
     */
    size_t maxDraws = 4;

    /**
     * @brief minDraws This is the minimum number of times a line must be drawn.
     * The optimizer won't produce solutions where the number of (re-)drawings is lower than this value.
     * Redrawings are counted as 1, drawings as @see drawFactor.
     */
    size_t minDraws = 3;

    /**
     * @brief drawFactor This is the cost factor between drawing a line and re-drawing a line.
     * 3 is a typical value, e.g. for triple stitch vs. running stitch in embroidery.
     */
    size_t drawFactor = 3;

public:
    LineOptimization();

    Solution greedySolution(size_t const start = 0, bool const start_reversed = false);
    Solution randomGreedySolution(size_t num_repeat = 1e3);

    size_t size() const {return lines.size();}

    /**
     * @brief addPaths Add a set of paths to the
     * @param in
     */
    void addPaths(Geom::PathVector const& in);

    void splitAtForks();


};

#endif // LINEOPTIMIZATION_H
