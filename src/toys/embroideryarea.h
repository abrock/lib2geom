#ifndef EMBROIDERYPATCH_H
#define EMBROIDERYPATCH_H
#include <vector>
#include "embroideryline.h"

class Hair;

/**
 * @brief The EmbroideryArea class represents a sub-area of a larger filled embroidery area
 * which can be filled without using connecting sub-stitches.
 */
class EmbroideryArea : public std::vector<EmbroideryLine> {


public:
    size_t min_level;
    std::vector<Point> forward_stitches;
    std::vector<Point> reverse_stitches;
    OutlineIntersection start;
    OutlineIntersection stop;
    OutlineIntersection r_start;
    OutlineIntersection r_stop;

    EmbroideryArea(EmbroideryLine first_line, const size_t _min_level)
        : min_level(_min_level)
        , start(first_line.startInter)
        , stop(first_line.endInter)
        , r_start(first_line.endInter)
        , r_stop (first_line.startInter)
    {
        push_back(first_line);
    }

    template<class Iter>
    void addPoints(std::vector<Point>& container, Iter start, Iter stop) {
        for(; start != stop; start++) {
            container.push_back(Point(start->x(), start->y()));
        }
    }

    friend bool operator< (const EmbroideryArea &a, const EmbroideryArea &b) {
        if (a.min_level < b.min_level) {
            return true;
        }
        if (a.min_level > b.min_level) {
            return false;
        }
        if (a.forward_stitches.size() < b.forward_stitches.size()) {
            return true;
        }
        if (a.forward_stitches.size() > b.forward_stitches.size()) {
            return false;
        }
        return a.size() < b.size();
    }

    /**
     * @brief finish Create the stitch lists for forward and reverse
     */
    void finish(const Hair& helper);
};

#endif // EMBROIDERYPATCH_H
