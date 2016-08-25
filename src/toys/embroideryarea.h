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
    const size_t min_level;
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

    /**
     * @brief finish Create the stitch lists for forward and reverse
     */
    void finish(const Hair& helper) {
        addPoints(forward_stitches, (*this)[0].begin(), (*this)[0].end());
        addPoints(reverse_stitches, (*this)[0].rbegin(), (*this)[0].rend());

        for (size_t ii = 1; ii < size(); ++ii) {
            Point current_forward_stop = stop.point();
            Point current_forward_last = forward_stitches.back();
            Point current_reverse_stop = r_stop.point();
            Point current_reverse_last = reverse_stitches.back();
            const double forward_error = (current_forward_stop - current_forward_last).length();
            const double reverse_error = (current_reverse_stop - current_reverse_last).length();
            assert(forward_error < 1e-4);
            assert(reverse_error < 1e-4);
        }
    }
};

#endif // EMBROIDERYPATCH_H
