#include "embroideryarea.h"

#include "hair.h"

void EmbroideryArea::finish(const Hair& helper) {
    addPoints(forward_stitches, (*this)[0].begin(), (*this)[0].end());
    addPoints(reverse_stitches, (*this)[0].rbegin(), (*this)[0].rend());

    for (size_t ii = 1; ii < size(); ++ii) {
        Point current_forward_stop = Point(stop.x(), stop.y());
        Point current_forward_last = forward_stitches.back();
        Point current_reverse_stop = Point(r_stop.x(), r_stop.y());
        Point current_reverse_last = reverse_stitches.back();
        const double forward_error = (current_forward_stop - current_forward_last).length();
        const double reverse_error = (current_reverse_stop - current_reverse_last).length();
        assert(forward_error < 1e-4);
        assert(reverse_error < 1e-4);

        const EmbroideryLine& current_line = (*this)[ii];
        std::vector<Point> forward_forward = helper.getShortestConnection(stop, current_line.startInter);
        std::vector<Point> forward_reverse = helper.getShortestConnection(stop, current_line.endInter);
        if (forward_forward.size() < forward_reverse.size()) {
            addPoints(forward_stitches, forward_forward.begin(), forward_forward.end());
            addPoints(forward_stitches, current_line.begin(), current_line.end());
            stop = current_line.endInter;

            std::vector<Point> reverse_bridge = helper.getShortestConnection(r_stop, current_line.endInter);
            addPoints(reverse_stitches, reverse_bridge.begin(), reverse_bridge.end());
            addPoints(reverse_stitches, current_line.rbegin(), current_line.rend());
            r_stop = current_line.startInter;
        } else {
            addPoints(forward_stitches, forward_reverse.begin(), forward_reverse.end());
            addPoints(forward_stitches, current_line.rbegin(), current_line.rend());
            stop = current_line.startInter;

            std::vector<Point> reverse_bridge = helper.getShortestConnection(r_stop, current_line.startInter);
            addPoints(reverse_stitches, reverse_bridge.begin(), reverse_bridge.end());
            addPoints(reverse_stitches, current_line.begin(), current_line.end());
            r_stop = current_line.endInter;
        }

    }
}
