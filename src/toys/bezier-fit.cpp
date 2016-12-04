#include <random>

#include <2geom/d2.h>
#include <2geom/sbasis.h>
#include <2geom/bezier-to-sbasis.h>
#include <2geom/sbasis-geometric.h>

#include <toys/path-cairo.h>
#include <toys/toy-framework-2.h>
#include "path-cairo.h"
#include <cairo/cairo.h>
#include <2geom/svg-path-writer.h>
#include <2geom/bezier-utils.h>

#include "../tests/bezierfit.h"

namespace experiment {
#include "../tests/bezierfit-a.h"
}

#include <vector>
using std::vector;
using namespace Geom;

static void dot_plot(cairo_t *cr, Piecewise<D2<SBasis> > const &M, double space=10){
    //double dt=(M[0].cuts.back()-M[0].cuts.front())/space;
    Piecewise<D2<SBasis> > Mperp = rot90(derivative(M)) * 2;
    for( double t = M.cuts.front(); t < M.cuts.back(); t += space) {
        Point pos = M(t), perp = Mperp(t);
        draw_line_seg(cr, pos + perp, pos - perp);
    }
    cairo_pw_d2_sb(cr, M);
    cairo_stroke(cr);
}

static void cross_plot(
        cairo_t *cr,
        std::vector<Geom::Point> const& data,
        Geom::Point offset = Geom::Point(0,0),
        double const size = 10){
    for (auto point : data) {
        Geom::Point dx(size/2, 0);
        Geom::Point dy(0, size/2);
        draw_line_seg(cr, offset + point + dx, offset + point - dx);
        draw_line_seg(cr, offset + point + dy, offset + point - dy);
    }
}

#define SIZE 4

class BezierFitTester: public Toy {
public:
    PointSetHandle b_handle;
    void draw(cairo_t *cr,
              std::ostringstream *notify,
              int width, int height, bool save, std::ostringstream *timer_stream) {

        if (first_time)
        {
            first_time = false;
            sliders[0].geometry(Point(50, 50), 100);
        }
        size_t const num_points = static_cast<size_t>(sliders[0].value());

        D2<SBasis> B1 = b_handle.asBezier();
        Piecewise<D2<SBasis> >B;
        B.concat(Piecewise<D2<SBasis> >(B1));

        // testing fuse_nearby_ends
        std::vector< Piecewise<D2<SBasis> > > pieces;
        pieces = fuse_nearby_ends(split_at_discontinuities(B),9);
        Piecewise<D2<SBasis> > C;
        for (unsigned i=0; i<pieces.size(); i++){
            C.concat(pieces[i]);
        }
        // testing fuse_nearby_ends

        cairo_set_line_width (cr, 2.);
        cairo_set_source_rgba (cr, 0., 0.5, 0., 1);
        //cairo_d2_sb(cr, B1);
        //cairo_pw_d2_sb(cr, C);
        //cairo_pw_d2_sb(cr, B);
        cairo_stroke(cr);

        Timer tm;
        Timer::Time als_time = tm.lap();

        cairo_set_source_rgba (cr, 0., 0., 0.9, 1);
        //dot_plot(cr,uniform_B);
        cairo_stroke(cr);

        std::cout << B[0] << std::endl;

        Geom::Affine translation;

        Geom::Path original_path;
        //original_bezier.append(B[0]);
        //original_bezier.appendNew<CubicBezier> (B[0]);
        CubicBezier original_bezier(b_handle.pts);
        original_path.append(original_bezier);

        std::vector<double> initial_t;
        std::vector<Geom::Point> curve_points;
        if (randomize_times) {
            std::uniform_real_distribution<double> dist_t(0,1);
            for (size_t ii = 0; ii < num_points; ++ii) {
                double const t = dist_t(generator);
                initial_t.push_back(t);
            }
            std::sort(initial_t.begin(), initial_t.end());
            double const min = initial_t.front();
            double const max = initial_t.back();
            for (auto& t : initial_t) {
                t = (t-min)/(max-min);
            }
            for (auto const t : initial_t) {
                curve_points.push_back(original_bezier.pointAt(t));
            }
        }
        else {
            for (size_t ii = 0; ii < num_points; ++ii) {
                double const t = static_cast<double>(ii) / (num_points-1);
                Geom::Point const p = original_bezier.pointAt(t);
                initial_t.push_back(t);
                curve_points.push_back(p);
            }
        }

        cairo_set_source_rgba (cr, 0., 0., .9, 1);
        cairo_path(cr, original_path);
        draw_text(cr, original_path.initialPoint(), "original curve and old fit");


        Geom::CubicBezier fitted_new;
        Geom::CubicBezier fitted_new_a;
        Geom::Point very_old_version_raw[4];
        bezier_fit_cubic(very_old_version_raw, curve_points.data(), curve_points.size(), 0.);
        Geom::CubicBezier very_old_bezier(
                    very_old_version_raw[0],
                very_old_version_raw[1],
                very_old_version_raw[2],
                very_old_version_raw[3]
                );

        Geom::Path very_old_version_path;
        very_old_version_path.append(very_old_bezier);

        cairo_set_source_rgba (cr, .7, .7, 0., 1);
        cairo_stroke(cr);
        cairo_path(cr, very_old_version_path);

        cairo_set_source_rgba (cr, 0., 0., .9, 1);
        cairo_stroke(cr);
        cross_plot(cr, curve_points);

        if(1) {
            Geom::CubicBezier combination(very_old_bezier);
            tm.ask_for_timeslice();
            tm.start();
            auto new_result_ig_a = experiment::fit_bezier(combination, curve_points);
            als_time = tm.lap();
            *notify << "Bezier fit a, old algorithm as initial guess, time = " << als_time << std::endl
                    << "Worst residual: " << new_result_ig_a.first << " at t=" << new_result_ig_a.second << std::endl;

            Geom::Path combination_path;
            translation.setTranslation(Geom::Point(300,300));
            combination_path.append(combination.transformed(translation));

            cairo_set_source_rgba (cr, .0, .0, .9, 1);
            cross_plot(cr, curve_points, translation.translation());
            cairo_path(cr, combination_path);
            draw_text(cr, combination_path.initialPoint(), "old fit as i.g.");
        }
        {
            tm.ask_for_timeslice();
            tm.start();
            auto new_result = fit_bezier(fitted_new, curve_points);
            als_time = tm.lap();
            *notify << "Bezier fit, time = " << als_time << std::endl
                    << "Worst residual: " << new_result.first << " at t=" << new_result.second << std::endl;


            Geom::Path fitted_new_path;
            translation.setTranslation(Geom::Point(300,0));
            fitted_new_path.append(fitted_new.transformed(translation));

            cairo_set_source_rgba (cr, .0, .9, .0, 1);
            cross_plot(cr, curve_points, translation.translation());
            cairo_path(cr, fitted_new_path);
            draw_text(cr, fitted_new_path.initialPoint(), "new fit");
        }

        {
            tm.ask_for_timeslice();
            tm.start();
            auto new_result_a = experiment::fit_bezier(fitted_new_a, curve_points);
            als_time = tm.lap();
            *notify << "Bezier fit a, time = " << als_time << std::endl
                    << "Worst residual: " << new_result_a.first << " at t=" << new_result_a.second << std::endl;




            Geom::Path fitted_new_a_path;
            translation.setTranslation(Geom::Point(0,300));
            fitted_new_a_path.append(fitted_new_a.transformed(translation));


            cairo_set_source_rgba (cr, .9, .0, .0, 1);
            cross_plot(cr, curve_points, translation.translation());
            cairo_path(cr, fitted_new_a_path);
            draw_text(cr, fitted_new_a_path.initialPoint(), "new fit (a)");
        }

        Geom::CubicBezier fixed_times_bezier;
        {
            tm.ask_for_timeslice();
            tm.start();
            auto fixed_times_result = experiment::fit_bezier_fixed_times(fixed_times_bezier, curve_points);
            als_time = tm.lap();
            *notify << "Bezier fit a (fixed times), time = " << als_time << std::endl
                    << "Worst residual: " << fixed_times_result.first << " at t=" << fixed_times_result.second << std::endl;

            Geom::Path fixed_times_path;
            translation.setTranslation(Geom::Point(600,300));
            fixed_times_path.append(fixed_times_bezier.transformed(translation));

            cairo_set_source_rgba (cr, .9, .0, .0, 1);
            cross_plot(cr, curve_points, translation.translation());
            cairo_path(cr, fixed_times_path);
            draw_text(cr, fixed_times_path.initialPoint(), "fixed t fit (a)");
        }

        Geom::CubicBezier fixed_times_ig_bezier = fixed_times_bezier;
        {
            tm.ask_for_timeslice();
            tm.start();
            auto fixed_times_ig_result = experiment::fit_bezier(fixed_times_ig_bezier, curve_points);
            als_time = tm.lap();
            *notify << "Bezier fit a (with fixed times as i.g.), time = " << als_time << std::endl
                    << "Worst residual: " << fixed_times_ig_result.first << " at t=" << fixed_times_ig_result.second << std::endl;

            Geom::Path fixed_times_path;
            translation.setTranslation(Geom::Point(900,300));
            fixed_times_path.append(fixed_times_ig_bezier.transformed(translation));

            cairo_set_source_rgba (cr, .9, .0, .0, 1);
            cross_plot(cr, curve_points, translation.translation());
            cairo_path(cr, fixed_times_path);
            draw_text(cr, fixed_times_path.initialPoint(), "new (a) with fixed t as i.g.");
        }

        std::cout << "original: " << write_svg_path(original_path) << std::endl;

        Geom::CubicBezier initial_guess(
                    curve_points.front(), curve_points.front(),
                    curve_points.back(), curve_points.back()
                    );
        {
            experiment::get_initial_guess(initial_guess, curve_points);

            Geom::Path initial_guess_path;
            translation.setTranslation(Geom::Point(600,0));
            initial_guess_path.append(initial_guess.transformed(translation));

            cairo_set_source_rgba (cr, .8, .0, .8, 1);
            cross_plot(cr, curve_points, translation.translation());
            cairo_path(cr, initial_guess_path);
            draw_text(cr, initial_guess_path.initialPoint(), "initial guess");
        }

        cairo_stroke(cr);

        Toy::draw(cr, notify, width, height, save,timer_stream);
    }

public:
    BezierFitTester(){
        for(int i = 0; i < SIZE; i++) {
            b_handle.push_back(150+uniform()*300,150+uniform()*300);
        }
        b_handle.pts[0] = Geom::Point(70,250);
        b_handle.pts[1] = Geom::Point(200,150);
        b_handle.pts[2] = Geom::Point(200,350);
        b_handle.pts[3] = Geom::Point(350,200);
        handles.push_back(&b_handle);
        // M 70 250 C 860 766 200 350 350 200
        // M 70 250 C 906 833 200 350 350 200
        // M 70 250 C 800 738 200 350 350 200
        sliders.push_back(Slider(2, 50, 1, 25, "number of points"));
        handles.push_back(&(sliders[0]));
    }
private:
    std::vector<Slider> sliders;
    bool first_time = true;
    bool randomize_times = false;
    std::random_device rd;
    std::default_random_engine generator;
};

int main(int argc, char **argv) {
    init(argc, argv, new BezierFitTester);
    return 0;
}

/*
  Local Variables:
  mode:c++
  c-file-style:"stroustrup"
  c-file-offsets:((innamespace . 0)(inline-open . 0)(case-label . +))
  indent-tabs-mode:nil
  fill-column:99
  End:
*/
//vim:filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:fileencoding=utf-8:textwidth=99:
