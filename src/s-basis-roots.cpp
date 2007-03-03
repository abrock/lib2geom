#include <cmath>
#include <map>

#include "s-basis.h"
#include "sbasis-to-bezier.h"
#include "solver.h"

namespace Geom{

void bounds(SBasis const & s,
	    double &lo, double &hi, int order) {
    int imax=s.size()-1;
    lo=0;
    hi=0;

    for(int i = imax; i >=order; i--) {
      double a=s[i][0];
      double b=s[i][1];
      double t;

      if (hi>0){t=((b-a)+hi)/2/hi;}
      if (hi<=0||t<0||t>1){
	hi=std::max(a,b);
      }else{
	hi=a*(1-t)+b*t+hi*t*(1-t);	  
      }

      if (lo<0){t=((b-a)+lo)/2/lo;}
      if (lo>=0||t<0||t>1){
	lo=std::min(a,b);
      }else{
	lo=a*(1-t)+b*t+lo*t*(1-t);	  
      }
    }
    if (order>0){
        lo*=pow(.25,order);
        hi*=pow(.25,order);
    }
}

void local_bounds(SBasis const & s,
                  double t0,double t1,
                  double &lo, double &hi,
                  int order) {
    int imax=s.size()-1;
    lo=0;
    hi=0;
    
    for(int i = imax; i >=order; i--) {
        double a=s[i][0];
        double b=s[i][1];
        double t;
        if (hi>0){t=((b-a)+hi)/2/hi;}
        if (hi<=0||t<t0||t>t1){
            hi=std::max(a*(1-t0)+b*t0+hi*t0*(1-t0),a*(1-t1)+b*t1+hi*t1*(1-t1));
        }else{
            hi=a*(1-t)+b*t+hi*t*(1-t);	  
        }
        if (lo<0){t=((b-a)+lo)/2/lo;}
        if (lo>=0||t<t0||t>t1){
            lo=std::min(a*(1-t0)+b*t0+lo*t0*(1-t0),a*(1-t1)+b*t1+lo*t1*(1-t1));
        }else{
            lo=a*(1-t)+b*t+lo*t*(1-t);	  
        }
    }
    if (order>0){
        lo*=pow(.25,order);
        hi*=pow(.25,order);
    }
}

//-- multi_roots ------------------------------------
// goal: solve f(t)=c for several c at once.
// algo: compute f at both ends of the given segment.
//       compute bounds for df on the segment.
//       deduce first and last times when a level can be reached. 
//       if not empty, do the same with this smaller (and devided) segment. 
//TODO: Make sure the code is "rounding-errors proof"!


static int upper_level(vector<double> const &levels,double x,double tol=0.){
    for (int i=0;i<levels.size();i++){
        if (x<=levels[i]+tol) return(i);
    }
    return(levels.size());
}

static void multi_roots_internal(SBasis const &f,
				 SBasis const &df,
				 std::vector<double> const &levels,
				 std::map<double,unsigned> &roots,
				 double tol,
				 double a,
				 double fa,
				 double b,
				 double fb){
    
    if (f.size()==0){
        int idx;
        idx=upper_level(levels,0,tol);
        if (idx<levels.size()&&fabs(levels[idx])<=tol){
            roots[a]=idx;
            roots[b]=idx;
        }
        return;
    }
////usefull? 
//     if (f.size()==1){
//         int idxa=upper_level(levels,fa);
//         int idxb=upper_level(levels,fb);
//         if (fa==fb){
//             if (fa==levels[idxa]){
//                 roots[a]=idxa;
//                 roots[b]=idxa;
//             }
//             return;
//         }
//         int idx_min=std::min(idxa,idxb);
//         int idx_max=std::max(idxa,idxb);
//         if (idx_max==levels.size()) idx_max-=1;
//         for(int i=idx_min;i<=idx_max; i++){
//             double t=a+(b-a)*(levels[i]-fa)/(fb-fa);
//             if(a<t&&t<b) roots[t]=i;
//         }
//         return;
//     }
    if ((b-a)<tol){
        //TODO: use different tol for t and f ?
        int idx=std::min(upper_level(levels,fa,tol),upper_level(levels,fb,tol));
        if (idx==levels.size()) idx-=1;
        double c=levels[idx];
        if((fa-c)*(fb-c)<=0||fabs(fa-c)<tol||fabs(fb-c)<tol){
            roots[(a+b)/2]=idx;
        }
        return;
    }
    
    int idxa=upper_level(levels,fa,tol);
    int idxb=upper_level(levels,fb,tol);
    double m,M;
    local_bounds(df,a,b,m,M);

    //first times when a level (higher or lower) can be riched from a or b.
    double ta_hi,tb_hi,ta_lo,tb_lo;
    ta_hi=ta_lo=b+1;//default values => no root there.
    tb_hi=tb_lo=a-1;//default values => no root there.

    if (fabs(fa-levels[idxa])<tol){
        ta_hi=ta_lo=a;
    }else{
        if (M>0 && idxa<levels.size())
            ta_hi=a+(levels[idxa  ]-fa)/M;
        if (m<0 && idxa>0)
            ta_lo=a+(levels[idxa-1]-fa)/m;
    }
    if (fabs(fb-levels[idxb])<tol){
        tb_hi=tb_lo=b;
    }else{
        if (m<0 && idxb<levels.size())
            tb_hi=b+(levels[idxb  ]-fb)/m;
        if (M>0 && idxb>0)
            tb_lo=b+(levels[idxb-1]-fb)/M;
    }
    
    double t0,t1;
    t0=std::min(ta_hi,ta_lo);    
    t1=std::max(tb_hi,tb_lo);
    if (t0>t1) return;
    
    double t,ft;
    if (t1-t0<tol){
        multi_roots_internal(f,df,levels,roots,tol,t0,f(t0),t1,f(t1));
    }else{
        t=(t0+t1)/2;
        ft=f(t);
        multi_roots_internal(f,df,levels,roots,tol,t0,f(t0),t ,ft   );
        multi_roots_internal(f,df,levels,roots,tol,t ,ft   ,t1,f(t1));
    }
}

std::map<double,unsigned> multi_roots(SBasis const &f,
                                      std::vector<double> const &levels,
                                      double tol,
                                      double a,
                                      double b){
    std::map<double,unsigned> roots;
    SBasis df=derivative(f);
    multi_roots_internal(f,df,levels,roots,tol,a,f(a),b,f(b));  
    return(roots);
}
//-------------------------------------

#if 0
double Laguerre_internal(SBasis const & p,
                         double x0,
                         double tol,
                         bool & quad_root) {
    double a = 2*tol;
    double xk = x0;
    double n = p.size();
    quad_root = false;
    while(a > tol) {
        //std::cout << "xk = " << xk << std::endl;
        BezOrd b = p.back();
        BezOrd d(0), f(0);
        double err = fabs(b);
        double abx = fabs(xk);
        for(int j = p.size()-2; j >= 0; j--) {
            f = xk*f + d;
            d = xk*d + b;
            b = xk*b + p[j];
            err = fabs(b) + abx*err;
        }
        
        err *= 1e-7; // magic epsilon for convergence, should be computed from tol
        
        double px = b;
        if(fabs(b) < err)
            return xk;
        //if(std::norm(px) < tol*tol)
        //    return xk;
        double G = d / px;
        double H = G*G - f / px;
        
        //std::cout << "G = " << G << "H = " << H;
        double radicand = (n - 1)*(n*H-G*G);
        //assert(radicand.real() > 0);
        if(radicand < 0)
            quad_root = true;
        //std::cout << "radicand = " << radicand << std::endl;
        if(G.real() < 0) // here we try to maximise the denominator avoiding cancellation
            a = - std::sqrt(radicand);
        else
            a = std::sqrt(radicand);
        //std::cout << "a = " << a << std::endl;
        a = n / (a + G);
        //std::cout << "a = " << a << std::endl;
        xk -= a;
    }
    //std::cout << "xk = " << xk << std::endl;
    return xk;
}
#endif

void subdiv_sbasis(SBasis const & s,
                   std::vector<double> & roots, 
                   double left, double right) {
    double lo, hi;
    bounds(s, lo, hi);
    if(lo > 0 || hi < 0)
        return; // no roots here
    if(s.tail_error(1) < 1e-7) {
        double t = s[0][0] / (s[0][0] - s[0][1]);
        roots.push_back(left*(1-t) + t*right);
        return;
    }
    double middle = (left + right)/2;
    subdiv_sbasis(compose(s, BezOrd(0, 0.5)), roots, left, middle);
    subdiv_sbasis(compose(s, BezOrd(0.5, 1.)), roots, middle, right);
}

// It is faster to use the bernstein root finder for small degree polynomials (<100?.

std::vector<double> roots(SBasis const & s) {
    if(s.size() == 0) return std::vector<double>();
    std::vector<double> b = sbasis_to_bezier(s), r;
    
    find_bernstein_roots(&b[0], b.size()-1, r, 0, 0., 1.);
    return r;
}

};

/*
  Local Variables:
  mode:c++
  c-file-style:"stroustrup"
  c-file-offsets:((innamespace . 0)(inline-open . 0)(case-label . +))
  indent-tabs-mode:nil
  fill-column:99
  End:
*/
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:encoding=utf-8:textwidth=99 :
