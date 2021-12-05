//
//  ccars.cpp
//
// Code adapted from
// https://github.com/hunzikp/arscpp
// Adaptive Rejection Sampler
// Philipp Hunziker June 3, 2018



#include "ccars.hpp"

#include "data_utils.hpp"
#include "gnu_plotter.hpp"


#include "Eigen/Dense"
#include "Eigen/Core"


#include <algorithm>
using namespace std;

double CCARS::logSumExp(double u, double v) {
    double rv  = max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, v)));
    return rv;
}
// returns the log of sum of exps
// computes ls = log(sum(exp(x),dim))
// but in a way that tries to avoid underflow/overflow
//
// basic idea: shift before exp and reshift back
// alpha  = max(values)
// log(sum(exp(x))) = alpha + log(sum(exp(x-alpha)));
double CCARS::logSumExp(const Eigen::VectorXd& values ) {
    double rv = 0;
    double maxValue = values.maxCoeff();
    //maxValueVector.Constant(maxValue);
    //get the indexes >-Inf in values and maxValueVector:
    //libigl::slice(A,indices,B)?
    rv  = maxValue + log((values.array() - maxValue).exp().sum()) ;
    return rv;
}
double CCARS::logDiffExp(double u, double v) {
    assert(u>v);
    double rv = max(u, v) + log(exp(u - max(u, v)) - exp(v - max(u, v)));
    return rv;
}


UpperHull::UpperHull(const Eigen::VectorXd& x_internals, const Eigen::VectorXd& hx_internals, const Eigen::VectorXd& hpx_internals, double _xlb, double _xrb, double hx_xlb, double hx_xrb,  double hpx_xlb, double hpx_xrb, int max_points):
hx_xlb(hx_xlb), hx_xrb(hx_xrb), hpx_xlb(hpx_xlb), hpx_xrb(hpx_xrb),
x(x_internals), hx(hx_internals), tangents_slopes(hpx_internals)  {
    // x assumed sorted!
    
    // No of segments
    n_segments = x_internals.size();
    
    // Set bounds
    xlb = _xlb;
    xrb = _xrb;
    
    x_intercepts.resize(x_internals.size()-1);
    hx_intercepts.resize(x_internals.size()-1);
    
    bool is_concave =  hpx_xlb<0;
    for (int i = 0; i < n_segments; ++i) {
        
        HullSegment segment;
        if (i==0){
            segment.z_left = xlb;
        }
        else if(i==(n_segments-1)){
            segment.z_right = xlb;
        }
        segment.x = x_internals(i);
        segment.hx = hx_internals(i);
        segment.hpx = hpx_internals(i);
        
        if (is_concave)
            assert(segment.hpx <0 );
        else
            assert(segment.hpx >=0 );
        segments.push_back(segment);
    }
    // Update z and scum
    update_intercepts_from_pos(0);
   // update_cumulative_from_pos(0);
};

void UpperHull::update_intercepts_from_pos(int pos) {
    
    for (int i = pos; i < n_segments; ++i) {
        
        HullSegment segment = segments.at(i); // We make a copy and reassign to segments below; inefficient, but safer
        
        // Assign left z
        if (i == 0) {
            segment.z_left = xlb;
        } else {
            segment.z_left = segments.at(i-1).z_right;
            
        }
        // Assign right z
        if (i < n_segments - 1) {
            HullSegment next_segment;
            next_segment = segments.at(i+1);
            double num = segment.hx - next_segment.hx + next_segment.hpx*(next_segment.x - segment.x);
            double denom = next_segment.hpx - segment.hpx;
            if (num*denom<=0){
                
               // std::cout<< "num*denom"<< num*denom<< std::endl;
                
            }
            if (num*denom<=0){
                //in this case the segment.z_right will be very close to segment.x
                
                segment.z_right = segment.x+0.5*(next_segment.x-segment.x);
                
            }
            else if (abs(denom) <1e-6 ){
                //the 2 tangents have the same slope
                segment.z_right = segment.x+0.5*(next_segment.x-segment.x);
                
            }
            else{
                segment.z_right = segment.x + (num/denom);
            }
            
            assert(segment.z_right >xlb);
            assert(segment.z_right <=xrb);
        } else {
            segment.z_right = xrb;
        }
        // Assign hull values
        segment.hu_left = segment.hx - segment.hpx*(segment.x - segment.z_left);
        segment.hu_right = segment.hx + segment.hpx*(segment.z_right - segment.x);
        
        // Reassign segment copy
        segments.at(i) = segment;
    }
};
void UpperHull::get_intercepts(Eigen::VectorXd &x_intercept, Eigen::VectorXd &hx_intercept){
    
    int k = 0;
    assert(x_intercept.size() ==(n_segments+1));
    assert(hx_intercept.size() ==(n_segments+1));
    for (int i = 0; i < n_segments; ++i) {
        if (i == 0) {
            x_intercept(k)=  segments[i].z_left;
            hx_intercept(k)=  segments[i].hu_left;
            k++;
        }
        
        x_intercept(k)=  segments[i].z_right;
        hx_intercept(k)=  segments[i].hu_right;
        k++;
        
    }
    
}
void UpperHull::update_cumulative_from_pos(int pos) {
    
    // Compute normalizer
    vector<double> log_cu_vec(n_segments);
    log_cu = 0;
    for (int i = pos; i < n_segments; ++i) {
        HullSegment segment;
        segment = segments.at(i);
        // double cu_i = (1/segment.hpx)*(exp(segment.hu_right) - exp(segment.hu_left));
        double log_cu_i;
        if (segment.hpx > 0) {
            log_cu_i = -log(segment.hpx) + CCARS::logDiffExp(segment.hu_right, segment.hu_left);
        } else {
            log_cu_i = CCARS::logDiffExp(segment.hu_left - log(-segment.hpx), segment.hu_right - log(-segment.hpx));
        }
        
        if (i == 0) {
            log_cu = log_cu_i;
        } else {
            log_cu = CCARS::logSumExp(log_cu, log_cu_i);
        }
        log_cu_vec.at(i) = log_cu_i;
    }
    
    // Compute and assign scum
    vector<double> scum_vec(n_segments);
    for (int i = 0; i < n_segments; ++i) {
        if (i == 0) {
            scum_vec.at(0) = exp(log_cu_vec.at(0) - log_cu);
            segments.at(0).scum_left = 0;
            segments.at(0).scum_right = scum_vec.at(0);
        } else {
            scum_vec.at(i) = scum_vec.at(i-1) + exp(log_cu_vec.at(i) - log_cu);
            segments.at(i).scum_left = segments.at(i-1).scum_right;
            segments.at(i).scum_right = scum_vec.at(i);
        }
    }
};

double UpperHull::sample(gsl_rng *rng) {
    
    double x_sample;
    double u = Random::randomUniformFromGsl2(rng);
    double min_scum = segments.at(0).scum_right;
    if (u < min_scum) {
        // Sample is on the left-most segment
        HullSegment s = segments.at(0);
        
        // Draw a sample from the upper hull
        
        // x_sample = s.z_right + (1/s.hpx)*log(1 + (s.hpx*exp(log_cu)*(u-s.scum_right))/exp(s.hu_right));
        // ^ This often leads to over/underflow
        
        if (s.hpx*(u-s.scum_right) > 0) { // LH term might be negative, so we can't compute the log naively
            x_sample = s.z_right + (1/s.hpx)*CCARS::logSumExp(0, log(s.hpx*(u-s.scum_right)) + log_cu - s.hu_right);
        } else {
            x_sample = s.z_right + (1/s.hpx)*CCARS::logDiffExp(0, log(-1*s.hpx*(u - s.scum_right)) + log_cu - s.hu_right);
        }
        
    } else {
        // Determine which segment the sample is on
        int segment_id;
        for (int i = 1; i < n_segments; ++i) {
            if (u > segments.at(i).scum_left && u < segments.at(i).scum_right) {
                segment_id = i;
            }
        }
        HullSegment s = segments.at(segment_id);
        
        // Draw a sample from the upper hull
        
        // x_sample = s.z_left + (1/s.hpx)*log(1 + (s.hpx*exp(log_cu)*(u-s.scum_left))/exp(s.hu_left));
        // ^ This is prone to over/underflow
        
        if (s.hpx > 0) { // s.hpx might be negative, so we can't compute the log naively!
            x_sample = s.z_left + (1/s.hpx)*CCARS::logSumExp(0, log(s.hpx) + log_cu + log(u - s.scum_left) - s.hu_left);
        } else {
            x_sample = s.z_left + (1/s.hpx)*CCARS::logDiffExp(0, log(-s.hpx) + log_cu + log(u - s.scum_left) - s.hu_left);
        }
    }
    
    return x_sample;
};

void UpperHull::add_segment(double new_x, double new_hx, double new_hpx) {
    
    // Determine segment position
    int pos = 0;
    while(pos < n_segments) {
        if (new_x < segments.at(pos).x) {
            break;
        }
        pos++;
    }
    
    add_segment(new_x,  new_hx,  new_hpx,  pos);
    
};
void UpperHull::add_segment(double new_x, double new_hx, double new_hpx, int pos) {
    
    HullSegment segment;
    segment.x = new_x;
    segment.hx = new_hx;
    segment.hpx = new_hpx;
    
    int internal_position = pos;
    //    if (pos ==1){
    //         internal_position = 0;
    //
    //    }
    //     else
    //         internal_position = pos;
    
    assert(internal_position <= x.size());
    assert(internal_position <= hx.size());
    assert(internal_position <= tangents_slopes.size());
    
    Eigen::VectorXd new_xs(x.size()+1);
    Eigen::VectorXd new_hxs(hx.size()+1);
    Eigen::VectorXd new_tangents_slopes(tangents_slopes.size()+1);
    
    new_xs.head(internal_position)= x.head(internal_position);
    new_hxs.head(internal_position)= hx.head(internal_position);
    new_tangents_slopes.head(internal_position)= tangents_slopes.head(internal_position);
    
    new_xs(internal_position)= new_x;
    new_hxs(internal_position)= new_hx;
    new_tangents_slopes(internal_position) =new_hpx;
    
    new_xs.tail(x.size()-internal_position)= x.tail(x.size()-internal_position);
    new_hxs.tail(x.size()-internal_position)= hx.tail(x.size()-internal_position);
    new_tangents_slopes.tail(x.size()-internal_position)= tangents_slopes.tail(x.size()-internal_position);
    
    
    x= new_xs;
    hx = new_hxs;
    
    tangents_slopes =new_tangents_slopes;
    assert(x(internal_position)==new_x);
    assert(hx(internal_position)==new_hx);
    assert(new_tangents_slopes(internal_position)==new_hpx);
    
    // Insert segment
    segments.insert(segments.begin()+internal_position, segment);
    n_segments = segments.size();
    
    // Update z and scum
    int idx_update_intercepts = (internal_position>0)?internal_position-1: 0;
    update_intercepts_from_pos(idx_update_intercepts);
    //update_cumulative_from_pos(0);
    
};

double UpperHull::get_hu(double x) {
    
    // Determine which segment x lies on
    int segment_id = 0;
    for (int i = 0; i < n_segments; ++i) {
        if (x > segments.at(i).z_left && x <= segments.at(i).z_right) {
            segment_id = i;
            break;
        }
    }
    HullSegment s = segments.at(segment_id);
    
    // Get hu
    double hu;
    if (segment_id == 0) {
        hu = s.hu_right - (s.z_right - x)*s.hpx;
    } else {
        hu = s.hu_left + (x - s.z_left)*s.hpx;
    }
    
    return hu;
};
int UpperHull::get_point_position(double new_x){
    
    // Determine point position
    int iter = 0;
    while (iter < x.size()) {
        if (new_x < x(iter)) {
            break;
        }
        iter++;
    }
    return iter;
}
///_____________________________________________________________________________________________
LowerHull::LowerHull(const Eigen::VectorXd&  x,const  Eigen::VectorXd&  hx, int max_points):x(x),hx(hx) {
    // x assumed sorted!
    
    n_points = x.size();
    assert(n_points>=2);
    
    slopes.resize(n_points);
    for(size_t i= 0; i <(n_points-1) ; i++){
        double x_left = x(i);
        double x_right = x(i+1);
        double hx_left = hx(i);
        double hx_right = hx(i+1);
        double d = (hx_right-hx_left)/(x_right-x_left);
        
        slopes(i)= d;
        
    }
    slopes(n_points-1)= slopes(n_points-2);//for having the same
    //length in vectors x, hx, hpx
};

int LowerHull::get_point_position(double new_x){
    
    // Determine point position
    int iter = 0;
    while (iter < n_points) {
        if (new_x < x(iter)) {
            break;
        }
        iter++;
    }
    return iter;
}
void LowerHull::add_segment(double new_x, double new_hx, double new_hpx) {
    
    // Determine point position
    int pos =get_point_position( new_x);
    
    add_segment(new_x,  new_hx,  new_hpx, pos);
}

void LowerHull::add_segment(double new_x, double new_hx, double new_hpx, int position) {
    
    
    Eigen::VectorXd new_xs(x.size()+1);
    Eigen::VectorXd new_hxs(hx.size()+1);
    Eigen::VectorXd new_slopes(slopes.size()+1);
    
    new_xs.head(position)= x.head(position);
    new_hxs.head(position)= hx.head(position);
    new_slopes.head(position)= slopes.head(position);
    
    new_xs(position)= new_x;
    new_hxs(position)= new_hx;
    
    new_xs.tail(x.size()-position)= x.tail(x.size()-position);
    new_hxs.tail(x.size()-position)= hx.tail(x.size()-position);
    new_slopes.tail(x.size()-position)= slopes.tail(x.size()-position);
    
    x= new_xs;
    hx =new_hxs;
    
    assert(x(position)==new_x);
    assert(hx(position)==new_hx);
    
    double hx_right = hx[position +1];
    double hx_left = hx[position -1];
    double x_right = x[position +1];
    double x_left = x[position -1];
    double new_hpx_left = (new_hx-hx_left)/(new_x-x_left);
    double new_hpx_right = (hx_right-new_hx)/(x_right-new_x);
    
    
    new_slopes(position-1)= new_hpx_left;
    new_slopes(position)= new_hpx_right;
    
    if (position==n_points-1)
        new_slopes(position+1)= new_slopes(position);
    //new_hpxs[iter-1]= new_hpx_left;
    //hpx.insert(hpx.begin()+iter, new_hpx_right);
    slopes = new_slopes;
    
    n_points = x.size();
}
double LowerHull::get_hl(double _x) {
    
    // Determine point position
    int iter = 0;
    while (iter < n_points) {
        if (_x < x(iter)) {
            break;
        }
        iter++;
    }
    
    double rv;
    if (iter == 0) {
        rv = -numeric_limits<double>::infinity();
    } else if (iter == n_points) {
        rv = -numeric_limits<double>::infinity();
    } else {
        double x_left = x(iter-1);
        double x_right = x(iter);
        double hx_left = hx(iter-1);
        double hx_right = hx(iter);
        double d = (hx_right-hx_left)/(x_right-x_left);
        rv = hx_left + (_x-x_left)*d;
    }
    
    return rv;
}
double LowerHull::get_slope(int segment_idx) {
    
    
    return slopes(segment_idx);
}
// the x of intercepts points are inside the segments
Eigen::VectorXd LowerHull::get_h_intercepts( const Eigen::VectorXd& x_intercepts) {
    
    Eigen::VectorXd h_intercepts(n_points-3);
    int num_intercepts = x_intercepts.size();
    
    h_intercepts = slopes(Eigen::seqN(1, num_intercepts)).array()*(x_intercepts-x(Eigen::seqN(1,num_intercepts))).array() + hx(Eigen::seqN(1, num_intercepts)).array();
    
    return h_intercepts;
}
///_____________________________________________________________________________________________
CCUpperHull::CCUpperHull(UpperHull &upper_hull_concave, LowerHull &upper_hull_convex, double xlb, double xrb, int max_points ) :
upper_hull_concave(&upper_hull_concave),
upper_hull_convex(&upper_hull_convex),
xlb(xlb),
xrb(xrb)
{
    //from upper_hull_concave  we get the even points of CCUpperHull
    //from upper_hull_convex we get the odd points of CCUpperHull
    x.resize(2*upper_hull_concave.n_segments+1);
    hx.resize(2*upper_hull_concave.n_segments+1);
    hpx.resize(2*upper_hull_concave.n_segments+1);
    
    //upper_hull_concave.n_segments is the number of internal points
    
    init(upper_hull_concave, upper_hull_convex,
         x, hx , hpx);
    
    
    update_cumulative();
    
    
}
void CCUpperHull::init(UpperHull &upper_hull_concave, LowerHull &upper_hull_convex,
                       Eigen::VectorXd& x_temp,Eigen::VectorXd& hx_temp,Eigen::VectorXd& hpx_temp){
    
    Eigen::VectorXd upper_hull_concave_x_intercepts(upper_hull_concave.n_segments+1);
    Eigen::VectorXd upper_hull_concave_hx_intercepts(upper_hull_concave.n_segments+1);
    
    upper_hull_concave.get_intercepts(upper_hull_concave_x_intercepts, upper_hull_concave_hx_intercepts);
    
    x_temp( Eigen::seqN(0,upper_hull_concave.n_segments+1,2))  =  upper_hull_concave_x_intercepts;
    
    x_temp(Eigen::seqN(1,upper_hull_convex.n_points-2,2))  =  upper_hull_convex.x(Eigen::seqN(1, Eigen::last-1));
    
    
    hx_temp(0) = upper_hull_concave.segments[0].hu_left + upper_hull_convex.hx(0);
    hx_temp(hx_temp.size()-1) = upper_hull_concave.segments[upper_hull_concave.n_segments-1].hu_right + upper_hull_convex.hx(upper_hull_convex.hx.size()-1);
    
    hpx_temp(0) =upper_hull_concave.tangents_slopes(0) + upper_hull_convex.slopes(0);
    hpx_temp(hpx_temp.size()-1) = upper_hull_concave.tangents_slopes(upper_hull_concave.n_segments-1)  + upper_hull_convex.slopes(upper_hull_convex.slopes.size()-1);
    
    
    
    //now  tangent points in concave upper bound(= internal points in convex upper bound)
    hx_temp(Eigen::seqN(1,upper_hull_convex.n_points-2,2))=upper_hull_concave.hx +
    upper_hull_convex.hx(Eigen::seqN(1, upper_hull_convex.n_points-2));
    
    hpx_temp(Eigen::seqN(1,upper_hull_convex.n_points-2,2))=upper_hull_concave.tangents_slopes +
    upper_hull_convex.slopes(Eigen::seqN(1, upper_hull_convex.n_points-2));
    
    //now intercept points of tangents in concave upper bound
    hx_temp(Eigen::seqN(2,upper_hull_concave.n_segments-1,2))= upper_hull_concave_hx_intercepts(Eigen::seqN(1, upper_hull_concave.n_segments-1)) + upper_hull_convex.get_h_intercepts(upper_hull_concave_x_intercepts(Eigen::seqN(1, upper_hull_concave.n_segments-1))) ;
    
    hpx_temp(Eigen::seqN(2,upper_hull_concave.n_segments-1,2))= upper_hull_concave.tangents_slopes(Eigen::seqN(1, upper_hull_concave.n_segments-1)) + upper_hull_convex.slopes(Eigen::seqN(1, upper_hull_concave.n_segments-1));
    
    
}
void CCUpperHull::update_cumulative(){
    int num_points = x.size();
    cumulative.resize(num_points-1);
    
    size_t current_num_tangent_points = upper_hull_concave->x.size();
    assert((num_points-1)==2*current_num_tangent_points);//num_points include tangent points, xlb, xrb and intercept points
    cumulative.setZero();
    for(size_t i = 0; i < num_points; ++i){
        if (hpx(i)==0){
            assert(x(i+1) >x(i));
            cumulative(i) = hx(i) + log(x(i+1)-x(i));
        }
        else if(hpx(i)>0){//positive slope
            if (i<num_points-1){
                cumulative(i) = hx(i+1) - log(hpx(i))+log(1-exp(hpx(i)*(x(i)-x(i+1))));
            }
            
            
        }
        else{//negative slope
            if (i<num_points-1)
                cumulative(i) = hx(i) - log(-hpx(i))+log(1-exp(hpx(i)*(x(i+1)-x(i))));
            else{
                
            }
        }
        
    }
    maxCumulative = cumulative.maxCoeff();
}
double CCUpperHull::sample(gsl_rng *rng) {
    
    return(cc_upper_hull.sample(rng));
    
}
double CCUpperHull::get_hu(double x) {
    
    return(cc_upper_hull.get_hu(x));
}

void CCUpperHull::add_segment(double x_new_point, double hx_concave, double hpx_concave,double hx_convex, double hpx_convex ) {
    
    int pos =  get_point_position( x_new_point);
    add_segment(x_new_point,  hx_concave,  hpx_concave, hx_convex,  hpx_convex,  pos );
    
}
void CCUpperHull::add_segment(double x_new_point, double hx_concave, double hpx_concave,double hx_convex, double hpx_convex, int pos_convex_hull ) {
    
    int idx_new_segment_concave ;//= pos/3 +1;
    
    idx_new_segment_concave = (pos_convex_hull+1)/2;// pos_convex_hull-1;
    //upper_hull_concave->add_segment(x_new_point, hx_concave,  hpx_concave, idx_new_segment_concave);
    upper_hull_concave->add_segment(x_new_point, hx_concave,  hpx_concave);
    
    int idx_new_segment_convex; //= pos/2 +1;
    idx_new_segment_convex = (pos_convex_hull<=1)? 1:(pos_convex_hull+1)/2 +1; //pos_convex_hull;
    //upper_hull_convex->add_segment(x_new_point, hx_convex,  hpx_convex, idx_new_segment_convex);
    upper_hull_convex->add_segment(x_new_point, hx_convex,  hpx_convex);
    
    Eigen::VectorXd upper_hull_concave_x_intercepts(upper_hull_concave->n_segments+1);
    Eigen::VectorXd upper_hull_concave_hx_intercepts(upper_hull_concave->n_segments+1);
    
    upper_hull_concave->get_intercepts(upper_hull_concave_x_intercepts, upper_hull_concave_hx_intercepts);
    
    Eigen::VectorXd new_x(2*upper_hull_concave->n_segments+1);
    Eigen::VectorXd new_hx(2*upper_hull_concave->n_segments+1);
    Eigen::VectorXd new_hpx(2*upper_hull_concave->n_segments+1);
    
    init(*upper_hull_concave, *upper_hull_convex,
         new_x, new_hx , new_hpx);
    
    
    
    x= new_x;
    hx = new_hx;
    hpx = new_hpx;
    update_cumulative();
    
}
int CCUpperHull::get_point_position(double new_x){
    
    // Determine point position
    int iter = 0;
    int n_points = x.size();
    while (iter < n_points) {
        if (new_x < x(iter)) {
            break;
        }
        iter++;
    }
    return iter;
}
//___________________________________________________________________________

CCLowerHull::CCLowerHull(LowerHull &lower_hull_concave, UpperHull &lower_hull_convex, double xlb, double xrb, int max_points):
lower_hull_concave(&lower_hull_concave),
lower_hull_convex(&lower_hull_convex),
xlb(xlb),
xrb(xrb)
{
    
    x.resize(2*lower_hull_convex.n_segments+1);
    hx.resize(2*lower_hull_convex.n_segments+1);
    hpx.resize(2*lower_hull_convex.n_segments+1);
    
    //upper_hull_concave.n_segments is the number of internal points
    init(lower_hull_concave, lower_hull_convex,  x,hx, hpx);
    
    //    Eigen::VectorXd lower_hull_convex_x_intercepts(lower_hull_convex.n_segments+1);
    //    Eigen::VectorXd lower_hull_convex_hx_intercepts(lower_hull_convex.n_segments+1);
    //
    //    lower_hull_convex.get_intercepts(lower_hull_convex_x_intercepts, lower_hull_convex_hx_intercepts);
    //
    //    x( Eigen::seqN(0,lower_hull_convex_x_intercepts.size(),2))  =  lower_hull_convex_x_intercepts;
    //
    //    x(Eigen::seqN(1,lower_hull_concave.n_points-2,2))  =  lower_hull_concave.x(Eigen::seqN(1, Eigen::last-1));
    //
    //
    //    hx(0) = lower_hull_convex.segments[0].hu_left + lower_hull_concave.hx(0);
    //    hx(hx.size()-1) = lower_hull_convex.segments[lower_hull_convex.n_segments-1].hu_right + lower_hull_concave.hx(lower_hull_concave.hx.size()-1);
    //
    //    hpx(0) =lower_hull_convex.tangents_slopes(0) + lower_hull_concave.slopes(0);
    //    hpx(hpx.size()-1) = lower_hull_convex.tangents_slopes(lower_hull_convex.n_segments-1)  + lower_hull_concave.slopes(lower_hull_concave.slopes.size()-1);
    //
    //
    //    //now  tangent points in concave upper bound(= internal points in convex upper bound)
    //
    //    hx(Eigen::seqN(1,lower_hull_concave.n_points-2,2))=lower_hull_convex.hx +
    //    lower_hull_concave.hx(Eigen::seqN(1, lower_hull_concave.n_points-2));
    //
    //    hpx(Eigen::seqN(1,lower_hull_concave.n_points-2,2))=lower_hull_convex.tangents_slopes +
    //    lower_hull_concave.slopes(Eigen::seqN(1, lower_hull_concave.n_points-2));
    //
    //    //now intercept points of tangents in concave upper bound
    //    hx(Eigen::seqN(2,lower_hull_convex.n_segments-1,2))= lower_hull_convex_hx_intercepts(Eigen::seqN(1, lower_hull_convex.n_segments-1)) + lower_hull_concave.get_h_intercepts(lower_hull_convex_x_intercepts(Eigen::seqN(1, lower_hull_convex.n_segments-1))) ;
    //
    //    hpx(Eigen::seqN(2,lower_hull_convex.n_segments-1,2))= lower_hull_convex.tangents_slopes(Eigen::seqN(1, lower_hull_convex.n_segments-1)) + lower_hull_concave.slopes(Eigen::seqN(1, lower_hull_convex.n_segments-1));
    //
    update_cumulative();
    
    
}
void CCLowerHull::init(LowerHull &lower_hull_concave, UpperHull &lower_hull_convex,   Eigen::VectorXd& x_temp,Eigen::VectorXd& hx_temp,Eigen::VectorXd& hpx_temp){
    
    
    Eigen::VectorXd lower_hull_convex_x_intercepts(lower_hull_convex.n_segments+1);
    Eigen::VectorXd lower_hull_convex_hx_intercepts(lower_hull_convex.n_segments+1);
    
    lower_hull_convex.get_intercepts(lower_hull_convex_x_intercepts, lower_hull_convex_hx_intercepts);
    
    x_temp( Eigen::seqN(0,lower_hull_convex_x_intercepts.size(),2))  =  lower_hull_convex_x_intercepts;
    
    x_temp(Eigen::seqN(1,lower_hull_concave.n_points-2,2))  =  lower_hull_concave.x(Eigen::seqN(1, Eigen::last-1));
    
    
    hx_temp(0) = lower_hull_convex.segments[0].hu_left + lower_hull_concave.hx(0);
    hx_temp(hx_temp.size()-1) = lower_hull_convex.segments[lower_hull_convex.n_segments-1].hu_right + lower_hull_concave.hx(lower_hull_concave.hx.size()-1);
    
    hpx_temp(0) =lower_hull_convex.tangents_slopes(0) + lower_hull_concave.slopes(0);
    hpx_temp(hpx_temp.size()-1) = lower_hull_convex.tangents_slopes(lower_hull_convex.n_segments-1)  + lower_hull_concave.slopes(lower_hull_concave.slopes.size()-1);
    
    
    //now  tangent points in concave upper bound(= internal points in convex upper bound)
    
    hx_temp(Eigen::seqN(1,lower_hull_concave.n_points-2,2))=lower_hull_convex.hx +
    lower_hull_concave.hx(Eigen::seqN(1, lower_hull_concave.n_points-2));
    
    hpx_temp(Eigen::seqN(1,lower_hull_concave.n_points-2,2))=lower_hull_convex.tangents_slopes +
    lower_hull_concave.slopes(Eigen::seqN(1, lower_hull_concave.n_points-2));
    
    //now intercept points of tangents in concave upper bound
    hx_temp(Eigen::seqN(2,lower_hull_convex.n_segments-1,2))= lower_hull_convex_hx_intercepts(Eigen::seqN(1, lower_hull_convex.n_segments-1)) + lower_hull_concave.get_h_intercepts(lower_hull_convex_x_intercepts(Eigen::seqN(1, lower_hull_convex.n_segments-1))) ;
    
    hpx_temp(Eigen::seqN(2,lower_hull_convex.n_segments-1,2))= lower_hull_convex.tangents_slopes(Eigen::seqN(1, lower_hull_convex.n_segments-1)) + lower_hull_concave.slopes(Eigen::seqN(1, lower_hull_convex.n_segments-1));
    
    
    
}
void CCLowerHull::update_cumulative(){
    int num_points = x.size();
    cumulative.resize(num_points-1);
    size_t current_num_tangent_points = lower_hull_convex->x.size();
    assert((num_points-1)==2*current_num_tangent_points);
    cumulative.setZero();
    for(size_t i = 0; i < num_points; ++i){
        if (hpx(i)==0){
            assert(x(i+1)>x(i));
            cumulative(i) = hx(i) + log(x(i+1)-x(i));
        }
        else if(hpx(i)>0){//positive slope
            //if (i> 1)
            if (i<num_points-1)
                cumulative(i) = hx(i+1) - log(hpx(i))+log(1-exp(hpx(i)*(x(i)-x(i+1))));
            
            
        }
        else{//negative slope
            if (i<num_points-1)
                cumulative(i) = hx(i) - log(-hpx(i))+log(1-exp(hpx(i)*(x(i+1)-x(i))));
            //  else{
            
            //  }
            
        }
        
    }
    minCumulative = cumulative.minCoeff();
}
void CCLowerHull::add_segment(double x_new_point, double hx_concave, double hpx_concave,double hx_convex, double hpx_convex, int pos ) {
    
    
    
    //int idx_new_segment_convex =(pos+1)/2;//pos-1;// pos/3 +1;
    // int idx_new_segment_concave=(pos<=1)? 1:(pos+1)/2 +1; //pos;// pos/2 +1;
    
    // lower_hull_concave->add_segment(x_new_point, hx_concave,  hpx_concave, idx_new_segment_concave);
    lower_hull_concave->add_segment(x_new_point, hx_concave,  hpx_concave);
    
    // lower_hull_convex->add_segment(x_new_point, hx_convex,  hpx_convex, idx_new_segment_convex);
    lower_hull_convex->add_segment(x_new_point, hx_convex,  hpx_convex);
    
    Eigen::VectorXd upper_hull_concave_x_intercepts(lower_hull_convex->n_segments+1);
    Eigen::VectorXd upper_hull_concave_hx_intercepts(lower_hull_convex->n_segments+1);
    
    lower_hull_convex->get_intercepts(upper_hull_concave_x_intercepts, upper_hull_concave_hx_intercepts);
    
    Eigen::VectorXd new_x(2*lower_hull_convex->n_segments+1);
    Eigen::VectorXd new_hx(2*lower_hull_convex->n_segments+1);
    Eigen::VectorXd new_hpx(2*lower_hull_convex->n_segments+1);
    
    init(*lower_hull_concave, *lower_hull_convex,
         new_x, new_hx , new_hpx);
    
    
    x= new_x;
    hx = new_hx;
    hpx = new_hpx;
    update_cumulative();
    //cc_lower_hull.update_cumulative_from_pos(0);
    //cc_lower_hull.update_scum();
    
}
int CCLowerHull::get_point_position(double new_x){
    
    // Determine point position
    int iter = 0;
    int n_points = x.size();
    while (iter < n_points) {
        if (new_x < x(iter)) {
            break;
        }
        iter++;
    }
    return iter;
}
//___________________________________________________________________________

CCARS::CCARS(CCLogDensity* const log_density,const Eigen::VectorXd& x, double xlb, double xrb, int max_points) :
log_density(log_density),
xlb(xlb),
xrb(xrb),
max_points(max_points)
{
    int num_points = x.size();
    // Check bounds
    if (xlb >= xrb) {
        std::cout<< "Upper bound is not larger than lower bound."<<std::endl;
    }
    
    // We need at least two starting points
    if (x.size() < 2) {
        std::cout<< "At least two starting points required."<<std::endl ;
    }
    // Order x
    // std::sort(x.data(), x.data()+x.size());
    // std::sort (x.begin(), x.end());
    
    assert(x[0]==xlb);
    assert(x[num_points-1]==xrb);
    
    
    // Check points in bounds
    if (x(0) < xlb || x(x.size()-1) > xrb) {
        std::cout<<"Starting point out of bound."<<std::endl;
    }
    
    Eigen::VectorXd internal_points= x.segment(1,num_points-2);
    
    // Evaluate funcs
    Eigen::VectorXd hx_concave_internals = log_density->h_concave(internal_points);
    Eigen::VectorXd  hpx_concave_internals= log_density->h_prime_concave(internal_points);
    
    Eigen::VectorXd hx_concave(hx_concave_internals.size()+2);
    hx_concave(0)=log_density->h_concave(xlb);
    hx_concave.segment(1,hx_concave.size()-2) = hx_concave_internals;
    hx_concave(hx_concave.size()-1) = log_density->h_concave(xrb);
    
    
    Eigen::VectorXd hpx_concave(hpx_concave_internals.size()+2);
    hpx_concave(0)=log_density->h_prime_concave(xlb);
    hpx_concave.segment(1,hpx_concave.size()-2) = hpx_concave_internals;
    hpx_concave(hpx_concave.size()-1) = log_density->h_prime_concave(xrb);
    
    
    Eigen::VectorXd hx_convex = log_density->h_convex(x);
    Eigen::VectorXd hpx_convex = log_density->h_prime_convex(x);
    
    
    //            // Try to make starting points valid (if they're invalid)
    //            int max_tries = 10;
    //            if (hpx_concave_internals(0) <= 0) {
    //                // Move left-most point left until valid
    //
    //                double hpx_left = hpx_concave_internals(0);
    //                double x_left = x(0);
    //                int tries = 0;
    //                while (hpx_left <= 0 && tries < max_tries) {
    //
    //                    if (isfinite(xlb)) {
    //                        x_left -= (x_left-xlb)/2; // Move half-way to the limit
    //                    } else {
    //                        x_left -= pow(2, tries); // Move left by power of 2
    //                    }
    //
    //                    hpx_left = log_density->h_prime_concave(x_left);
    //
    //                    tries++;
    //                }
    //
    //                if (tries < max_tries) {
    //                    hx_concave_internals(0) = hpx_left;
    //                    internal_points(0) = x_left;
    //                    hx_concave_internals(0) = log_density->h_concave(x_left);
    //                } else {
    //                    std::cout<< "Could not find valid lower starting point."<<std::endl;
    //                }
    //            }
    //
    //            int last_ind = hx_concave_internals.size() - 1;
    //            if (hpx_concave_internals(last_ind) >= 0) {
    //                // Move right-most point right until valid
    //
    //                double hpx_right = hpx_concave_internals(last_ind);
    //                double x_right = x(last_ind);
    //                int tries = 0;
    //                while (hpx_right >= 0 && tries < max_tries) {
    //
    //                    if (isfinite(xrb)) {
    //                        x_right += (xrb - x_right)/2; // Move half-way to the limit
    //                    } else {
    //                        x_right += pow(2, tries); // Move right by power of 2
    //                    }
    //
    //                    hpx_right = log_density->h_prime_concave(x_right);
    //
    //                    tries++;
    //                }
    //
    //                if (tries < max_tries) {
    //                    hx_concave_internals(last_ind) = hpx_right;
    //                    internal_points(last_ind) = x_right;
    //                    hx_concave_internals(last_ind) = log_density->h_concave(x_right);
    //                } else {
    //                    // std::cout << "x candidates: " << x << std::endl;
    //                    std::cout <<"Could not find valid upper starting point."<< std::endl;
    //                }
    //            }
    //

    // Create the hulls
    upper_hull_concave = UpperHull(internal_points, hx_concave_internals, hpx_concave_internals, xlb, xrb, hx_concave(0), hx_concave(hx_concave.size()-1), hpx_concave(0), hpx_concave(hpx_concave.size()-1), max_points);
    lower_hull_concave = LowerHull(x, hx_concave, max_points);
    
    lower_hull_convex = UpperHull(internal_points, hx_convex(Eigen::seqN(1, num_points-2)), hpx_convex(Eigen::seqN(1, num_points-2)), xlb, xrb, hx_convex(0), hx_convex(hx_convex.size()-1), hpx_convex(0), hpx_convex(hpx_convex.size()-1), max_points  );
    upper_hull_convex = LowerHull(x, hx_convex, max_points);
    
    upper_hull = CCUpperHull(upper_hull_concave, upper_hull_convex, xlb, xrb,  max_points);
    
    lower_hull = CCLowerHull(lower_hull_concave, lower_hull_convex, xlb, xrb,  max_points);
    
    //ratio of bound areas
     Eigen::ArrayXd diff1= upper_hull.cumulative.array()-upper_hull.maxCumulative;
    Eigen::ArrayXd diff2= lower_hull.cumulative.array()-upper_hull.maxCumulative;
    Eigen::ArrayXd  diff3= upper_hull.cumulative-lower_hull.cumulative;
    //using std::exp;
    
    Eigen::VectorXd  res1 =(upper_hull.cumulative.array()-upper_hull.maxCumulative).exp();
    Eigen::VectorXd  res2 =(lower_hull.cumulative.array()-upper_hull.maxCumulative).exp();
   
    //ratioAreas = (upper_hull.cumulative.array()-lower_hull.cumulative.array()).abs();
    ratioAreas = (upper_hull.cumulative.array()-upper_hull.maxCumulative).exp()-(lower_hull.cumulative.array()-upper_hull.maxCumulative).exp();
};

double CCARS::approximateLogIntegral(gsl_rng *rng,  double  max_number_tangent_points) {
    
    double result;
    
    size_t current_num_tangent_points = upper_hull_concave.x.size();
    
    Eigen::VectorXd boundedArea(current_num_tangent_points+1);
    Eigen::VectorXd absDiffAreas;
    double tmp, u_bound, l_bound;
    size_t i, j, pos;
    double new_x;
    bool useDiffAreas = true;
    double logCumAreaInterval = 0.0;
    double hconcave, hconvex, h;
    size_t trials = upper_hull_concave.x.size();
    //std::cout<< ratioAreas << std::endl;
    while(trials < max_number_tangent_points && exp(logSumExp(lower_hull.cumulative)-logSumExp(upper_hull.cumulative))< 0.999 )//exp(logSumExp(lower_hull.cumulative)-logSumExp(upper_hull.cumulative)) should be always be < 1
    {
        //choose a segment
        current_num_tangent_points = upper_hull_concave.x.size();
        boundedArea.resize(current_num_tangent_points+1);
        absDiffAreas = (upper_hull.cumulative.array()-lower_hull.cumulative.array()).abs();
        
        if (!useDiffAreas){
        
        boundedArea(0) =ratioAreas(0);
        boundedArea(boundedArea.size()-1) = ratioAreas(ratioAreas.size()-1);
        boundedArea(Eigen::seqN(1,boundedArea.size()-2 ))=  ratioAreas(Eigen::seqN(1, current_num_tangent_points-1, 2 ))+ratioAreas(Eigen::seqN(2,current_num_tangent_points-1, 2 ));
        }
        else{
        
        
        boundedArea(0) =absDiffAreas(0);
        boundedArea(boundedArea.size()-1) = absDiffAreas(absDiffAreas.size()-1);
        boundedArea(Eigen::seqN(1,boundedArea.size()-2 ))=  absDiffAreas(Eigen::seqN(1, current_num_tangent_points-1, 2 ))+absDiffAreas(Eigen::seqN(2,current_num_tangent_points-1, 2 ));
            
        }
        
//        std::cout << "boundedArea \n" << boundedArea << std::endl;
//        std::cout << "ratioAreas \n" << ratioAreas << std::endl;
//        std::cout << "upper x \n" << upper_hull.x << std::endl;
        tmp = boundedArea.maxCoeff();
        
    
        if (boundedArea(0)==tmp){
            //std::cout << " first segment " << std::endl;
            pos= 0;
            new_x = sampleNewX(rng, xlb, upper_hull.x(1), upper_hull.hx(0),upper_hull.hx(1), upper_hull.hpx(0), logCumAreaInterval );
        }
        else if (boundedArea(boundedArea.size()-1)==tmp){//the max is in the last position
            //std::cout << " last segment " << std::endl;
            i= 2*current_num_tangent_points-1;
            pos = boundedArea.size()-1;
            new_x = sampleNewX(rng, upper_hull.x(i), xrb, upper_hull.hx(i),upper_hull.hx(i+1), upper_hull.hpx(i), logCumAreaInterval);
            
            
        }
        else{
            // std::cout << " intermediate segment " << std::endl;
            j = 1;
            while (boundedArea(j)!=tmp)
                j++;
            
            pos = j;
            i= 2*j;
            j= (upper_hull.x(i)< lower_hull.x(i))?1:0;
            
            
            double term1 = upper_hull.hx(i) -lower_hull.hpx(i-j)*(upper_hull.x(i)-lower_hull.x(i-j))-lower_hull.hx(i-j);
            double term2 = upper_hull.hpx(i+j-1)*(lower_hull.x(i)-upper_hull.x(i+j-1))+upper_hull.hx(i+j-1)-lower_hull.hx(i);
            if (term1 >term2){
                
                new_x=upper_hull.x(i);
                //std::cout << "new x from in the upper hull " <<  std::endl;
            }
            else{
                new_x=lower_hull.x(i);
                //std::cout << "new x is from the lower hull " <<  std::endl;
            }
            
        }
        assert(i< upper_hull.hx.size());
        u_bound =upper_hull.hpx(i)*(new_x-upper_hull.x(i))+upper_hull.hx(i);
        
        if (lower_hull.x(i)< new_x){
            l_bound = lower_hull.hpx(i)*(new_x-lower_hull.x(i))+ lower_hull.hx(i) ;
            
        }else{
            l_bound = lower_hull.hpx(i-1)*(new_x-lower_hull.x(i-1))+ lower_hull.hx(i-1) ;
            
        }
        if( (upper_hull.x(i+1)- upper_hull.x(i)) > 1e-5){
            
            current_num_tangent_points++;
            hconcave =log_density->h_concave(new_x);
            hconvex =log_density->h_convex(new_x);;
            h = hconcave +hconvex;
           // std::cout << " New x "<< new_x <<" segment position: "<< pos << " out of " << boundedArea.size()-1 << std::endl;
            
            addPoint(new_x, hconcave, hconvex, i);
        }
        else{
            //dont include this point
            // std::cout << "point rejected " << std::endl;
             trials++;
        }
    }
    //log of upper bound of the integral area
    //result = upper_hull.cumulative.sum();
    result = logSumExp(upper_hull.cumulative);
//    std::cout<< " current_num_tangent_points "<< current_num_tangent_points << std::endl;
//    std::cout<< " log integral "<< result << std::endl;
//     std::cout<< " final diff "<< exp(logSumExp(lower_hull.cumulative)-logSumExp(upper_hull.cumulative)) << "\n" <<std::endl;
    return result;
};


double CCARS::sampleNewX(gsl_rng *rng,double x1, double x2, double y1, double y2, double m, double &logCumArea ){
    //samples new_x from  distribution with cdf
    //p(x) = (exp((y2-y1)/(x2-x1)*(x-x1))-1)/(exp(y2-y1)-1)
    double new_x;
    double m1;
    logCumArea= 0.0;
    double allIntervalArea ;
    double u = Random::randomUniformFromGsl2(rng);
    if (y2>y1){
        u= 1-log1p(u*expm1(y1-y2))/(y1-y2);
        new_x = u*x2+(1-u)*x1;
        m1 = (y2-y1)/(x2-x1);
        logCumArea = - log(m1) +y2 +log(1-exp(m1*(x1-x2)));;
        
       // allIntervalArea =  log((exp(m1*(x2-x1))-1)/(exp(y2-y1)-1));
       // assert(allIntervalArea == 0);

    }
    else if(y1>y2)
    {
        u= log1p(u*expm1(y2-y1))/(y2-y1);
        new_x = u*x2+(1-u)*x1;
        m1 = (y2-y1)/(x2-x1);
        logCumArea = y1 - log(-m1)+log(1-exp(m1*(x2-x1)));
        //allIntervalArea =  log((exp(m1*(x2-x1))-1)/(exp(y2-y1)-1));
        //assert(allIntervalArea == 0);
    }
    else{
        
        new_x = u*x2+(1-u)*x1;
        logCumArea = y1+log(new_x-x1);
        allIntervalArea =  0.0;
    }
    assert(!isnan(logCumArea));
    
//
//    if (hpx(i)==0){
//               assert(x(i+1) >x(i));
//               cumulative(i) = hx(i) + log(x(i+1)-x(i));
//           }
//           else if(hpx(i)>0){//positive slope
//               if (i<num_points-1){
//                   cumulative(i) = hx(i+1) - log(hpx(i))+log(1-exp(hpx(i)*(x(i)-x(i+1))));
//               }
//
//
//           }
//           else{//negative slope
//               if (i<num_points-1)
//                   cumulative(i) = hx(i) - log(-hpx(i))+log(1-exp(hpx(i)*(x(i+1)-x(i))));
//               else{
//
//               }
//           }
    return new_x;
}
double CCARS::sample(gsl_rng *rng, double& y, double& logComplementArea){
    double new_x;
    bool accepted = false;
    size_t i, j;
    size_t maxCumulative;
    logComplementArea = 0.0;
    double logcumArea = 0.0;
    double logIntegral;
    while (!accepted){
        
        maxCumulative = upper_hull.cumulative.maxCoeff();
        size_t num_segments = upper_hull.cumulative.size();
        Eigen::VectorXd result = (upper_hull.cumulative.array()-maxCumulative).exp();
        Eigen::VectorXd cumulativeSum(result.size());
        partial_sum(result.begin(), result.end(), cumulativeSum.begin(), plus<double>());
        size_t count = 0;
        double rnd = Random::randomUniformFromGsl2(rng);
        double term = rnd*cumulativeSum(Eigen::last);
        
        for(size_t i = 0; i<num_segments ; ++i){
            if (term >cumulativeSum(i))
                count++;
        }
        
      //  i = 1 + count;
        if (count >0)
            i = count -1;
        else
            i = 0;

        j= i;
//        if (j ==cumulativeSum.size()){
//            j = j -1;
//        }
        if (i >0 && term <cumulativeSum(j))
            std::cout << "oops " << std::endl;
        
        if (i >0 && term >cumulativeSum(j+1))
                 std::cout << "oops " << std::endl;
        
        assert(i  ==0 || term >= cumulativeSum(j));
        assert(i  ==0 || term <= cumulativeSum(j+1));
        
        logcumArea = maxCumulative +log(cumulativeSum(j));
        logIntegral = maxCumulative +log(cumulativeSum(Eigen::last));
        assert(!isnan(logIntegral));
        
        if (logcumArea <logIntegral)
           logComplementArea = CCARS::logDiffExp(logIntegral, logcumArea);
        else
            logComplementArea = 0.0;
        
        double logCumAreaInterval;
        i = i+1; //the first element in upper_hull.x is xlb and does not
        //correspond to a segment
        new_x = sampleNewX(rng, upper_hull.x(i), upper_hull.x(i+1), upper_hull.hx(i),upper_hull.hx(i+1), upper_hull.hpx(i), logCumAreaInterval );
        
//        double extraLogArea = 0.0;
//        if (upper_hull.hpx(i)==0){
//            assert(new_x >upper_hull.x(i));
//            extraLogArea+= upper_hull.hx(i) + log(new_x -upper_hull.x(i));
//        }
//        else if(upper_hull.hpx(i)>0){//positive slope
//
//            extraLogArea+= upper_hull.hx(i+1) - log(upper_hull.hpx(i))+log(1-exp(upper_hull.hpx(i)*(upper_hull.x(i)-new_x)));
//        }
//        else{//negative slope
//            extraLogArea+= upper_hull.hx(i) - log(-upper_hull.hpx(i))+log(1-exp(upper_hull.hpx(i)*(new_x-upper_hull.x(i))));
//
//        }
        y= logCumAreaInterval;
        assert(CCARS::logSumExp(logcumArea,logComplementArea) ==logIntegral );
        logcumArea  = CCARS::logSumExp(logcumArea, logCumAreaInterval) ;
      // assert((logComplementArea-logCumAreaInterval)>1e-6);
        if (logComplementArea>logCumAreaInterval)
           logComplementArea = logDiffExp(logComplementArea, logCumAreaInterval);
        else
            logComplementArea = 1.0;
//
        //assert(CCARS::logSumExp(logcumArea,logComplementArea) -logIntegral <1e-6 );
        //check for acceptance
        double yupperBoundValue = upper_hull.hpx(i)*(new_x-upper_hull.x(i)) + upper_hull.hx(i);
        double ylowerBoundValue;
        if (lower_hull.x(i) < new_x){
            ylowerBoundValue = lower_hull.hpx(i)*(new_x-lower_hull.x(i)) + lower_hull.hx(i);
        }
        else{
            ylowerBoundValue =  lower_hull.hpx(i-1)*(new_x-lower_hull.x(i-1)) + lower_hull.hx(i-1);
            
        }
        double u = Random::randomUniformFromGsl2(rng);
        if ( exp(ylowerBoundValue-yupperBoundValue)> u)//accept squeezing step
        {
            //double hconcave =log_density->h_concave(new_x);
            //double hconvex =log_density->h_convex(new_x);;
           // y = hconcave +hconvex;
            accepted = true;
            
        }
        else{
            double hconcave =log_density->h_concave(new_x);
            double hconvex =log_density->h_convex(new_x);;
            //y = hconcave +hconvex;
            if (exp(y- yupperBoundValue)>u){//accept sample at rejection step
                accepted = true;
            }
            else{//reject and update bounds
               // std::cout << " rejected " << std::endl;
                double dconcave =log_density->h_prime_concave(new_x);
                double dconvex =log_density->h_prime_convex(new_x);
                int curr_num_points = upper_hull.x.size();
                upper_hull.add_segment(new_x, hconcave,dconcave,  hconvex, dconvex, i);
                lower_hull.add_segment(new_x, hconcave,dconcave,  hconvex, dconvex, i);
                assert(upper_hull.x.size() == curr_num_points+2 );
                assert(lower_hull.x.size() == curr_num_points+2 );
               // ratioAreas = (upper_hull.cumulative.array()-lower_hull.cumulative.array()).abs();
                ratioAreas = exp(upper_hull.cumulative.array()-upper_hull.maxCumulative)-exp(lower_hull.cumulative.array()-upper_hull.maxCumulative);
            }
        }
    }
    
    return new_x;
}
void CCARS::addPoint(double new_x){
    
    int pos_new_x = lower_hull.lower_hull_concave->get_point_position(new_x);
    double  hconcave_new_x =log_density->h_concave(new_x);
    double hconvex_new_x =log_density->h_convex(new_x);;
    addPoint(new_x,  hconcave_new_x,  hconvex_new_x,  pos_new_x);
    
}
void CCARS::addPoint(double new_x, double hconcave_new_x, double hconvex_new_x, int pos_new_x){
    size_t current_num_tangent_points = upper_hull_concave.x.size();
    if (pos_new_x == 2*current_num_tangent_points-1 ){//2*current_num_tangent_points //(current_num_tangent_points)+1
        //new point is at the left of last tangent point
        double dconcave = log_density->h_prime_concave(new_x);
        double dconvex = log_density->h_prime_convex(new_x);
        
        if (xrb >new_x){
            current_num_tangent_points++;
            upper_hull.add_segment(new_x, hconcave_new_x,dconcave,  hconvex_new_x, dconvex, pos_new_x);
            
            lower_hull.add_segment(new_x, hconcave_new_x,dconcave,  hconvex_new_x, dconvex, pos_new_x);
            
          
        }
        
    }
    else if(pos_new_x==1){
        //new point is at the right of first tangent point
        double dconcave = log_density->h_prime_concave(new_x);
        double dconvex = log_density->h_prime_convex(new_x);
        if (xlb <new_x){
            current_num_tangent_points++;
            upper_hull.add_segment(new_x, hconcave_new_x,dconcave,  hconvex_new_x, dconvex, pos_new_x);
            
            lower_hull.add_segment(new_x, hconcave_new_x,dconcave,  hconvex_new_x, dconvex, pos_new_x);
            
        }
        
    }
    else{
        
        double dconcave = log_density->h_prime_concave(new_x);
        double dconvex = log_density->h_prime_convex(new_x);
        if (xlb <new_x){
            current_num_tangent_points++;
            upper_hull.add_segment(new_x, hconcave_new_x,dconcave,  hconvex_new_x, dconvex, pos_new_x);
            
            lower_hull.add_segment(new_x, hconcave_new_x,dconcave,  hconvex_new_x, dconvex, pos_new_x);
            
        }
        else{
            
            std::cout << "Error!. The new x cannot be greater than the left limit of interval" << std::endl;
        }
        
    }
    //ratioAreas = (upper_hull.cumulative.array()-lower_hull.cumulative.array()).abs();
    ratioAreas = (upper_hull.cumulative.array()-upper_hull.maxCumulative).exp()-(lower_hull.cumulative.array()-upper_hull.maxCumulative).exp();
    
}
void CCARS::plot( GNUPlotter &plotter, int num_grid_points, int iter,  int idxFirstID, int idxSecondId){
    
    
    Eigen::VectorXd xserieEigen = Eigen::VectorXd::LinSpaced(num_grid_points,xlb,xrb);
    std::vector<double> xseries(xserieEigen.data(),xserieEigen.data()+ xserieEigen.size() );
    Eigen::VectorXd yserieEigenConcave = log_density->h_concave(xserieEigen);
    Eigen::VectorXd yserieEigenConvex = log_density->h_convex(xserieEigen);
    Eigen::VectorXd yserieEigenFull= yserieEigenConcave+yserieEigenConvex;
    
    Eigen::VectorXd derivativeEigenConcave = log_density->h_prime_concave(xserieEigen);
    Eigen::VectorXd derivativeEigenConvex = log_density->h_prime_convex(xserieEigen);
    
    std::vector<double> yseriesConcave(yserieEigenConcave.data(),yserieEigenConcave.data()+ yserieEigenConcave.size());
    std::vector<double> yseriesConvex(yserieEigenConvex.data(),yserieEigenConvex.data()+ yserieEigenConvex.size());
    std::vector<double> yseriesFull(yserieEigenFull.data(),yserieEigenFull.data()+ yserieEigenFull.size());
    
    std::vector<double> derivativeConcave(derivativeEigenConcave.data(),derivativeEigenConcave.data()+ derivativeEigenConcave.size());
    std::vector<double> derivativeConvex(derivativeEigenConvex.data(),derivativeEigenConvex.data()+ derivativeEigenConvex.size());
    
    
    std::vector<double> xSegmentsUpper(upper_hull.x.data(),upper_hull.x.data()+ upper_hull.x.size());
    std::vector<double> ySegmentsUpper(upper_hull.hx.data(),upper_hull.hx.data()+ upper_hull.hx.size());
    std::vector<double> xSegmentsLower(lower_hull.x.data(),lower_hull.x.data()+ lower_hull.x.size());
    std::vector<double> ySegmentsLower(lower_hull.hx.data(),lower_hull.hx.data()+ lower_hull.hx.size());
    
    string filename ="iter_"+ std::to_string(iter)+"_log_proposal_distribution_pair_" + std::to_string(idxFirstID)+ "_"+ std::to_string(idxSecondId)+ ".png";
    string title ="iter "+ std::to_string(iter)+" log proposal distribution for pair (" + std::to_string(idxFirstID)+ ", "+ std::to_string(idxSecondId)+ ")";
    plotter.plot2dSerieWithSegments(xseries,  yseriesFull,xSegmentsUpper,ySegmentsUpper,xSegmentsLower,ySegmentsLower, true, filename,title);
    
    
}
//----------------------------------------------------------------------------
GenotypeJCPairLogProposal::GenotypeJCPairLogProposal(int numSites,
                                                     int numStates, double Torigin, double delta, double theta, double pair_creation_time,
                                                     double time_left_child,
                                                     double time_right_child, const double* left_clv, const double* right_clv,
                                                     bool normalizedCLVs):
numSites(numSites),numStates(numStates), Torigin(Torigin),delta(delta),theta(theta),pair_creation_time(pair_creation_time),time_left_child(time_left_child),
time_right_child(time_right_child), normalizedCLVs(normalizedCLVs)
{
    double stat_prob = 1.0 / numStates;
    size_t span = numStates * numSites;
    
    
    num_concave_sites = 0;
    num_convex_sites = 0;
    int num_zero_sites = 0;
    oneMinusSumTermConcave.resize(numSites);
    oneMinusSumTermConvex.resize(numSites);
    
    firstTermConcave.resize(numSites);
    firstTermConvex.resize(numSites);
    for(size_t i= 0; i< numSites; ++i){
        
        
        double sum_for_site = 0.0;
        double sum_for_left_clv = 0.0;
        double sum_for_right_clv = 0.0;
        for(size_t j= 0; j< numStates; ++j){
            sum_for_site+= stat_prob * left_clv[j]* right_clv[j];
            sum_for_left_clv+= stat_prob*left_clv[j];
            sum_for_right_clv+= stat_prob*right_clv[j];
            // std::cout << " left_pclv " << left_clv[j] << std::endl;
            // std::cout << " right_pclv " << right_clv[j] << std::endl;
        }
        //assert(1-sum_for_site >= -numStates);
        //assert(1-sum_for_site <= numStates );
        if (normalizedCLVs){
            assert(abs(sum_for_left_clv -1.0)< 1e-4);
            assert(abs(sum_for_right_clv -1.0)< 1e-4);
            
            if (1-sum_for_site > 0){
                oneMinusSumTermConcave(num_concave_sites)= 1-sum_for_site;
                num_concave_sites++;
            }
            else if(1-sum_for_site < 0) {
                oneMinusSumTermConvex(num_convex_sites)=1-sum_for_site;
                num_convex_sites++;
            }
            else{
                
                num_zero_sites++;
            }
        }
        else{
//            std::cout << " sum_left*sum_right term site i "<< i <<" "<< sum_for_left_clv*sum_for_right_clv<< std::endl;
//            std::cout << " sum for site i " << i <<" "<<  sum_for_site<< std::endl;
//            std::cout << "one minus ratio for site i " << i <<" "<< 1- sum_for_site/(sum_for_left_clv*sum_for_right_clv) << std::endl;
//            std::cout << " coeff of exp site i " << i <<" "<<  sum_for_left_clv*sum_for_right_clv-sum_for_site<< std::endl;
            
            if (sum_for_site < sum_for_left_clv*sum_for_right_clv){
                
                oneMinusSumTermConcave(num_concave_sites)= (sum_for_left_clv*sum_for_right_clv)-sum_for_site;
                firstTermConcave(num_concave_sites) =sum_for_left_clv*sum_for_right_clv;
                
                assert( oneMinusSumTermConcave(num_concave_sites)>0);
                assert( firstTermConcave(num_concave_sites)>0);
                num_concave_sites++;
            }
            else if(sum_for_site > (sum_for_left_clv*sum_for_right_clv)) {
                oneMinusSumTermConvex(num_convex_sites)=(sum_for_left_clv*sum_for_right_clv)-sum_for_site;
                firstTermConvex(num_convex_sites) =sum_for_left_clv*sum_for_right_clv;
                
                assert( oneMinusSumTermConvex(num_convex_sites)<0);
                assert( firstTermConvex(num_convex_sites)>0);
                num_convex_sites++;
            }
            else{
                
                num_zero_sites++;
            }
            
            
        }
        assert(num_zero_sites==0);
        // std::cout << "Number sites with zero term"<< num_zero_sites << std::endl;
        
        left_clv+=numStates;
        right_clv+=numStates;
    }
    left_clv -= span;
    right_clv -= span;
    
}

GenotypeJCPairLogProposal::GenotypeJCPairLogProposal(int numSites,
                                                     int numStates, double Torigin, double delta, double theta, double pair_creation_time,
                                                     double time_left_child,
                                                     double time_right_child,std::vector<double>& oneMinusSumTermConcaveVector, std::vector<double>& oneMinusSumTermConvexVector):
numSites(numSites),numStates(numStates), Torigin(Torigin),delta(delta),theta(theta),pair_creation_time(pair_creation_time),time_left_child(time_left_child),
time_right_child(time_right_child)
{
    
    
    
    num_concave_sites = 0;
    num_convex_sites = 0;
    normalizedCLVs = true;
    
    num_concave_sites= oneMinusSumTermConcaveVector.size();
    num_convex_sites= oneMinusSumTermConvexVector.size();
    
    assert(num_concave_sites+num_convex_sites==numSites);
    
    oneMinusSumTermConcave= Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(oneMinusSumTermConcaveVector.data(), oneMinusSumTermConcaveVector.size());
    
    oneMinusSumTermConvex= Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(oneMinusSumTermConvexVector.data(), oneMinusSumTermConvexVector.size());
    
    
}
Eigen::VectorXd GenotypeJCPairLogProposal::clean( const Eigen::VectorXd& x) {
    // ensure x is below 700 (max value for which we can compute exp)
    Eigen::VectorXd result(x.size());
    result = x;
    for (int i = 0; i < result.size(); ++i) {
        if (result(i) > 700) {
            result(i) = 700;
        }
    }
    return result;
}
Eigen::VectorXd GenotypeJCPairLogProposal::h_concave(const Eigen::VectorXd& x) {
    Eigen::VectorXd result = clean(x);
    
    result = result.unaryExpr([this](double d) {
        return h_concave(d);
    });
    
    return result;
}
Eigen::VectorXd GenotypeJCPairLogProposal::h_convex(const Eigen::VectorXd& x) {
    Eigen::VectorXd result = clean(x);
    
    result = result.unaryExpr([this](double d) {
        return h_convex(d);
    });
    
    return result;
}

Eigen::VectorXd GenotypeJCPairLogProposal::h_prime_concave(const Eigen::VectorXd& x) {
    Eigen::VectorXd result = clean(x);
    
    
    result = result.unaryExpr([this](double d) {
        return h_prime_concave(d);
    });
    
    return result;
    
}

Eigen::VectorXd GenotypeJCPairLogProposal::h_prime_convex(const Eigen::VectorXd& x) {
    Eigen::VectorXd result = clean(x);
    
    result = result.unaryExpr([this](double d) {
        return h_prime_convex(d);
    });
    
    return result;
    
}

double GenotypeJCPairLogProposal::h_concave(double x) {
    std::vector<double> input = {x};
    
    double  result = 0.0;
    result+= -Population::FmodelTstandard(x, Torigin, delta, K)+ Population::FmodelTstandard(pair_creation_time, Torigin, delta,K);
    
    Eigen::VectorXd transformed;
    if  (normalizedCLVs){
        transformed = oneMinusSumTermConcave.head(num_concave_sites).unaryExpr([this,x](const double&  d) {
            double term = 1.0-exp(-2*x+time_left_child+time_right_child)*d;
            
            if (term < 1e-20){
                term = 1e-20;
                std::cout<<"term 0 in h_concave for x= " << x << std::endl;
            }
            return log(term);
        });
    }
    else {
        transformed = firstTermConcave.head(num_concave_sites) -exp(-2*x+time_left_child+time_right_child)*oneMinusSumTermConcave.head(num_concave_sites);
        
        transformed.unaryExpr([]( double d) {
            
            if (d < 1e-20){
                
                std::cout<<"term 0 in h_convex for x= " <<  std::endl;
                return log(1e-20);
            }
            else{
                return log(d);
                
            }
        });
        
    }
    assert(transformed.size()==num_concave_sites);
    result+= transformed.sum();
    
    return result;
}
Eigen::VectorXd  GenotypeJCPairLogProposal::minus_lambda(const Eigen::VectorXd& x) {
    
    Eigen::VectorXd result = clean(x);
    
    result = result.unaryExpr([this](double d) {
           double term =  -Population::FmodelTstandard(d, Torigin, delta, K)+ Population::FmodelTstandard(pair_creation_time, Torigin, delta,K);
            return term;
       });
       
    return result;
}
double GenotypeJCPairLogProposal::h_convex(double x) {
    std::vector<double> input = {x};
    
    double  result = 0.0;
    //    result+= -Population::FmodelTstandard(x, Torigin, delta, K)+ Population::FmodelTstandard(pair_creation_time, Torigin, delta,K);
    
    Eigen::VectorXd transformed;
    if  (normalizedCLVs){
        transformed= oneMinusSumTermConvex.head(num_convex_sites).unaryExpr([this,x](const double&  d) {
            double term =1.0-exp(-2*x+time_left_child+time_right_child)*d;
            if (term < 1e-20){
                term = 1e-20;
                std::cout<<"term 0 in h_convex for x= " << x << std::endl;
            }
            return log(term);
        });
    }
    else{
        transformed = firstTermConvex.head(num_convex_sites) -exp(-2*x+time_left_child+time_right_child)*oneMinusSumTermConvex.head(num_convex_sites);
        
        transformed.unaryExpr([]( double d) {
            
            if (d < 1e-20){
                
                std::cout<<"term 0 in h_convex for x= " <<  std::endl;
                return log(1e-20);
            }
            else{
                return log(d);
                
            }
        });
        
    }
    assert(transformed.size()==num_convex_sites);
    result+= transformed.sum();
    
    return result;
}

double GenotypeJCPairLogProposal::h_prime_concave(double x) {
    std::vector<double> input = {x};
    double  result = 0.0;
    result+= -2.0/Population::CalculateH(x, Torigin, delta,K);// this part has second derivative >0 then it  is convex
    
    Eigen::VectorXd transformed;
    if  (normalizedCLVs){
        
        transformed = oneMinusSumTermConcave.head(num_concave_sites).unaryExpr([this,x](const double&  d) {
            double term= exp(-2*x+time_left_child+time_right_child)*d;
            return (2.0* term  / (1-term));
        });
        
    }
    else{
        double term = exp(-2.0*x+time_left_child+time_right_child);
        transformed = 2.0*term*oneMinusSumTermConcave.head(num_concave_sites).array()
                           /(firstTermConcave.head(num_concave_sites).array() -term*oneMinusSumTermConcave.head(num_concave_sites).array() );
//        transformed = 2.0
//        /(firstTermConcave.head(num_concave_sites).array()/(term*oneMinusSumTermConcave.head(num_concave_sites).array()) -1.0 );
        
    }
    
    assert(transformed.size()==num_concave_sites);
    result+= transformed.sum();
    
    
    return result;
}
double GenotypeJCPairLogProposal::h_prime_convex(double x) {
    std::vector<double> input = {x};
    double  result = 0.0;
    // result+= -2.0/Population::CalculateH(x, Torigin, delta,K)  ;
    
    Eigen::VectorXd transformed;
    if  (normalizedCLVs){
        transformed= oneMinusSumTermConvex.head(num_convex_sites).unaryExpr([this,x](const double& d) {
            double term= exp(-2*x+time_left_child+time_right_child)*d;
            return (2.0* term  / (1-term));
        });
    }
    else{
        double term = exp(-2.0*x+time_left_child+time_right_child);
        transformed = 2.0*term*oneMinusSumTermConvex.head(num_convex_sites).array()
                    /(firstTermConvex.head(num_convex_sites).array() -term*oneMinusSumTermConvex.head(num_convex_sites).array() );
       // transformed = 2.0
       //      /(firstTermConvex.head(num_convex_sites).array()/(term*oneMinusSumTermConvex.head(num_convex_sites).array()) -1.0 );

    }
    assert(transformed.size()==num_convex_sites);
    result+= transformed.sum();
    
    return result;
}
void GenotypeJCPairLogProposal::plot(GNUPlotter &plotter, int num_grid_points, int iter,  int idxFirstID, int idxSecondId ){
    
    Eigen::VectorXd xserieEigen = Eigen::VectorXd::LinSpaced(num_grid_points,pair_creation_time,0.999*Torigin);
    std::vector<double> xseries(xserieEigen.data(),xserieEigen.data()+ xserieEigen.size() );
    Eigen::VectorXd yserieEigenConcave = h_concave(xserieEigen);
    Eigen::VectorXd yserieEigenConvex = h_convex(xserieEigen);
    Eigen::VectorXd yMinusLambda = minus_lambda(xserieEigen);
    Eigen::VectorXd yserieEigenFull= yserieEigenConcave+yserieEigenConvex;
    
    Eigen::VectorXd derivativeEigenConcave = h_prime_concave(xserieEigen);
    Eigen::VectorXd derivativeEigenConvex = h_prime_convex(xserieEigen);
    
    std::vector<double> yseriesConcave(yserieEigenConcave.data(),yserieEigenConcave.data()+ yserieEigenConcave.size());
    std::vector<double> yseriesConvex(yserieEigenConvex.data(),yserieEigenConvex.data()+ yserieEigenConvex.size());
    std::vector<double> yseriesFull(yserieEigenFull.data(),yserieEigenFull.data()+ yserieEigenFull.size());
     std::vector<double> yseriesMinusLambda(yMinusLambda.data(),yMinusLambda.data()+ yMinusLambda.size());
    
    std::vector<double> derivativeConcave(derivativeEigenConcave.data(),derivativeEigenConcave.data()+ derivativeEigenConcave.size());
    std::vector<double> derivativeConvex(derivativeEigenConvex.data(),derivativeEigenConvex.data()+ derivativeEigenConvex.size());
    
    plotter.plot2dSerie(xseries,  yseriesMinusLambda, true, "iter_"+ std::to_string(iter)+"_log_lambda_pair_" + std::to_string(idxFirstID)+ "_"+ std::to_string(idxSecondId)+ ".png"," Log distribution of coal. time for pair "+std::to_string(idxFirstID)+ ", "+ std::to_string(idxSecondId) );
    
    plotter.plot2dSerie(xseries,  yseriesConcave, true, "iter_"+ std::to_string(iter)+"_concave_part_" + std::to_string(idxFirstID)+ "_"+ std::to_string(idxSecondId)+ ".png","concave part" );
    plotter.plot2dSerie(xseries,  yseriesConvex, true, "iter_"+ std::to_string(iter)+ "_convex_part_" + std::to_string(idxFirstID)+ "_"+ std::to_string(idxSecondId)+ ".png","convex part" );
    plotter.plot2dSerie(xseries,  yseriesFull, true, "iter_"+ std::to_string(iter)+"_log_proposal_pair_" + std::to_string(idxFirstID)+ "_"+ std::to_string(idxSecondId)+ ".png"," log proposal distrib" );
    
    
}
