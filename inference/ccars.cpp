//
//  ccars.cpp
//
//  Code from
// https://github.com/hunzikp/arscpp
// Adaptive Rejection Sampler

//Philipp Hunziker June 3, 2018

#include "ccars.hpp"
//#include "data_types.hpp"
#include "data_utils.hpp"
#include "random.h"
#include "random.h"
#include <Eigen/Dense>

#include <algorithm>
using namespace std;
double logSumExp(double u, double v) {
    double rv  = max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, v)));
    return rv;
}
double logDiffExp(double u, double v) {
    double rv = max(u, v) + log(exp(u - max(u, v)) - exp(v - max(u, v)));
    return rv;
}


struct HullSegment {
public:
    double x;//x of tangent point where the segment and log(g(x)) meet
    double hx;//log(g(x)) at x
    double hpx;//derivative of log(g(x)) at x
    double z_left;//left abscissae of segment
    double z_right;//right abscissae of segment
    double hu_left;//ordinate of left end point of  segment
    double hu_right;//ordinates of right end point of segment
    double scum_left;//cumulative distribution of the area under upper hull
    double scum_right;
};

class UpperHull {
public:
    int n_segments;
    vector<HullSegment> segments;
    double log_cu;
    double xlb;
    double xrb;
    
    UpperHull() {}; // Default constructor
    
    UpperHull(std::vector<double> x, std::vector<double> hx, std::vector<double> hpx, double _xlb, double _xrb) {
        // x assumed sorted!
        
        // No of segments
        n_segments = x.size();
        
        // Set bounds
        xlb = _xlb;
        xrb = _xrb;
        
        // Ensure left (right) most point is left (right) of mode
        if (hpx.at(0) < 0) {
            std::cout<< "Smallest starting point is right of mode."<<std::endl;
        }
        if (hpx.at(hpx.size()-1) > 0) {
            std::cout<<"Largest starting point is left of mode."<<std::endl;
        }
        
        // Initialize segments
        for (int i = 0; i < n_segments; ++i) {
            
            HullSegment segment;
            segment.x = x.at(i);
            segment.hx = hx.at(i);
            segment.hpx = hpx.at(i);
            
            segments.push_back(segment);
        }
        
        // Update z and scum
        update_z();
        update_scum();
        
    };
    
    void update_z() {
        
        for (int i = 0; i < n_segments; ++i) {
            
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
                segment.z_right = segment.x + (num/denom);
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
    
    void update_scum() {
        
        // Compute normalizer
        vector<double> log_cu_vec(n_segments);
        log_cu = 0;
        for (int i = 0; i < n_segments; ++i) {
            HullSegment segment;
            segment = segments.at(i);
            // double cu_i = (1/segment.hpx)*(exp(segment.hu_right) - exp(segment.hu_left));
            double log_cu_i;
            if (segment.hpx > 0) {
                log_cu_i = -log(segment.hpx) + logDiffExp(segment.hu_right, segment.hu_left);
            } else {
                log_cu_i = logDiffExp(segment.hu_left - log(-segment.hpx), segment.hu_right - log(-segment.hpx));
            }
            
            if (i == 0) {
                log_cu = log_cu_i;
            } else {
                log_cu = logSumExp(log_cu, log_cu_i);
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
    
    double sample(gsl_rng *rng) {
        
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
                x_sample = s.z_right + (1/s.hpx)*logSumExp(0, log(s.hpx*(u-s.scum_right)) + log_cu - s.hu_right);
            } else {
                x_sample = s.z_right + (1/s.hpx)*logDiffExp(0, log(-1*s.hpx*(u - s.scum_right)) + log_cu - s.hu_right);
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
                x_sample = s.z_left + (1/s.hpx)*logSumExp(0, log(s.hpx) + log_cu + log(u - s.scum_left) - s.hu_left);
            } else {
                x_sample = s.z_left + (1/s.hpx)*logDiffExp(0, log(-s.hpx) + log_cu + log(u - s.scum_left) - s.hu_left);
            }
        }
        
        return x_sample;
    };
    
    void add_segment(double x, double hx, double hpx) {
        
        HullSegment segment;
        segment.x = x;
        segment.hx = hx;
        segment.hpx = hpx;
        
        // Determine segment position
        int iter = 0;
        while(iter < n_segments) {
            if (x < segments.at(iter).x) {
                break;
            }
            iter++;
        }
        
        // Insert segment
        segments.insert(segments.begin()+iter, segment);
        n_segments = segments.size();
        
        // Update z and scum
        update_z();
        update_scum();
        
    };
    
    double get_hu(double x) {
        
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
};

class LowerHull {
public:
    vector<double> x;
    vector<double> hx;
    vector<double> hpx;
    int n_points;
    
    LowerHull() {}; // Default constructor
    
    LowerHull(std::vector<double> x, std::vector<double> hx, std::vector<double> hpx):x(x),hx(hx), hpx(hpx) {
        // x assumed sorted!
        
        n_points = x.size();
        
    };
    
    int get_point_position(double new_x){
        
        // Determine point position
        int iter = 0;
        while (iter < n_points) {
            if (new_x < x.at(iter)) {
                break;
            }
            iter++;
        }
        return iter;
    }
    void add_segment(double new_x, double new_hx, double new_hpx) {
        
        // Determine point position
        int iter =get_point_position( new_x);
        
        // Assign
        x.insert(x.begin()+iter, new_x);
        hx.insert(hx.begin()+iter, new_hx);
        hpx.insert(hpx.begin()+iter, new_hpx);
        
        n_points = x.size();
    }
    
    double get_hl(double _x) {
        
        // Determine point position
        int iter = 0;
        while (iter < n_points) {
            if (_x < x.at(iter)) {
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
            double x_left = x.at(iter-1);
            double x_right = x.at(iter);
            double hx_left = hx.at(iter-1);
            double hx_right = hx.at(iter);
            double d = (hx_right-hx_left)/(x_right-x_left);
            rv = hx_left + (_x-x_left)*d;
        }
        
        return rv;
    }
    double get_slope(int segment_idx) {
        
         double x_left = x.at(segment_idx);
         double x_right = x.at(segment_idx+1);
         double hx_left = hx.at(segment_idx);
         double hx_right = hx.at(segment_idx+1);
         double d = (hx_right-hx_left)/(x_right-x_left);
         return d;
    }
    
};
class CCUpperHull{
public:
    UpperHull *upper_hull_concave;
    LowerHull *upper_hull_convex;
    UpperHull cc_upper_hull;
    double xlb;
    double xrb;
    
    CCUpperHull(){};
    
    CCUpperHull(UpperHull &upper_hull_concave, LowerHull &upper_hull_convex, double xlb, double xrb ) :
      upper_hull_concave(&upper_hull_concave),
      upper_hull_convex(&upper_hull_convex),
      xlb(xlb),
      xrb(xrb)
    {
        for (int i = 0; i < upper_hull_concave.n_segments; ++i) {
                 
                 HullSegment segment;
                 segment.x = upper_hull_concave.segments.at(i).x;
                 segment.z_left = upper_hull_concave.segments.at(i).z_left;
                 segment.z_right = upper_hull_concave.segments.at(i).z_right;
            
                 segment.hu_left = upper_hull_concave.segments.at(i).hu_left +upper_hull_convex.hx.at(i);
            
                 segment.hu_right = upper_hull_concave.segments.at(i).hu_right +upper_hull_convex.hx.at(i+1);
            
                 segment.hpx = upper_hull_concave.segments.at(i).hpx + upper_hull_convex.get_slope(i) ;
                 
                 cc_upper_hull.segments.push_back(segment);
             }
        
        cc_upper_hull.update_scum();
        
    }
    double sample(gsl_rng *rng) {
        
        return(cc_upper_hull.sample(rng));
        
    }
    double get_hu(double x) {
        
       return(cc_upper_hull.get_hu(x));
    }
    
    void add_segment(double x, double hx_concave, double hpx_concave,double hx_convex, double hpx_convex ) {
        
        
        upper_hull_concave->add_segment(x, hx_concave,  hpx_concave);
        
        int pos =  upper_hull_convex->get_point_position( x);
        
        upper_hull_convex->add_segment(x, hx_convex,  hpx_convex);
        
        cc_upper_hull.add_segment(x,hx_concave+ hx_convex,  hpx_concave + upper_hull_convex->get_slope(pos) );
        
        cc_upper_hull.update_scum();
        
    }

};
//___________________________________________________________________________
class CCARS {
    // CCAR sampling
public:
    LowerHull upper_hull_convex;
    UpperHull lower_hull_convex;
    
    UpperHull upper_hull_concave;
    LowerHull lower_hull_concave;
    
    CCUpperHull upper_hull;
    
  

    
    double xlb;
    double xrb;
    int max_points;
    
    const static int MAX_REJECTIONS = 500;
    
    CCARS() {}; // Default constructor
    
    CCARS(CCLogDensity* const log_density, std::vector<double>  x, double xlb, double xrb, int max_points) :
    xlb(xlb),
    xrb(xrb),
    max_points(max_points)
    {
        
        // Check bounds
        if (xlb >= xrb) {
            std::cout<< "Upper bound is not larger than lower bound."<<std::endl;
        }
        
        // We need at least two starting points
        if (x.size() < 2) {
            std::cout<< "At least two starting points required."<<std::endl ;
        }
        
        // Order x
        std::sort (x.begin(), x.end());
        
        // Check points in bounds
        if (x.at(0) < xlb || x.at(x.size()-1) > xrb) {
            std::cout<<"Starting point out of bound."<<std::endl;
        }
        
        // Evaluate funcs
        std::vector<double> hx_concave = log_density->h_concave(x);
        std::vector<double>  hpx_concave = log_density->h_prime_concave(x);
        
        std::vector<double> hx_convex = log_density->h_convex(x);
        std::vector<double>  hpx_convex = log_density->h_prime_convex(x);
        
        // Try to make starting points valid (if they're invalid)
        int max_tries = 10;
        if (hpx_concave.at(0) <= 0) {
            // Move left-most point left until valid
            
            double hpx_left = hpx_concave.at(0);
            double x_left = x.at(0);
            int tries = 0;
            while (hpx_left <= 0 && tries < max_tries) {
                
                if (isfinite(xlb)) {
                    x_left -= (x_left-xlb)/2; // Move half-way to the limit
                } else {
                    x_left -= pow(2, tries); // Move left by power of 2
                }
                
                hpx_left = log_density->h_prime_concave(x_left);
                
                tries++;
            }
            
            if (tries < max_tries) {
                hpx_concave.at(0) = hpx_left;
                x.at(0) = x_left;
                hx_concave.at(0) = log_density->h_concave(x_left);
            } else {
                std::cout<< "Could not find valid lower starting point."<<std::endl;
            }
        }
        
        int last_ind = hpx_concave.size() - 1;
        if (hpx_concave.at(last_ind) >= 0) {
            // Move right-most point right until valid
            
            double hpx_right = hpx_concave.at(last_ind);
            double x_right = x.at(last_ind);
            int tries = 0;
            while (hpx_right >= 0 && tries < max_tries) {
                
                if (isfinite(xrb)) {
                    x_right += (xrb - x_right)/2; // Move half-way to the limit
                } else {
                    x_right += pow(2, tries); // Move right by power of 2
                }
                
                hpx_right = log_density->h_prime_concave(x_right);
                
                tries++;
            }
            
            if (tries < max_tries) {
                hpx_concave.at(last_ind) = hpx_right;
                x.at(last_ind) = x_right;
                hx_concave.at(last_ind) = log_density->h_concave(x_right);
            } else {
                // std::cout << "x candidates: " << x << std::endl;
                std::cout <<"Could not find valid upper starting point."<< std::endl;
            }
        }
        
        // Create the hull
        upper_hull_concave = UpperHull(x, hx_concave, hpx_concave, xlb, xrb);
        lower_hull_concave = LowerHull(x, hx_concave, hpx_concave);
        
        lower_hull_convex = UpperHull(x, hx_concave, hpx_concave, xlb, xrb);
        upper_hull_convex = LowerHull(x, hx_concave, hpx_concave);
        
        upper_hull = CCUpperHull(upper_hull_concave, upper_hull_convex, xlb, xrb);
        
    };
    
    std::vector<double> sample(gsl_rng *rng, unsigned N, CCLogDensity* const log_density) {
        
        vector<double> samples;
        int rejections = 0;
        
        while(samples.size() < N) {
            
            double x_sample = upper_hull.sample(rng);
            double u = Random::randomUniformFromGsl2(rng);
            
           
            //if (u < exp( upper_hull.get_hu(x_sample) -hx_concave-hx_convex))
            if (u < exp(lower_hull_concave.get_hl(x_sample)+lower_hull_convex.get_hu(x_sample)- upper_hull.get_hu(x_sample))) {
                // Accept!
                samples.push_back(x_sample);
                rejections = 0;
            } else {
                
                double hx_concave = log_density->h_concave(x_sample);
                double hx_convex = log_density->h_convex(x_sample);


                if (u < exp(hx_concave+hx_convex - upper_hull.get_hu(x_sample))) {
                    // Accept!
                    samples.push_back(x_sample);
                    rejections = 0;
                } else {
                    // Reject!
                    rejections++;
                }
                
                // Add hull segment
                int points = lower_hull_concave.x.size();
                if (points < max_points) {
                    double hpx_concave = log_density->h_prime_concave(x_sample);
                    double hpx_convex = log_density->h_prime_convex(x_sample);
                    
                    upper_hull.add_segment(x_sample,  hx_concave,  hpx_concave, hx_convex,  hpx_convex );
    
                }
            }
            
            if (rejections > MAX_REJECTIONS) {
                std::cout << "Warning: Maximum number of rejections reached. Returning zero sample." <<  std::endl;
                samples.push_back(0);
            }
            
        }
        
        return samples;
    };
};
//----------------------------------------------------------------------------
GenotypeJCPairLogProposal::GenotypeJCPairLogProposal(int numSites,
                                                     int numStates, double Torigin, double delta, double theta, double pair_creation_time,
                                                     double time_left_child,
                                                     double time_right_child, double* left_clv, double* right_clv):
numSites(numSites),numStates(numStates), Torigin(Torigin),delta(delta),theta(theta),pair_creation_time(pair_creation_time),time_left_child(time_left_child),
time_right_child(time_right_child),
left_clv(left_clv),right_clv(right_clv)
{
    double stat_prob = 1.0 / numStates;
    size_t span = numStates * numSites;
    
    
    for(size_t i= 0; i< numSites; ++i){
        
        
        double sum_for_site = 0.0;
        for(size_t j= 0; j< numStates; ++j){
            sum_for_site+= stat_prob * left_clv[j]* right_clv[j];
            std::cout << " left_pclv " << left_clv[j] << std::endl;
            std::cout << " right_pclv " << right_clv[j] << std::endl;
        }
        if (1-sum_for_site >= 0)
            oneMinusSumTermConcave.push_back(1-sum_for_site);
        else
            oneMinusSumTermConvex.push_back(1-sum_for_site);
        
        left_clv+=numStates;
        right_clv+=numStates;
    }
    left_clv -= span;
    right_clv -= span;
    
}

std::vector<double> GenotypeJCPairLogProposal::clean(std::vector<double> x) {
    // ensure x is below 700 (max value for which we can compute exp)
    for (int i = 0; i < x.size(); ++i) {
        if (x.at(i) > 700) {
            x.at(i) = 700;
        }
    }
    return x;
}
std::vector<double> GenotypeJCPairLogProposal::h_concave(std::vector<double> x) {
    x = clean(x);
    std::vector<double> result(x.size());
    
    std::transform(x.begin(), x.end(), std::back_inserter(result),
                   [this](const double& c){ return h_concave(c); });
    
    return result;
}
std::vector<double> GenotypeJCPairLogProposal::h_convex(std::vector<double> x) {
    x = clean(x);
    std::vector<double> result(x.size());
    
    std::transform(x.begin(), x.end(), std::back_inserter(result),
                   [this](const double& c){ return h_convex(c); });
    
    return result;
}

std::vector<double> GenotypeJCPairLogProposal::h_prime_concave(std::vector<double> x) {
    x = clean(x);
    std::vector<double> result(x.size());
    
    std::transform(x.begin(), x.end(), std::back_inserter(result),
                   [this](const double& c){ return h_prime_concave(c); });
    return result;
}

std::vector<double> GenotypeJCPairLogProposal::h_prime_convex(std::vector<double> x) {
    x = clean(x);
    std::vector<double> result(x.size());
    
    std::transform(x.begin(), x.end(), std::back_inserter(result),
                   [this](const double& c){ return h_prime_convex(c); });
    return result;
}

double GenotypeJCPairLogProposal::h_concave(double x) {
    std::vector<double> input = {x};
    
    double  result = 0.0;
    result+= -Population::FmodelTstandard(x, Torigin, delta, K)+ Population::FmodelTstandard(pair_creation_time, Torigin, delta,K);
    size_t numSitesConcave = oneMinusSumTermConcave.size();
    for(size_t i = 0; i < numSitesConcave ; ++i){
        result+=log(1-exp(-2*x+time_left_child+time_right_child)*oneMinusSumTermConcave[i]);
        
    }
    return result;
}
double GenotypeJCPairLogProposal::h_convex(double x) {
    std::vector<double> input = {x};
    
    double  result = 0.0;
    result+= -Population::FmodelTstandard(x, Torigin, delta, K)+ Population::FmodelTstandard(pair_creation_time, Torigin, delta,K);
    size_t numSitesConvex = oneMinusSumTermConvex.size();
    for(size_t i = 0; i < numSitesConvex ; ++i){
        result+=log(1-exp(-2*x+time_left_child+time_right_child)*oneMinusSumTermConvex[i]);
        
    }
    return result;
}

double GenotypeJCPairLogProposal::h_prime_concave(double x) {
    std::vector<double> input = {x};
    double  result = 0.0;
    result+= -2.0/Population::CalculateH(x, Torigin, delta,K);// this part has second derivative <0 then it  is concave
    size_t numSitesConcave = oneMinusSumTermConcave.size();
    for(size_t i = 0; i < numSitesConcave ; ++i){
        double term= exp(-2*x+time_left_child+time_right_child)*oneMinusSumTermConcave[i];
        result+= 2.0* term  / (1-term);
        
    }
    return result;
}
double GenotypeJCPairLogProposal::h_prime_convex(double x) {
    std::vector<double> input = {x};
    double  result = 0.0;
    //result+= -2.0/Population::CalculateH(x, Torigin, delta,K);
    size_t numSitesConvex = oneMinusSumTermConvex.size();
    for(size_t i = 0; i < numSitesConvex ; ++i){
        double term= exp(-2*x+time_left_child+time_right_child)*oneMinusSumTermConvex[i];
        result+= 2.0* term  / (1-term);
        
    }
   
    return result;
}
