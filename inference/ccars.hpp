//
//  ccars.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 12/10/2021.
//

#ifndef ccars_hpp
#define ccars_hpp

#include "data_types.hpp"
#include "random.h"
#include <Eigen/Dense>

using namespace std;
class GNUPlotter;
class CCLogDensity {
    
public:
    
    
    virtual Eigen::VectorXd h_concave(const Eigen::VectorXd& x) = 0;
    virtual Eigen::VectorXd h_convex(const Eigen::VectorXd& x) = 0;
    virtual Eigen::VectorXd h_prime_concave(const Eigen::VectorXd&  x) = 0;
    virtual Eigen::VectorXd h_prime_convex(const Eigen::VectorXd& x) = 0;
    virtual double h_concave(double x) = 0;
    virtual double h_convex(double x) = 0;
    virtual double h_prime_concave(double x) = 0;
    virtual double h_prime_convex(double x) = 0;
    virtual void plot(GNUPlotter &plotter, int num_grid_points,int iter,  int idxFirstID, int idxSecondId )= 0;
};

class GenotypeJCPairLogProposal : public CCLogDensity {
    
public:
    int numSites;
    int numStates;
    double Torigin;
    double delta;
    double theta;
    double pair_creation_time;
    double time_left_child;
    double time_right_child;
    bool normalizedCLVs;
    
    double K = 0.8;
    Eigen::VectorXd oneMinusSumTermConcave;
    Eigen::VectorXd oneMinusSumTermConvex;
    Eigen::VectorXd firstTermConcave;
    Eigen::VectorXd firstTermConvex;
    int num_concave_sites ;
    int num_convex_sites ;
    
    GenotypeJCPairLogProposal( int numSites,
                              int numStates,double Torigin, double delta, double theta, double pair_creation_time,
                              double time_left_child,
                              double time_right_child, const double* left_clv, const double* right_clv,bool normalizedCLVs);
    
    GenotypeJCPairLogProposal(int numSites,
                              int numStates, double Torigin, double delta, double theta, double pair_creation_time,
                              double time_left_child,
                              double time_right_child,std::vector<double>& oneMinusSumTermConcaveVector, std::vector<double>& oneMinusSumTermConvexVector);
    
    Eigen::VectorXd clean(const Eigen::VectorXd& x);
    
    Eigen::VectorXd h_concave(const Eigen::VectorXd& x);
    Eigen::VectorXd h_convex(const Eigen::VectorXd& x);
    Eigen::VectorXd h_prime_concave(const Eigen::VectorXd& x);
    Eigen::VectorXd h_prime_convex(const Eigen::VectorXd& x);
    double h_concave(double x);
    double h_convex(double x);
    double h_prime_concave(double x);
    double h_prime_convex(double x);
    
    void plot(GNUPlotter &plotter, int num_grid_points,int iter, int idxFirstID, int idxSecondId );
    Eigen::VectorXd  minus_lambda(const Eigen::VectorXd& x);
};
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
    std::vector<HullSegment> segments;
    double log_cu;
    double xlb;
    double xrb;
    double hx_xlb;
    double hx_xrb;
    double hpx_xlb;
    double hpx_xrb;
    
    Eigen::VectorXd x;
    Eigen::VectorXd hx;
    Eigen::VectorXd tangents_slopes;
    Eigen::VectorXd x_intercepts;
    Eigen::VectorXd hx_intercepts;
    
    UpperHull() {}; // Default constructor
    
    UpperHull(const Eigen::VectorXd& x,const Eigen::VectorXd& hx,const  Eigen::VectorXd& hpx, double _xlb, double _xrb, double hx_xlb, double hx_xrb, double hpx_xlb, double hpx_xrb) ;
    
    
    double get_hu(double x) ;
    void add_segment(double x, double hx, double hpx);
    void add_segment(double x, double hx, double hpx, int pos);
    double sample(gsl_rng *rng) ;
    void update_cumulative_from_pos(int pos);
    void update_intercepts_from_pos(int pos);
    void get_intercepts(Eigen::VectorXd &x_intercept, Eigen::VectorXd &hx_intercept);
    int get_point_position(double new_x);
    
};
class LowerHull {
public:
    //    std::vector<double> x;
    //    std::vector<double> hx;
    //    std::vector<double> hpx;
    
    Eigen::VectorXd x;
    Eigen::VectorXd hx;
    Eigen::VectorXd slopes;
    int n_points;
    
    LowerHull() {}; // Default constructor
    
    LowerHull(const Eigen::VectorXd&  x,const Eigen::VectorXd&  hx) ;
    
    int get_point_position(double new_x);
    double get_slope(int segment_idx);
    double get_hl(double _x);
    void add_segment(double new_x, double new_hx, double new_hpx);
    void add_segment(double new_x, double new_hx, double new_hpx, int position);
    Eigen::VectorXd get_h_intercepts( const Eigen::VectorXd& x_intercepts);
    
};
class CCUpperHull{
public:
    UpperHull *upper_hull_concave;
    LowerHull *upper_hull_convex;
    UpperHull cc_upper_hull;
    
    Eigen::VectorXd x;
    Eigen::VectorXd hx;
    Eigen::VectorXd hpx;
    double xlb;
    double xrb;
    double maxCumulative;
    Eigen::VectorXd cumulative;
    CCUpperHull(){};
    
    CCUpperHull(UpperHull &upper_hull_concave, LowerHull &upper_hull_convex, double xlb, double xrb );
    
    
    double get_hu(double x) ;
    void add_segment(double x, double hx_concave, double hpx_concave,double hx_convex, double hpx_convex );
    void add_segment(double x, double hx_concave, double hpx_concave,double hx_convex, double hpx_convex, int pos );
    void init(UpperHull &upper_hull_concave, LowerHull &upper_hull_convex,
              Eigen::VectorXd& x_temp,Eigen::VectorXd& hx_temp,Eigen::VectorXd& hpx_temp);
    
    double sample(gsl_rng *rng);
    void update_cumulative();
    int get_point_position(double new_x);
    
    
};
class CCLowerHull{
public:
    LowerHull *lower_hull_concave;
    UpperHull *lower_hull_convex;
    LowerHull cc_lower_hull;
    
    Eigen::VectorXd x;
    Eigen::VectorXd hx;
    Eigen::VectorXd hpx;
    double xlb;
    double xrb;
    double minCumulative;
    Eigen::VectorXd cumulative;
    CCLowerHull(){};
    
    CCLowerHull(LowerHull &lower_hull_concave, UpperHull &lower_hull_convex, double xlb, double xrb );
    
    
    double get_hu(double x) ;
    void add_segment(double x, double hx_concave, double hpx_concave,double hx_convex, double hpx_convex,  int position );
    
    void add_segment(double new_x, double new_hx, double new_hpx, int position);
    
    double sample(gsl_rng *rng);
    void update_cumulative();
    int get_point_position(double new_x);
    void init(LowerHull &lower_hull_concave, UpperHull &lower_hull_convex,   Eigen::VectorXd& x_temp,Eigen::VectorXd& hx_temp,Eigen::VectorXd& hpx_temp);
    
    
};
class CCARS {
    // CCAR sampling
public:
    LowerHull upper_hull_convex;
    UpperHull lower_hull_convex;
    
    UpperHull upper_hull_concave;
    LowerHull lower_hull_concave;
    
    CCUpperHull upper_hull;
    CCLowerHull lower_hull;
    
    CCLogDensity*  log_density;
    
    double xlb;
    double xrb;
    int max_points;
    Eigen::VectorXd ratioAreas;
    
    const static int MAX_REJECTIONS = 500;
    
    CCARS() {}; // Default constructor
    
    CCARS(CCLogDensity* const log_density,const Eigen::VectorXd&  x, double xlb, double xrb, int max_points);
    
    std::vector<double> sample(gsl_rng *rng, unsigned N, CCLogDensity* const log_density) ;
    double approximateLogIntegral(gsl_rng *rng,  double  max_number_internal_points);
    
    double sampleNewX(gsl_rng *rng,double x1, double x2, double y1, double y2, double m, double &logCumArea );
    void addPoint(double new_x, double hconcave_new_x, double hconvex_new_x, int pos_new_x);
    void addPoint(double new_x);
    double sample(gsl_rng *rng,  double& y, double& logComplementArea);
    
    static  double logSumExp(const Eigen::VectorXd& values);
    static double logSumExp(double u, double v);
    static double logDiffExp(double u, double v);
    void plot( GNUPlotter &plotter, int num_grid_points,int iter,  int idxFirstID, int idxSecondId);
};
#endif /* ccars_hpp */
