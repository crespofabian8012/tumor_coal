//
//  gnu_plotter.cpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 14/09/2021.
//

#include "gnu_plotter.hpp"
#include <Eigen/Core>
//using namespace boost;
//using namespace boost::iostreams;
GNUPlotter::GNUPlotter(){
    
    gp = new Gnuplot("/usr/local/bin/gnuplot");
    
}
void GNUPlotter::setToSavePng(std::string path){
    *gp << "set terminal png\n";
    *gp << "set output 'my_graph_1.png'\n";
};
int GNUPlotter::plot2dSeveralSeries(std::vector<double> &x, std::vector<std::vector<double>> &yseries ){
    
    int result = 0;
    if (yseries.size()==0)
        return result;
    if (x.size()!=yseries[0].size())
    {
        std::cout << " The lengt of x series is not equal to the y series"  << std::endl;
        result = 1;
        return result;
    }
    
    auto plots = (*gp).splotGroup();
    
    Eigen::Map<Eigen::VectorXd> xSeries(&x[0], x.size());
    *gp << "set xrange  [" << xSeries.minCoeff()<< ":" << xSeries.maxCoeff()<<"]\n";
    
    *gp << "plot '-' with vectors title 'logLikelihoodRatio'\n";
    (*gp).send(yseries);
    
    
    return 0;
}
int GNUPlotter::plot2dSeveralSeries(std::vector<std::vector<double>> &xseries, std::vector<std::vector<double>> &yseries, std::vector<std::string> &labels, bool savePng, std::string filePath, double verticalLine1, double verticalLine2 ){
    
    int result = 0;
    if (yseries.size()==0)
        return result;
    
    if (savePng){
        *gp << "set terminal png size 900,900\n";
        *gp << "set output  '"<< filePath << "'\n";
    }
    else{
        *gp << "set term qt persist size 900,900\n";
        
    }
    *gp << "set title  'Log Likelihood ratios for differents pairs'"<<"\n";
    *gp << "set xlabel 'time'"<<"\n";
    *gp <<  "set ylabel 'Log Likelihood Ratios'"<<"\n";
    
    
    *gp <<"set key outside top right\n";
    
    size_t nYSeries = yseries.size();
    size_t nXSeries = xseries.size();
    assert(nXSeries == nYSeries);
    std::vector<std::vector<std::pair<double, double>>> xy_pts_A;
    std::vector<std::pair<double, double>> temp;
    int i = 0;
    double y;
    double x;
    
    bool include = true;
    i= 0;
    double maxAxisValue = 10000;
    double xmax = -maxAxisValue;
    double xmin = maxAxisValue;
    double ymax = -maxAxisValue;
    double ymin = maxAxisValue;
    auto plots = (*gp).plotGroup();
    size_t size;
    for(size_t j=0 ; j< nYSeries; ++j) {
        include = true;
        
        i= 0;
        temp.clear();
        size = xseries[j].size();
        
        for(size_t i=0 ; i< size; ++i) {
            y= yseries[j].at(i);
            x= xseries[j].at(i);
            
            if (y< -maxAxisValue){
                include = false;
                break;
            }
            if (y > ymax)
                ymax= y;
            if (y < ymin)
                ymin= y;
            if (x > xmax)
                xmax= x;
            if (x < xmin)
                xmin= x;
            
            // std::cout << "x "<< x<< " y "<< y  << std::endl;
            temp.push_back(std::make_pair(x, y));
        }
        if (include){
            std::string cmd = "with points title 'pair" + labels[j]+"'";
            plots.add_plot1d(temp, cmd);
            
            xy_pts_A.push_back(temp);
        }
    }
    
    *gp << "set xrange  [" << xmin<< ":" << xmax<<"]\n";
    *gp << "set yrange  [" << ymin<< ":" << ymax<<"]\n";
    
    *gp << "set arrow from "<< verticalLine1<< ","<< ymin <<" to "<< verticalLine1 << "," <<ymax <<" nohead lc rgb \'black\'\n";
    *gp << "set arrow from "<< verticalLine2<< ","<< ymin <<" to "<< verticalLine2 << "," <<ymax <<" nohead lc rgb \'red\'\n";
    *gp << plots;
    
    return 0;
}

int GNUPlotter::plot2dSerie(std::vector<double> &xseries,  std::vector<double> &yseries, std::vector<std::string> &labels,bool doSavePng, std::string fileNameIncrements, double verticalLine1   ){
    
    int result = 0;
    if (yseries.size()==0)
        return result;
    
    if (doSavePng){
        *gp << "set terminal png size 900,900\n";
        *gp << "set output  '"<< fileNameIncrements << "'\n";
    }
    else{
        *gp << "set term qt persist size 900,900\n";
        
    }
    
    *gp << "set title  'Norm of Log Likelihood ratios for differents increments'"<<"\n";
    *gp << "set xlabel 'Increments'"<<"\n";
    *gp <<  "set ylabel 'Norm'"<<"\n";
    
    
    //size_t nYSeries = yseries.size();
    std::vector<std::pair<double, double>> xy_pts_A;
    std::vector<std::pair<double, double>> temp;

    
    size_t nYSeries = yseries.size();
    size_t nXSeries = xseries.size();
    assert(nXSeries == nYSeries);
    
    double y;
    double x;
    
    bool include = true;

    double maxAxisValue = 10000;
     double xmax = -maxAxisValue;
     double xmin = maxAxisValue;
     double ymax = -maxAxisValue;
     double ymin = maxAxisValue;
    auto plots = (*gp).plotGroup();
    size_t size =nXSeries ;
    
    
    for(size_t i=0 ; i< size; ++i) {
        y= yseries.at(i);
        x= xseries.at(i);
        
        if (y< -maxAxisValue){
            include = false;
            break;
        }
        if (y > ymax)
            ymax= y;
        if (y < ymin)
            ymin= y;
        if (x > xmax)
            xmax= x;
        if (x < xmin)
            xmin= x;
        
        // std::cout << "x "<< x<< " y "<< y  << std::endl;
        temp.push_back(std::make_pair(x, y));
    }
    if (include){
        std::string cmd = "with points title ''";
        plots.add_plot1d(temp, cmd);
        
       // xy_pts_A.push_back(temp);
    }
    
    
    *gp << "set xrange  [" << xmin<< ":" << xmax<<"]\n";
    *gp << "set yrange  [" << ymin<< ":" << ymax<<"]\n";
    *gp << "set arrow from "<< verticalLine1<< ","<< ymin <<" to "<< verticalLine1 << "," <<ymax <<" nohead lc rgb \'red\'\n";
    *gp << plots;
    
    return 0;
}
int GNUPlotter::plot2dSerie(std::vector<double> &xseries,  std::vector<double> &yseries,bool doSavePng, std::string fileName,std::string plotTitle ){
    
    int result = 0;
    if (yseries.size()==0)
        return result;
    
    if (doSavePng){
        *gp << "set terminal png size 900,900\n";
        *gp << "set output  '"<< fileName << "'\n";
    }
    else{
        *gp << "set term qt persist size 900,900\n";
        
    }
    
    *gp << "set title  '"<< plotTitle <<"'\n";
    *gp << "set xlabel 'time'"<<"\n";
    *gp <<  "set ylabel ''"<<"\n";
    
    
    //size_t nYSeries = yseries.size();
    std::vector<std::pair<double, double>> xy_pts_A;
    std::vector<std::pair<double, double>> temp;

    
    size_t nYSeries = yseries.size();
    size_t nXSeries = xseries.size();
    assert(nXSeries == nYSeries);
    
    double y;
    double x;
    
    bool include = true;

    double maxAxisValue = 10000;
     double xmax = -maxAxisValue;
     double xmin = maxAxisValue;
     double ymax = -maxAxisValue;
     double ymin = maxAxisValue;
    auto plots = (*gp).plotGroup();
    size_t size =nXSeries ;
    
    
    for(size_t i=0 ; i< size; ++i) {
        y= yseries.at(i);
        x= xseries.at(i);
        
        if (y< -maxAxisValue){
            include = false;
            break;
        }
        if (y > ymax)
            ymax= y;
        if (y < ymin)
            ymin= y;
        if (x > xmax)
            xmax= x;
        if (x < xmin)
            xmin= x;
        
        // std::cout << "x "<< x<< " y "<< y  << std::endl;
        temp.push_back(std::make_pair(x, y));
    }
    if (include){
        std::string cmd = "with points title ''";
        plots.add_plot1d(temp, cmd);
        
       // xy_pts_A.push_back(temp);
    }
    
    
    *gp << "set xrange  [" << xmin<< ":" << xmax<<"]\n";
    *gp << "set yrange  [" << ymin<< ":" << ymax<<"]\n";

    *gp << plots;
    
    return 0;
}
int GNUPlotter::plot2dSerieWithSegments(std::vector<double> &xseries,  std::vector<double> &yseries, std::vector<double> &xSegmentsUpper,std::vector<double> &ySegmentsUpper,std::vector<double> &xSegmentsLower,std::vector<double> &ySegmentsLower, bool doSavePng, std::string fileName,std::string plotTitle ){
    
    int result = 0;
    if (yseries.size()==0)
        return result;
    
    if (doSavePng){
        *gp << "set terminal png size 900,900\n";
        *gp << "set output  '"<< fileName << "'\n";
    }
    else{
        *gp << "set term qt persist size 900,900\n";
        
    }
    
    *gp << "set title  '"<< plotTitle <<"'\n";
    *gp << "set xlabel 'time'"<<"\n";
    *gp <<  "set ylabel ''"<<"\n";
    
    
    //size_t nYSeries = yseries.size();
    std::vector<std::pair<double, double>> xy_pts_A;
    std::vector<std::pair<double, double>> temp;

    
    size_t nYSeries = yseries.size();
    size_t nXSeries = xseries.size();
    assert(nXSeries == nYSeries);
    
    double y;
    double x;
    
    bool include = true;

    double maxAxisValue = 10000;
     double xmax = -maxAxisValue;
     double xmin = maxAxisValue;
     double ymax = -maxAxisValue;
     double ymin = maxAxisValue;
    auto plots = (*gp).plotGroup();
    size_t size =nXSeries ;
    
   // std::vector<boost::tuple<double, double, double, double> > segments_upper;
   // std::vector<boost::tuple<double, double, double, double> > segments_lower;
    std::vector<std::pair<double, double>> points_upper;
    std::vector<std::pair<double, double>> points_lower;
    
    assert(xSegmentsUpper.size()==ySegmentsUpper.size());
    for(size_t i=0 ; i< xSegmentsUpper.size()-1; ++i) {
        
       points_upper.emplace_back(xSegmentsUpper[i],ySegmentsUpper[i]);
//        double dx_upper = xSegmentsUpper[i+1]-xSegmentsUpper[i];
//        double dy_upper = ySegmentsUpper[i+1]-ySegmentsUpper[i];
//        segments_upper.push_back(boost::make_tuple(
//             xSegmentsUpper[i],
//             ySegmentsUpper[i],
//             dx_upper,
//             dy_upper
//        ));
    }
    
    assert(xSegmentsLower.size()==ySegmentsLower.size());
    for(size_t i=0 ; i< xSegmentsLower.size()-1; ++i) {
 
        points_lower.emplace_back(xSegmentsLower[i],ySegmentsLower[i]);
//         double dx_lower = xSegmentsLower[i+1]-xSegmentsLower[i];
//         double dy_lower = ySegmentsLower[i+1]-ySegmentsLower[i];
//         segments_upper.push_back(boost::make_tuple(
//              xSegmentsLower[i],
//              ySegmentsLower[i],
//              dx_lower,
//              dy_lower
//         ));
     }
    
    for(size_t i=0 ; i< size; ++i) {
        y= yseries.at(i);
        x= xseries.at(i);
        
        if (y< -maxAxisValue){
            include = false;
            break;
        }
        if (y > ymax)
            ymax= y;
        if (y < ymin)
            ymin= y;
        if (x > xmax)
            xmax= x;
        if (x < xmin)
            xmin= x;
        
        // std::cout << "x "<< x<< " y "<< y  << std::endl;
        temp.push_back(std::make_pair(x, y));
    }
    if (include){
        std::string cmd = "with points title ''";
        plots.add_plot1d(temp, cmd);
        
       // xy_pts_A.push_back(temp);
    }
    
    
    *gp << "set xrange  [" << xmin<< ":" << xmax<<"]\n";
    *gp << "set yrange  [" << ymin<< ":" << ymax<<"]\n";
    
    *gp << "set linetype 1 lw 1 lc rgb \'red\' \n";
    *gp << "set linetype 2 lw 2 lc rgb \'blue\' \n";
     *gp << "set linetype cycle 2 \n";
    
   // *gp << "plot '-' with lines  title 'upper bound', '-' with vectors nohead title 'lower bound'\n";
    //*gp << cmd;
     plots.add_plot1d(points_upper,"with lines  title 'upper bound' ");
    //plots.add_plot1d(segments_upper, cmd);
    //*gp->send1d(points_upper);
    //cmd = "with vectors nohead title 'lower bound'  lc rgb \'blue\'\n";
    //*gp << cmd;
    //plots.add_plot1d(segments_lower, cmd);
    plots.add_plot1d(points_lower,"with lines  title 'lower bound' ");
   // *gp->send1d(points_lower);

    *gp << plots;
    
    return 0;
}
