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
     double xmax = -10000;
     double xmin = 100000;
     double ymax = -10000;
     double ymin = 100000;
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
             
             if (y< -10000000){
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

int GNUPlotter::plot2dSerie(double xfrom, double xto, int n, std::vector<double> &yseries, std::vector<std::string> &labels ){
    
    int result = 0;
    if (yseries.size()==0)
        return result;
   
   *gp << "set title  'Log Likelihood ratios for differents pairs'"<<"\n";
    *gp << "set xlabel 'Increments'"<<"\n";
   *gp <<  "set ylabel 'Log Likelihood Ratios'"<<"\n";
    *gp << "set xrange  [" << xfrom<< ":" << xto<<"]\n";
    
    //size_t nYSeries = yseries.size();
    std::vector<std::pair<double, double>> xy_pts_A;
    std::vector<std::pair<double, double>> temp;
    int i = 0;
    double y;
    
    double d = (xto-xfrom)/(1.0*n);
    bool include = true;
    i= 0;
    double ymax = -10000;
    double ymin = 100000;
    auto plots = (*gp).plotGroup();

    for(double x=xfrom ; x<= xto; x+=d) {
            y= yseries.at(i);
            
            if (y< -10000){
                include = false;
                continue;
            }
            if (y > ymax)
                ymax= y;
            if (y < ymin)
               ymin= y;
            i++;
           // std::cout << "x "<< x<< " y "<< y  << std::endl;
            xy_pts_A.push_back(std::make_pair(x, y));
        }
   
        std::string cmd = "with points title 'loglik ratio increment" + labels[i]+"'";             plots.add_plot1d(temp, cmd);
        
             // xy_pts_A.push_back(temp);
          
    
       *gp << "set yrange  [" << ymin<< ":" << ymax<<"]\n";
      *gp << plots;

    return 0;
}
