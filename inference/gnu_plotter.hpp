//
//  gnu_plotter.hpp
//  run
//
//  Created by Fausto Fabian Crespo Fernandez on 14/09/2021.
//

#ifndef gnu_plotter_hpp
#define gnu_plotter_hpp

#include "gnuplot-iostream.h"


class GNUPlotter{
    
public:
    
    Gnuplot *gp;
    
    GNUPlotter();
    
    void setToSavePng(std::string path);
    int plot2dSeveralSeries(std::vector<double> &x, std::vector<std::vector<double>> &yseries );
    int plot2dSeveralSeries(std::vector<std::vector<double>> &xseries, std::vector<std::vector<double>> &yseries,std::vector<std::string> &labels,  bool savePng, std::string filePath, double verticalLine1, double verticalLine2 );
    int plot2dSerie(std::vector<double> &xseries,  std::vector<double> &yseries, std::vector<std::string> &labels,bool doSavePng, std::string fileNameIncrements, double verticalLine1   );
    
    int plot2dSerie(std::vector<double> &xseries,  std::vector<double> &yseries,bool doSavePng, std::string fileName,std::string plotTitle );
    
    int plot2dSerieWithSegments(std::vector<double> &xseries,  std::vector<double> &yseries, std::vector<double> &xSegmentsUpper,std::vector<double> &ySegmentsUpper,std::vector<double> &xSegmentsLower,std::vector<double> &ySegmentsLower, bool doSavePng, std::string fileName,std::string plotTitle);
};

#endif /* gnu_plotter_hpp */
