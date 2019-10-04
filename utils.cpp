//
//  utils.cpp
//  tumor_coal
//
//  Created by Fausto Fabian Crespo Fernandez on 5/10/19.
//

#include "utils.hpp"
int compare (const void * a, const void * b)

{
    
    double *p1 = (double *)a;
    
    double *p2 = (double *)b;
    
    
    
    if (*p1 > *p2) return 1;
    
    else if (*p2 > *p1 )return -1;
    
    else return 0;
    
}
