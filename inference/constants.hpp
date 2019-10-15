//
//  constants.h
//  tumor_coal
//
//  Created by Fausto Fabian Crespo Fernandez on 2019-10-08.
//

#ifndef constants_h
#define constants_h

const int NUM_COLS = 6;

/*                                        AA CC GG TT AC AG AT CG CT GT            */
static const char mut_dist[10][10] = {
    {  0, 2, 2, 2, 1, 1, 1, 2, 2, 2 },   /* AA */
    {  2, 0, 2, 2, 1, 2, 2, 1, 1, 2 },   /* CC */
    {  2, 2, 0, 2, 2, 1, 2, 1, 2, 1 },   /* GG */
    {  2, 2, 2, 0, 2, 2, 1, 2, 1, 1 },   /* TT */
    {  1, 1, 2, 2, 0, 1, 1, 1, 1, 2 },   /* AC */
    {  1, 2, 1, 2, 1, 0, 1, 1, 2, 1 },   /* AG */
    {  1, 2, 2, 1, 1, 1, 0, 2, 1, 1 },   /* AT */
    {  2, 1, 1, 2, 1, 1, 2, 0, 1, 1 },   /* CG */
    {  2, 1, 2, 1, 1, 2, 1, 1, 0, 1 },   /* CT */
    {  2, 2, 1, 1, 2, 1, 1, 1, 1, 0 }    /* GT */
};



#endif /* constants_h */
