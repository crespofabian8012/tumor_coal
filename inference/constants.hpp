/*
 * constants definition
 */

#ifndef constants_hpp
#define constants_hpp

const int NUM_COLS = 6;


/*                                        AA CC GG TT AC AG AT CG CT GT            */
static const char mut_dist[10][10] = {  {  0, 2, 2, 2, 1, 1, 1, 2, 2, 2 },   /* AA */
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

/*                                        AA CC GG TT  AC AG AT CG CT GT  CA GA TA GC TC TG           */
static const char mut_dist16[16][16] = {{  0, 2, 2, 2,  1, 1, 1, 2, 2, 2,  1, 1, 1, 2, 2, 2 },  /*AA*/
                                        {  2, 0, 2, 2,  1, 2, 2, 1, 1, 2,  1, 2, 2, 1, 1, 2 },  /*CC*/
                                        {  2, 2, 0, 2,  2, 1, 2, 1, 2, 1,  2, 1, 2, 1, 2, 1 },  /*GG*/
                                        {  2, 2, 2, 0,  2, 2, 1, 2, 1, 1,  2, 2, 1, 2, 1, 1 },  /*TT*/
                                        {  1, 1, 2, 2,  0, 1, 1, 2, 2, 2,  2, 2, 2, 1, 1, 2 },  /*AC*/
                                        {  1, 2, 1, 2,  1, 0, 1, 1, 2, 2,  2, 2, 2, 2, 2, 1 },  /*AG*/
                                        {  1, 2, 2, 1,  1, 1, 0, 2, 1, 1,  2, 2, 2, 2, 2, 2 },  /*AT*/
                                        {  2, 1, 1, 2,  2, 1, 2, 0, 1, 2,  1, 2, 2, 2, 2, 1 },  /*CG*/
                                        {  2, 1, 2, 1,  2, 2, 1, 1, 0, 1,  1, 2, 2, 2, 2, 2 },  /*CT*/
                                        {  2, 2, 1, 1,  2, 2, 1, 2, 1, 0,  2, 1, 2, 1, 2, 2 },  /*GT*/
                                        {  1, 1, 2, 2,  2, 2, 2, 1, 1, 2,  0, 1, 1, 2, 2, 2 },  /*CA*/
                                        {  1, 2, 1, 2,  2, 2, 2, 2, 2, 1,  1, 0, 1, 1, 2, 2 },  /*GA*/
                                        {  1, 2, 2, 1,  2, 2, 2, 2, 2, 2,  1, 1, 0, 2, 1, 1 },  /*TA*/
                                        {  2, 1, 1, 2,  1, 2, 2, 2, 2, 1,  2, 1, 2, 0, 1, 2 },  /*GC*/
                                        {  2, 1, 2, 1,  1, 2, 2, 2, 2, 2,  2, 2, 1, 1, 0, 1 },  /*TC*/
                                        {  2, 2, 1, 1,  2, 1, 2, 1, 2, 2,  2, 2, 1, 2, 1, 0 }   /*TG*/
                                        };


#endif /* constants_h */
