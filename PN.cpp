#include"astroSimulation.h"
#include"math.h"
#include "stdafx.h"
inline double GetForce(const VEC& p1,const VEC &p2,const VEC &v1,const VEC &v2,const double &m1,const double &m2,const int index)
{
    static double r,r2,r3,r4,r5,r6;
    static double dx,dy,dz,v1s,v2s,v12s,cv,nv1,nv2,nv12,v_component,n_component;
    static VEC n,v12;
    dx = p1.x - p2.x , dy = p1.y - p2.y ;
    r = sqrt(dx*dx + dy*dy);
    r2 = r*r , r3=r2*r , r4=r2*r2 , r5=r2*r3 , r6=r3*r3;
    n.x = dx/r , n.y = dy/r ;
    v12.x = v1.x - v2.x ,v12.y = v1.y - v2.y;
    v1s = v1.x*v1.x + v1.y*v1.y;
    v2s = v2.x*v2.x + v2.y*v2.y;
    v12s = v12.x*v12.x + v12.y*v12.y;
    cv = v1.x*v2.x + v1.y*v2.y;
    nv1 = n.x*v1.x + n.y*v1.y;
    nv2 = n.x*v2.x + n.y*v2.y;
    nv12 = n.x*v12.x + n.y*v12.y;
    v_component = *((double*)&v12+index-1);
    n_component = *((double*)&n+index-1);
    return -G*m2*n_component/r2\
    +c2Recip*(\
              (5*G2*m1*m2/r3 + 4*G2*m2*m2/r3 + G*m2/r2*(1.5*nv2*nv2 - v1s + 4*cv -2*v2s) )*n_component \
               + G*m2/r2*(4*nv1-3*nv2)*v_component\
              );/*\
    +c4Recip*(\
              (\
               -14.25*G3*m1*m1*m2/r4 - 34.5*G3*m1*m2*m2/r4 - 9*G3*m2*m2*m2/r4\
               + G*m2/r2*(-1.875*nv2*nv2*nv2*nv2 + 1.5*nv2*nv2*v1s - 6*nv2*nv2*cv - 2*cv*cv + 4.5* nv2*nv2*v2s + 4*cv*v2s - 2*v2s*v2s)\
               + G2*m1*m2/r3*(19.5*nv1*nv1 - 39*nv1*nv2+8.5*nv2*nv2 - 3.75*v1s - 2.5*cv + 1.25*v2s)\
               + G2*m2*m2/r3*(2*nv1*nv1 - 4*nv1*nv2 - 6*nv2*nv2 - 8*cv + 4*v2s )\
              )*n_component \
             +(\
               G2*m2*m2/r3*(-2*nv1-2*nv2) + G2*m1*m2/r3*(-15.75*nv1+13.75*nv2)\
               +G*m2/r2*(-6*nv1*nv2*nv2 + 4.5*nv2*nv2*nv2 + nv2*v1s - 4*nv1*cv + 4*nv2*cv + 4*nv1*v2s - 5*nv2*v2s)\
              )*v_component \
             )\
    +c5Recip*(\
              (\
               208.0*G3*m1*m2*m2*nv12/15/r4 - 4.8* G3*m1*m1*m2/r4*nv12 + 2.4*G2*m1*m2/r3*nv12*v12s\
              )*n_component\
             +(\
               1.6*G3*m1*m1*m2/r4 - 6.4*G3*m1*m2*m2/r4 - 0.8*G2*m1*m2/r3*v12s\
              )*v_component\
             )\
    +c6Recip*(\
              (\
               G*m2/r2*(2.1785*nv2*nv2*nv2*nv2*nv2*nv2 - 1.875*nv2*nv2*nv2*nv2*v1s + 7.5*nv2*nv2*nv2*nv2*cv + 3*nv2*nv2*cv*cv\
                        -7.5*nv2*nv2*nv2*nv2*v2s + 1.5*nv2*nv2*v1s*v2s - 12*nv2*nv2*cv*v2s - 2*cv*cv*v2s + 7.5*nv2*nv2*v2s*v2s\
                        +4*cv*v2s*v2s - 2*v2s*v2s*v2s\
                        )\
               +G2*m1*m2/r3*(\
                             -21.375*nv1*nv1*nv1*nv1 + 85.5*nv1*nv1*nv1*nv2 - 180.75*nv1*nv1*nv2*nv2 + 191.5*nv1*nv2*nv2*nv2\
                             -56.875*nv2*nv2*nv2*nv2 + 57.25*nv1*nv1*v1s -102.5*nv1*nv2*v1s + 47.75*nv2*nv2*v1s - 11.375*v1s*v1s\
                             -114.5*nv1*nv1*cv + 244*nv1*nv2*cv - 112.5*nv2*nv2*cv + 45.5*v1s*cv - 44.25*cv*cv + 57.25*nv1*nv1*v2s\
                             -141.5*nv1*nv2*v2s + 64.75*nv2*nv2*v2s - 22.75*v1s*v2s + 43*cv*v2s - 10.125*v2s*v2s\
                            )\
               +G2*m2*m2/r3*(\
                             -6*nv1*nv1*nv2*nv2 + 12*nv1*nv2*nv2*nv2 + 6*nv2*nv2*nv2*nv2 + 4*nv1*nv2*cv + 12*nv2*nv2*cv\
                             +4*cv*cv - 4*nv1*nv2*v2s - 12*nv2*nv2*v2s - 8*cv*v2s + 4*v2s\
                            )\
               +G3*m2*m2*m2/r4*(-nv1*nv1 + 2*nv1*nv2 + 21.5*nv2*nv2 + 18*cv - 9*v2s )\
               +G3*m1*m2*m2/r4*(51.875*nv1*nv1 - 93.75*nv1*nv2 + 139.125*nv2*nv2 - 9.609375*nv12*PI*PI + 18*v1s +1.921875*PI*PI*v12s\
                                +33*cv - 16.5*v2s)\
               +G3*m1*m1*m2/r4*(\
                                -273.1369048*nv1*nv1 + 572.0238095*nv1*nv2 - 249.2619048*nv2*nv2 + 57.37738095*v1s - 86.2547619*cv\
                                +43.12738095*v2s + 110*nv12*nv12*log(r/g_r1)-22*v12s*log(r/g_r1)\
                               )\
               +16*G4*m2*m2*m2*m2/r5 + 149.7091387*G4*m1*m1*m2*m2/r5 + G4*m1*m1*m1*m2/r5*(-2.529365079+44.0/3*log(r/g_r1)) + G4*m1*m2*m2*m2/r5*(150.4885038 - 44.0/3*log(r/g_r2))\
              )*n_component\
             +(\
               G*m2/r2*(\
                        7.5*nv1*nv2*nv2*nv2*nv2 - 5.625*nv2*nv2*nv2*nv2*nv2 - 1.5*nv2*nv2*nv2*v1s + 6*nv1*nv2*nv2*cv - 6*nv2*nv2*nv2*cv\
                        -2*nv2*cv*cv - 12*nv1*nv2*nv2*v2s + 12*nv2*nv2*nv2*v2s + nv2*v1s*v2s - 4*nv1*cv*v2s + 8*nv2*cv*v2s + 4*nv1*v2s*v2s -7*nv2*v2s*v2s\
                       )\
               +G2*m2*m2/r3*(\
                             -2*nv1*nv1*nv2 + 8*nv1*nv2*nv2 + 2*nv2*nv2*nv2 + 2*nv1*cv + 4*nv2*cv - 2*nv1*v2s - 4* nv2*v2s\
                            )\
               +G2*m1*m2/r3*(\
                             -60.75*nv1*nv1*nv1 + 141.25*nv1*nv1*nv2 - 67.25*nv1*nv2*nv2 - 95.0/12*nv2*nv2*nv2 + 25.875*nv1*v1s - 17.125*nv2*v1s -36*nv1*cv\
                             +6.75*nv2*cv + 10.125*nv1*v2s + 10.375*nv2*v2s\
                            )\
               +G3*m2*m2*m2/r4*(4*nv1+5*nv2)\
               +G3*m1*m2*m2/r4*(-38.375*nv1 + 59.875*nv2 + 3.84375*nv12*PI*PI)\
               +G3*m1*m1*m2/r4*(31397.0/420*nv1 - 36227.0/420*nv2 - 44*nv12*log(r/g_r1))\
              )*v_component\
             )\
    +c7Recip*(\
              (\
              G4*m1*m1*m1*m2/r5*(3992.0/105*nv1 - 4328.0/105*nv2)\
              +G4*m1*m1*m2*m2/r6*(-13576.0/105*nv1 + 2872.0/21*nv2)\
              -3172.0/21*G4*m1*m2*m2*m2/r6*nv12\
              +G3*m1*m1*m2/r4*(\
                               48*nv1*nv1*nv1 - 139*nv1*nv1*nv2 + 148.8*nv1*nv2*nv2 - 57.6*nv2*nv2*nv2 - 4888.0/105*nv1*v1s\
                               +5056.0/105*nv2*v1s + 2056.0/21*nv1*cv - 2224.0/21*nv2*cv - 1028.0/21*nv1*v2s + 5812.0/105*nv2*v2s\
                              )\
              +G3*m1*m2*m2/r4*(\
                               -116.4*nv1*nv1*nv1 + 349.2*nv1*nv1*nv2 - 390.8*nv1*nv2*nv2 + 158*nv2*nv2*nv2 + 3568.0/105*nv12*v1s\
                               -2864.0/35*nv1*cv + 10048.0/105*nv2*cv + 1432.0/35*nv1*v2s - 5752.0/105*nv2*v2s\
                              )\
              +G2*m1*m2/r3*(\
                            -56*nv12*nv12*nv12*nv12*nv12 + 60*nv1*nv1*nv1*v12s - 180*nv1*nv1*nv2*v12s + 174*nv1*nv2*nv2*v12s\
                            -54*nv2*nv2*nv2*v12s - 246.0/35*nv12*v1s*v1s + 1068.0/35*nv1*v1s*cv - 984.0/35*nv2*v1s*cv - 1068.0/35*nv1*cv*cv\
                            +180.0/7*nv2*cv*cv - 534.0/35*nv1*v1s*v2s + 90.0/7*nv2*v1s*v2s + 984.0/35*nv1*cv*v2s - 732.0/35*nv2*cv*v2s\
                            -204.0/35*nv1*v2s*v2s + 24.0/7*nv2*v2s*v2s\
                           )\
              )*n_component\
             +(\
                -184.0/21*G4*m1*m1*m1*m2/r5 + 6224.0/105*G4*m1*m1*m2*m2/r6 + 6388.0/105*G4*m1*m2*m2*m2/r6\
                +G3*m1*m1*m2/r4*(\
                                 52.0/15*nv1*nv1 - 56.0/15*nv1*nv2 - 44.0/15*nv2*nv2 - 132.0/35*v1s + 152.0/35*cv - 48.0/35*v2s\
                                )\
                +G3*m1*m2*m2/r4*(\
                                 454.0/15*nv1*nv1 - 74.4*nv1*nv2 + 854.0/15*nv2*nv2 - 152.0/21*v1s + 2864.0/105*cv - 1768.0/105*v2s\
                                )\
                +G2*m1*m2/r3*(\
                              60*nv12*nv12*nv12*nv12 - 69.6*nv1*nv1*v12s + 136.8*nv1*nv2*v12s - 66*nv2*nv2*v12s + 334.0/35*v1s*v1s\
                              -1336.0/35*v1s*cv + 1308.0/35*cv*cv + 654.0/35*v1s*v2s - 1252.0/35*cv*v2s + 292.0/35*v2s*v2s\
                             )\
               )*v_component\
             );*/
}
