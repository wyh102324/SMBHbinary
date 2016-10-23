#include"astroSimulation.h"

inline double GetForce(const VEC &p1, const VEC &p2, const double &mass,const int index)
{
    static double r, r2;
    static double dx, dy, dz, n_component;
    static VEC n;
    dx = p1.x - p2.x , dy = p1.y - p2.y , dz= p1.z - p2.z;
    r = sqrt(dx*dx + dy*dy + dz*dz);
    r2 = r*r ;
    n.x = dx/r , n.y = dy/r ,n.z= dz/r;

    n_component = *((double*)&n + index-1);
    return -G*mass*n_component/r2;

}

inline double GetPNForce(const VEC &p1, const VEC &p2, const double &mass1，const double &mass2,const int index)
{
    static double r, r2, r3;
    static double dx, dy, dz, n_component;
    static VEC n;
    dx = p1.x - p2.x , dy = p1.y - p2.y , dz= p1.z - p2.z;
    r = sqrt(dx*dx + dy*dy + dz*dz);
    r2 = r*r ,r3 = r*r2;
    n.x = dx/r , n.y = dy/r ,n.z= dz/r;
    
    n_component = *((double*)&n + index-1);
    return (-G*mass2*/r2-c2Recip*6*G*G*mass2*(mass1+mass2)/r3) *n_component;
    
}

