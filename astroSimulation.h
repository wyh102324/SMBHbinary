#ifndef ASTROSIMULATION_H_INCLUDED
#define ASTROSIMULATION_H_INCLUDED

#define G (6.67259e-11)
#define c (299792458)
#define m_sun (1.989e30)
#define R_sun (0.005*AU)
#define AU (149597870700.0)
#define PC (206256.0*AU)
#define UNEVENTFUL 0
#define MERGER -1
#define M0_TDE 1
#define M1_TDE 2
#define INITIAL_INVALID -2
#define END_EVOLUTION 1
#define KEEP_EVOLUTION 0
#define EVOLUTION_TRACK_ON 1
#define EVOLUTION_TRACK_OFF 0

const double G2 = G*G, G3 = G*G*G, G4 = G*G*G*G;
const double cRecip = 1.0/c;
const double c2Recip = cRecip*cRecip;
const double c3Recip = c2Recip*cRecip;
const double c4Recip = c2Recip*c2Recip;
const double c5Recip = c2Recip*c3Recip;
const double c6Recip = c3Recip*c3Recip;
const double c7Recip = c3Recip*c4Recip;
const double PI = acos(-1.0);

struct VEC
{
    double x, y, z;
};

struct IndexLabel
{
    int* re_seq;
    int* seq;
};

struct DynamicSystem
{
    IndexLabel ParticleLabel;
    int ParticleNumber;
    VEC*position;
    VEC*velocity;
    VEC*acceleration;
    double*Mass;

};

#endif // ASTROSIMULATION_H_INCLUDED
