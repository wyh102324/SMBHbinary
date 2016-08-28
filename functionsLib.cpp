#ifndef FUNCTIONSLIB
#define FUNCTIONSLIB

#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"astroSimulation.h"

double GetDistance(VEC& p1, VEC& p2)
{
    static double dx, dy, dz;
    dx = p1.x - p2.x, dy = p1.y - p2.y, dz = p1.z - p2.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

double GetMergeLimit(double m1, double m2)
{
    return (pow(m1/m_sun, 3.0/4) + pow(m2/m_sun, 3.0/4))*R_sun;
}

double GetTEDRadius(double alpha, double MC, double m)
{
 return alpha*pow(m/m_sun, 3.0/4)*R_sun*pow(MC/m, 1.0/3);
}

void MatterSet(double*m_des, double*m_src, int ParticleNumber)
{
    for(int i = 0;i < ParticleNumber;i++)
    {
        m_des[i] = m_src[i];
    }
}

void SetNatureParticleLabel(const IndexLabel &Label, const int ParticleNumber)
{
    for(int i = 0; i < ParticleNumber;i++)
        Label.re_seq[i] = i, Label.seq[i] = i;
}

void Integrator_Leapfrog(DynamicSystem &DS, const double StepLength, const DynamicSystem &DS_src)
{
     double factor = StepLength*StepLength*0.5, hdt = 0.5*StepLength;
     for(int i = 0;i < DS_src.ParticleNumber;i++)
     {
        DS.position[i].x = DS_src.position[i].x + (DS_src.velocity[i].x * StepLength + DS_src.acceleration[i].x * factor);
        DS.position[i].y = DS_src.position[i].y + (DS_src.velocity[i].y * StepLength + DS_src.acceleration[i].y * factor);
        DS.position[i].z = DS_src.position[i].z + (DS_src.velocity[i].z * StepLength + DS_src.acceleration[i].z * factor);
     }
     for(int i = 0;i < DS_src.ParticleNumber;i++)
     {
        DS.velocity[i].x = DS_src.velocity[i].x + (DS_src.acceleration[i].x * hdt);
        DS.velocity[i].y = DS_src.velocity[i].y + (DS_src.acceleration[i].y * hdt);
        DS.velocity[i].z = DS_src.velocity[i].z + (DS_src.acceleration[i].z * hdt);
     }
     for(int i = 0;i < DS_src.ParticleNumber;i++)
     {
         DS.acceleration[i].x = DS.acceleration[i].y = DS.acceleration[i].z = 0;
         for(int j = 0;j < DS_src.ParticleNumber;j++)
         {
             if(i!=j)
             {
                 DS.acceleration[i].x += GetForce(DS.position[i], DS.position[j], DS_src.Mass[j], 1);
                 DS.acceleration[i].y += GetForce(DS.position[i], DS.position[j], DS_src.Mass[j], 2);
                 DS.acceleration[i].z += GetForce(DS.position[i], DS.position[j], DS_src.Mass[j], 3);
             }
         }
     }
     for(int i = 0;i < DS_src.ParticleNumber;i++)
     {
        DS.velocity[i].x += (DS.acceleration[i].x * hdt);
        DS.velocity[i].y += (DS.acceleration[i].y * hdt);
        DS.velocity[i].z += (DS.acceleration[i].z * hdt);
     }
}

void StatusCopy(DynamicSystem &DS_dest, const DynamicSystem &DS_src)
{
memcpy(DS_dest.Mass, DS_src.Mass, sizeof(double)*DS_src.ParticleNumber);
memcpy(DS_dest.position, DS_src.position, sizeof(VEC)*DS_src.ParticleNumber);
memcpy(DS_dest.velocity, DS_src.velocity, sizeof(VEC)*DS_src.ParticleNumber);
memcpy(DS_dest.acceleration, DS_src.acceleration, sizeof(VEC)*DS_src.ParticleNumber);
memcpy(DS_dest.ParticleLabel.seq, DS_src.ParticleLabel.seq, sizeof(int)*DS_src.ParticleNumber);
memcpy(DS_dest.ParticleLabel.re_seq, DS_src.ParticleLabel.re_seq, sizeof(int)*DS_src.ParticleNumber);
DS_dest.ParticleNumber = DS_src.ParticleNumber;
}

void EvolutionTracking(FILE* EvolutionTrackFile, const double &time, const DynamicSystem &DS, const int Number)
{
    fprintf(EvolutionTrackFile, "%lf ", time);
    for(int i = 0 ; i < Number; i++)
    {
        fprintf(EvolutionTrackFile, "%lf %lf %lf %lf %lf %lf ", DS.position[DS.ParticleLabel.seq[i]].x/AU, DS.position[DS.ParticleLabel.seq[i]].y/AU, DS.position[DS.ParticleLabel.seq[i]].z/AU, DS.velocity[DS.ParticleLabel.seq[i]].x, DS.velocity[DS.ParticleLabel.seq[i]].y, DS.velocity[DS.ParticleLabel.seq[i]].z);
    }
    fprintf(EvolutionTrackFile, "\r\n");
}

void WriteToFileBit(FILE* DataRecord_file, const double &time, VEC* position, VEC* velocity, const int ParticleNumber)
{
    fwrite(&time, sizeof(double), 1, DataRecord_file);
    fwrite(position, sizeof(VEC), ParticleNumber, DataRecord_file);
    fwrite(velocity, sizeof(VEC), ParticleNumber, DataRecord_file);
}

void GetCentroid(const VEC &position1, const  VEC &position2, const double* M, VEC &centroidv)
{
    centroidv.x = (position1.x*M[0] + position2.x*M[1])/(M[0] + M[1]);
    centroidv.y = (position1.y*M[0] + position2.y*M[1])/(M[0] + M[1]);
    centroidv.z = (position1.z*M[0] + position2.z*M[1])/(M[0] + M[1]);
}

const double ERR = 1e-7;
bool errCheck_Orbits(const DynamicSystem &DS_try, const DynamicSystem &DS_htry, const int OrbitNumber)
{
    bool e = 0, te;
    double dr;
    double r1, tmp;
    VEC centroid_try, centroid_htry;
    if(OrbitNumber == 3)
    {
        dr = GetDistance(DS_try.position[0], DS_try.position[1]);
        r1 = GetDistance(DS_htry.position[0], DS_htry.position[1]);
        e =  fabs(dr - r1)/dr > ERR;

        dr = GetDistance(DS_try.position[0], DS_try.position[2]);
        r1 = GetDistance(DS_htry.position[0], DS_htry.position[2]);
        te = fabs(dr - r1)/dr > ERR;
        e = e||te;
    }
    else if(OrbitNumber == 4)
    {
        dr = GetDistance(DS_try.position[DS_try.ParticleLabel.seq[3]], DS_try.position[DS_try.ParticleLabel.seq[2]]);
        r1 = GetDistance(DS_htry.position[DS_htry.ParticleLabel.seq[3]], DS_htry.position[DS_htry.ParticleLabel.seq[2]]);
        e = fabs(dr - r1)/dr > ERR;

        dr = GetDistance(DS_try.position[DS_try.ParticleLabel.seq[0]], DS_try.position[DS_try.ParticleLabel.seq[1]]);
        r1 = GetDistance(DS_htry.position[DS_htry.ParticleLabel.seq[0]], DS_htry.position[DS_htry.ParticleLabel.seq[1]]);
        tmp = fabs(dr - r1)/dr;
        if(tmp!=0)
        {
            te =  tmp > ERR;
            e = e||te;
        }

        GetCentroid(DS_try.position[DS_try.ParticleLabel.seq[0]], DS_try.position[DS_try.ParticleLabel.seq[1]], DS_try.Mass, centroid_try);
        GetCentroid(DS_htry.position[DS_htry.ParticleLabel.seq[0]], DS_htry.position[DS_htry.ParticleLabel.seq[1]], DS_htry.Mass, centroid_htry);
        dr = GetDistance(centroid_try, DS_try.position[DS_try.ParticleLabel.seq[2]]);
        r1 = GetDistance(centroid_htry, DS_htry.position[DS_htry.ParticleLabel.seq[2]]);
        tmp = fabs(dr - r1)/dr;
        if(tmp!=0)
        {
            te =  tmp > ERR;
            e = e||te;
        }
    }
    else
    {
        e = 0;
    }
    return e;
}
#endif
