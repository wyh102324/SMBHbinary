#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"astroSimulation.h"
#include"functionsLib.cpp"

void DynamicSystemInitial(DynamicSystem &DS, double* M, const int IniParticleNumber)
{
    DS.ParticleNumber = IniParticleNumber;
    DS.position = new VEC[IniParticleNumber];
    DS.velocity = new VEC[IniParticleNumber];
    DS.acceleration = new VEC[IniParticleNumber];
    DS.Mass = new double[IniParticleNumber];
    MatterSet(DS.Mass, M, IniParticleNumber);
    DS.ParticleLabel.seq = new int[IniParticleNumber];
    DS.ParticleLabel.re_seq = new int[IniParticleNumber];
    SetNatureParticleLabel(DS.ParticleLabel, IniParticleNumber);
}

void DestroyDynamicSystem(DynamicSystem &DS)
{
    delete DS.Mass;
    delete DS.ParticleLabel.seq;
    delete DS.ParticleLabel.re_seq;
    delete DS.position;
    delete DS.velocity;
    delete DS.acceleration;
    return;
}

void ResetParticle(DynamicSystem &DS, const int TDE_Particle)//resent the label of particles after TDE merger.
{
    VEC swap_p, swap_v, swap_a;
    double swap_m;
    swap_p = DS.position[DS.ParticleLabel.seq[TDE_Particle]];
    swap_v = DS.velocity[DS.ParticleLabel.seq[TDE_Particle]];
    swap_a = DS.acceleration[DS.ParticleLabel.seq[TDE_Particle]];
    swap_m = DS.Mass[DS.ParticleLabel.seq[TDE_Particle]];
    int i = 0;
    for(i = DS.ParticleLabel.seq[TDE_Particle];i < DS.ParticleNumber - 1;i++)
    {
        DS.position[i] = DS.position[i+1];
        DS.velocity[i] = DS.velocity[i+1];
        DS.acceleration[i] = DS.acceleration[i+1];
        DS.Mass[i] = DS.Mass[i+1];
        DS.ParticleLabel.seq[ DS.ParticleLabel.re_seq[i+1] ] = i;
        DS.ParticleLabel.re_seq[i] = DS.ParticleLabel.re_seq[i+1];
    }
    DS.position[i] = swap_p;
    DS.velocity[i] = swap_v;
    DS.acceleration[i] = swap_a;
    DS.Mass[i] = swap_m;
    DS.ParticleLabel.re_seq[i] = TDE_Particle;
    DS.ParticleLabel.seq[TDE_Particle] = DS.ParticleNumber - 1;
    DS.ParticleNumber--;
}

void LastStatusKeep(FILE* DestinyInformatioFile, const DynamicSystem &DS, const int Number)
{
    for(int i = 0 ; i < Number; i++)
    {
        fprintf(DestinyInformatioFile, "%lf %lf %lf %lf %lf %lf ", DS.position[DS.ParticleLabel.seq[i]].x/AU, DS.position[DS.ParticleLabel.seq[i]].y/AU, DS.position[DS.ParticleLabel.seq[i]].z/AU, DS.velocity[DS.ParticleLabel.seq[i]].x, DS.velocity[DS.ParticleLabel.seq[i]].y, DS.velocity[DS.ParticleLabel.seq[i]].z);
    }
    fprintf(DestinyInformatioFile, "\r\n");
}

int JudgeDestiny(FILE* DestinyInformatioFile, const double &CurrentTime, DynamicSystem &DS, const int IniParticleNumber)
{
    static double DistanceMeasure;
    static double TEDR1 = GetTEDRadius(3, DS.Mass[2], DS.Mass[0]);
    static double TEDR2 = GetTEDRadius(3, DS.Mass[2], DS.Mass[1]);
    static double MergeLimit = GetMergeLimit(DS.Mass[0], DS.Mass[1]);
    static int destiny = 0;
    if(CurrentTime==0)
        destiny = UNEVENTFUL;
    if(destiny == UNEVENTFUL)
    {
        DistanceMeasure = GetDistance(DS.position[0], DS.position[1]);
        if(DistanceMeasure < MergeLimit)
        {
            fprintf(DestinyInformatioFile, "%d %lf %d %lf ", MERGER, CurrentTime, UNEVENTFUL, CurrentTime);
            LastStatusKeep(DestinyInformatioFile, DS, IniParticleNumber);
            return END_EVOLUTION;
        }
    }
    if(DS.ParticleLabel.seq[0] < DS.ParticleNumber)//this means m0 has not been TDE or merger(TDE, merger particles will be moved to the end of the label sequence)
    {
        DistanceMeasure = GetDistance(DS.position[DS.ParticleLabel.seq[0]], DS.position[DS.ParticleLabel.seq[2]]);
        if(DistanceMeasure < TEDR1)
        {
            if(CurrentTime == 0)
            {
                fprintf(DestinyInformatioFile, "%d %lf %d %lf ", INITIAL_INVALID, CurrentTime, INITIAL_INVALID, CurrentTime);
                LastStatusKeep(DestinyInformatioFile, DS, IniParticleNumber);
                return END_EVOLUTION;
            }
            else
            {
                destiny += M0_TDE;
                fprintf(DestinyInformatioFile, "%d %lf ", M0_TDE, CurrentTime);
                ResetParticle(DS, 0);
                if(DS.ParticleNumber <= IniParticleNumber - 2)
                {
                    LastStatusKeep(DestinyInformatioFile, DS, IniParticleNumber);
                    return END_EVOLUTION;
                }
            }
        }
    }
    if(DS.ParticleLabel.seq[1] < DS.ParticleNumber)//this means m1 has not been TDE or merger(TDE, merger particles will be moved to the end of the label sequence)
    {
        DistanceMeasure = GetDistance(DS.position[DS.ParticleLabel.seq[1]], DS.position[DS.ParticleLabel.seq[2]]);
        if(DistanceMeasure < TEDR2)
        {
            if(CurrentTime == 0)
            {
                fprintf(DestinyInformatioFile, "%d %lf %d %lf ", INITIAL_INVALID, CurrentTime, INITIAL_INVALID, CurrentTime);
                LastStatusKeep(DestinyInformatioFile, DS, IniParticleNumber);
                return END_EVOLUTION;
            }
            else
            {
                destiny += M1_TDE;
                fprintf(DestinyInformatioFile, "%d %lf ", M1_TDE, CurrentTime);
                ResetParticle(DS, 1);
                if(DS.ParticleNumber <= IniParticleNumber - 2)
                {
                    LastStatusKeep(DestinyInformatioFile, DS, IniParticleNumber);
                    return END_EVOLUTION;
                }
            }
        }
    }
    return KEEP_EVOLUTION;
}

void DestinyFinish(FILE* DestinyInformatioFile, const double &CurrentTime, DynamicSystem &DS, const int IniParticleNumber)
{
    if(DS.ParticleNumber < IniParticleNumber)
    {
        fprintf(DestinyInformatioFile, "%d %lf ", UNEVENTFUL, CurrentTime);
        LastStatusKeep(DestinyInformatioFile, DS, IniParticleNumber);
    }
    else
    {
        fprintf(DestinyInformatioFile, "%d %lf %d %lf ", UNEVENTFUL, CurrentTime, UNEVENTFUL, CurrentTime);
        LastStatusKeep(DestinyInformatioFile, DS, IniParticleNumber);
    }
}

void initParameter_4Orbits(DynamicSystem &DS, const double* Major_Semi_Axis, const double* E, const double* Tilt)
{
    double r = Major_Semi_Axis[0]*(1 + E[0]), R = Major_Semi_Axis[1]*(1 + E[1]), R2 = Major_Semi_Axis[2]*(1 + E[2]);
    double M_total = DS.Mass[0] + DS.Mass[1] + DS.Mass[2] + DS.Mass[3];
    double vin = sqrt((1 - E[0])/(1 + E[0])*G*(DS.Mass[0] + DS.Mass[1])/Major_Semi_Axis[0]);
    double vout = sqrt((1 - E[1])/(1 + E[1])*G*(DS.Mass[0] + DS.Mass[1] + DS.Mass[2])/Major_Semi_Axis[1]);
    double vout1 = sqrt((1 - E[2])/(1 + E[2])*G*(M_total)/Major_Semi_Axis[2]);
    
    DS.position[0].x = 0 ;
    DS.position[0].y = (DS.Mass[2]/(DS.Mass[0] + DS.Mass[1] + DS.Mass[2]))*R*cos(Tilt[1]) + r*DS.Mass[1]/(DS.Mass[0] + DS.Mass[1])*cos(Tilt[0]);
    DS.position[0].z = -(DS.Mass[2]/(DS.Mass[0] + DS.Mass[1] + DS.Mass[2]))*R*sin(Tilt[1]) + r*DS.Mass[1]/(DS.Mass[0] + DS.Mass[1])*sin(Tilt[0]);
    DS.position[1].x = 0 ;
    DS.position[1].y = (DS.Mass[2]/(DS.Mass[0] + DS.Mass[1] + DS.Mass[2]))*R*cos(Tilt[1]) - r*DS.Mass[0]/(DS.Mass[0] + DS.Mass[1])*cos(Tilt[0]);
    DS.position[1].z = -(DS.Mass[2]/(DS.Mass[0] + DS.Mass[1] + DS.Mass[2]))*R*sin(Tilt[1]) - r*DS.Mass[0]/(DS.Mass[0] + DS.Mass[1])*sin(Tilt[0]);
    DS.position[2].x = 0 ;
    DS.position[2].y = -((DS.Mass[0] + DS.Mass[1])/(DS.Mass[0] + DS.Mass[1] + DS.Mass[2]))*R*cos(Tilt[1]) ;
    DS.position[2].z = ((DS.Mass[0] + DS.Mass[1])/(DS.Mass[0] + DS.Mass[1] + DS.Mass[2]))*R*sin(Tilt[1]);
    DS.position[3].x = 0 ;
    DS.position[3].y = R2*cos(Tilt[2]);
    DS.position[3].z = R2*sin(Tilt[2]) ;
    
    DS.velocity[0].y = 0 , DS.velocity[0].x = -(DS.Mass[2] + DS.Mass[3])/M_total*vout + DS.Mass[1]/(DS.Mass[0] + DS.Mass[1])*vin - DS.Mass[3]/M_total*vout1, DS.velocity[0].z = 0;
    DS.velocity[1].y = 0 , DS.velocity[1].x = -(DS.Mass[2] + DS.Mass[3])/M_total*vout - DS.Mass[0]/(DS.Mass[0] + DS.Mass[1])*vin - DS.Mass[3]/M_total*vout1, DS.velocity[1].z = 0;
    DS.velocity[2].y = 0 , DS.velocity[2].x = ((DS.Mass[0] + DS.Mass[1])*vout - DS.Mass[3]*vout1)/M_total , DS.velocity[2].z = 0;
    DS.velocity[3].y = 0 , DS.velocity[3].x = ((DS.Mass[0] + DS.Mass[1])*vout + (DS.Mass[0] + DS.Mass[1] + DS.Mass[2])*vout1)/M_total, DS.velocity[3].z = 0;
    for(int i=0;i < 4;i++)
    {
        DS.acceleration[i].x = GetForce(DS.position[i], DS.position[(i+1)%4], DS.Mass[(i+1)%4], 1)\
        + GetForce(DS.position[i], DS.position[(i + 2)%4], DS.Mass[(i+2)%4], 1)\
        + GetForce(DS.position[i], DS.position[(i+3)%4], DS.Mass[(i+3)%4], 1);
        DS.acceleration[i].y = GetForce(DS.position[i], DS.position[(i+1)%4], DS.Mass[(i+1)%4], 2)\
        + GetForce(DS.position[i], DS.position[(i+2)%4], DS.Mass[(i+2)%4], 2)\
        + GetForce(DS.position[i], DS.position[(i+3)%4], DS.Mass[(i+3)%4], 2);
        DS.acceleration[i].z = GetForce(DS.position[i], DS.position[(i+1)%4], DS.Mass[(i+1)%4], 3)\
        + GetForce(DS.position[i], DS.position[(i+2)%4], DS.Mass[(i+2)%4], 3)\
        + GetForce(DS.position[i], DS.position[(i+3)%4], DS.Mass[(i+3)%4], 3);
    }
    
}

void initParameter_3Orbits(DynamicSystem &DS, const double* Major_Semi_Axis, const double* E, const double* Tilt)
{
    double r = Major_Semi_Axis[0]*(1 + E[0]), R = Major_Semi_Axis[1]*(1 + E[1]);
    DS.position[0].x = 0 ;
    DS.position[0].y = (DS.Mass[2]/(DS.Mass[0] + DS.Mass[1] + DS.Mass[2]))*R*cos(Tilt[1]) + r*DS.Mass[1]/(DS.Mass[0] + DS.Mass[1])*cos(Tilt[0]);
    DS.position[0].z = -(DS.Mass[2]/(DS.Mass[0] + DS.Mass[1] + DS.Mass[2]))*R*sin(Tilt[1]) + r*DS.Mass[1]/(DS.Mass[0] + DS.Mass[1])*sin(Tilt[0]);
    DS.position[1].x = 0 ;
    DS.position[1].y = (DS.Mass[2]/(DS.Mass[0] + DS.Mass[1] + DS.Mass[2]))*R*cos(Tilt[1]) - r*DS.Mass[0]/(DS.Mass[0] + DS.Mass[1])*cos(Tilt[0]);
    DS.position[1].z = -(DS.Mass[2]/(DS.Mass[0] + DS.Mass[1] + DS.Mass[2]))*R*sin(Tilt[1]) - DS.Mass[0]/(DS.Mass[0] + DS.Mass[1])*r*sin(Tilt[0]);
    DS.position[2].x = 0 ;
    DS.position[2].y = -((DS.Mass[0] + DS.Mass[1])/(DS.Mass[0] + DS.Mass[1] + DS.Mass[2]))*R*cos(Tilt[1]);
    DS.position[2].z = ((DS.Mass[0] + DS.Mass[1])/(DS.Mass[0] + DS.Mass[1] + DS.Mass[2]))*R*sin(Tilt[1]);
    
    DS.velocity[0].y = 0;
    DS.velocity[0].x = -(DS.Mass[2]/(DS.Mass[0] + DS.Mass[1] + DS.Mass[2]))*sqrt((1 - E[1])/(1 + E[1])*G*(DS.Mass[2] + DS.Mass[1] + DS.Mass[0])/Major_Semi_Axis[1]) - DS.Mass[1]/(DS.Mass[0] + DS.Mass[1])*sqrt((1 - E[0])/(1 + E[0])*G*(DS.Mass[0] + DS.Mass[1])/Major_Semi_Axis[0]);
    DS.velocity[0].z = 0;
    DS.velocity[1].y = 0;
    DS.velocity[1].x = -(DS.Mass[2]/(DS.Mass[0] + DS.Mass[1] + DS.Mass[2]))*sqrt((1 - E[1])/(1 + E[1])*G*(DS.Mass[2] + DS.Mass[1] + DS.Mass[0])/Major_Semi_Axis[1]) + DS.Mass[0]/(DS.Mass[0] + DS.Mass[1])*sqrt((1 - E[0])/(1 + E[0])*G*(DS.Mass[0] + DS.Mass[1])/Major_Semi_Axis[0]);
    DS.velocity[1].z = 0;
    DS.velocity[2].x = ((DS.Mass[0] + DS.Mass[1])/(DS.Mass[0] + DS.Mass[1] + DS.Mass[2]))*sqrt((1 - E[1])/(1 + E[1])*G*(DS.Mass[2] + DS.Mass[1] + DS.Mass[0])/Major_Semi_Axis[1]);
    DS.velocity[2].y = 0;
    DS.velocity[2].z = 0;
    for(int i=0;i < 3;i++)
    {
        DS.acceleration[i].x = GetForce(DS.position[i], DS.position[(i+1)%3], DS.Mass[(i+1)%3], 1)\
        + GetForce(DS.position[i], DS.position[(i+2)%3], DS.Mass[(i+2)%3], 1);
        DS.acceleration[i].y = GetForce(DS.position[i], DS.position[(i+1)%3], DS.Mass[(i+1)%3], 2)\
        + GetForce(DS.position[i], DS.position[(i+2)%3], DS.Mass[(i+2)%3], 2);
        DS.acceleration[i].z = GetForce(DS.position[i], DS.position[(i+1)%3], DS.Mass[(i+1)%3], 3)\
        + GetForce(DS.position[i], DS.position[(i+2)%3], DS.Mass[(i+2)%3], 3);
    }
    
}

