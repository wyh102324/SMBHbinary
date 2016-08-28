#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"astroSimulation.h"

double GetDistance(VEC& p1, VEC& p2)
{
    static double dx,dy,dz;
    dx = p1.x - p2.x , dy = p1.y - p2.y , dz = p1.z - p2.z;
    return sqrt(dx*dx+dy*dy+dz*dz);
}

double GetMergeLimit(double m1,double m2)
{
    return (pow(m1/m_sun,3.0/4)+pow(m2/m_sun,3.0/4))*R_sun;
}

double GetTEDRadius(double alpha,double MC,double m)
{
 return alpha*pow(m/m_sun,3.0/4)*R_sun*pow(MC/m,1.0/3);
}

void MatterSet(double*m_des,double*m_src,int ParticleNumber)
{
    for(int i = 0;i<ParticleNumber;i++)
    {
        m_des[i]=m_src[i];
    }
}

void SetNatureParticleLabel(const IndexLabel& Label,const int ParticleNumber)
{
    for(int i=0; i<ParticleNumber;i++)
        Label.re_seq[i]=i,Label.seq[i]=i;
}

void Integrator_Leapfrog(DynamicSystem&DS,const double StepLength,const DynamicSystem&DS_src)
{
     double factor=StepLength*StepLength*0.5,hdt=0.5*StepLength;
     for(int i=0;i<DS_src.ParticleNumber;i++)
     {
        DS.position[i].x = DS_src.position[i].x + (DS_src.velocity[i].x * StepLength + DS_src.acceleration[i].x * factor);
        DS.position[i].y = DS_src.position[i].y + (DS_src.velocity[i].y * StepLength + DS_src.acceleration[i].y * factor);
        DS.position[i].z = DS_src.position[i].z + (DS_src.velocity[i].z * StepLength + DS_src.acceleration[i].z * factor);
     }
     for(int i=0;i<DS_src.ParticleNumber;i++)
     {
        DS.velocity[i].x = DS_src.velocity[i].x + (DS_src.acceleration[i].x * hdt);
        DS.velocity[i].y = DS_src.velocity[i].y + (DS_src.acceleration[i].y * hdt);
        DS.velocity[i].z = DS_src.velocity[i].z + (DS_src.acceleration[i].z * hdt);
     }
     for(int i=0;i<DS_src.ParticleNumber;i++)
     {
         DS.acceleration[i].x=DS.acceleration[i].y=DS.acceleration[i].z=0;
         for(int j=0;j<DS_src.ParticleNumber;j++)
         {
             if(i!=j)
             {
                 DS.acceleration[i].x += GetForce(DS.position[i],DS.position[j],DS_src.Mass[j],1);
                 DS.acceleration[i].y += GetForce(DS.position[i],DS.position[j],DS_src.Mass[j],2);
                 DS.acceleration[i].z += GetForce(DS.position[i],DS.position[j],DS_src.Mass[j],3);
             }
         }
     }
     for(int i=0;i<DS_src.ParticleNumber;i++)
     {
        DS.velocity[i].x += (DS.acceleration[i].x * hdt);
        DS.velocity[i].y += (DS.acceleration[i].y * hdt);
        DS.velocity[i].z += (DS.acceleration[i].z * hdt);
     }
}

void initParameter_4Orbits(DynamicSystem&DS,const double *Major_Semi_Axis,const double* E,const double*Tilt)
{
    double r = Major_Semi_Axis[0]*(1+E[0]),R=Major_Semi_Axis[1]*(1+E[1]),R2=Major_Semi_Axis[2]*(1+E[2]);
    double M_total=DS.Mass[0]+DS.Mass[1]+DS.Mass[2]+DS.Mass[3];
    double vin=sqrt((1-E[0])/(1+E[0])*G*(DS.Mass[0]+DS.Mass[1])/Major_Semi_Axis[0]);
    double vout=sqrt((1-E[1])/(1+E[1])*G*(DS.Mass[0]+DS.Mass[1]+DS.Mass[2])/Major_Semi_Axis[1]);
    double vout1=sqrt((1-E[2])/(1+E[2])*G*(M_total)/Major_Semi_Axis[2]);

    DS.position[0].x = 0 ;
    DS.position[0].y = (DS.Mass[2]/(DS.Mass[0]+DS.Mass[1]+DS.Mass[2]))*R*cos(Tilt[1])+r*DS.Mass[1]/(DS.Mass[0]+DS.Mass[1])*cos(Tilt[0]);
    DS.position[0].z = -(DS.Mass[2]/(DS.Mass[0]+DS.Mass[1]+DS.Mass[2]))*R*sin(Tilt[1])+r*DS.Mass[1]/(DS.Mass[0]+DS.Mass[1])*sin(Tilt[0]);
    DS.position[1].x = 0 ;
    DS.position[1].y = (DS.Mass[2]/(DS.Mass[0]+DS.Mass[1]+DS.Mass[2]))*R*cos(Tilt[1])-r*DS.Mass[0]/(DS.Mass[0]+DS.Mass[1])*cos(Tilt[0]);
    DS.position[1].z = -(DS.Mass[2]/(DS.Mass[0]+DS.Mass[1]+DS.Mass[2]))*R*sin(Tilt[1])-DS.Mass[0]/(DS.Mass[0]+DS.Mass[1])*r*sin(Tilt[0]);
    DS.position[2].x = 0 ;
    DS.position[2].y = -((DS.Mass[0]+DS.Mass[1])/(DS.Mass[0]+DS.Mass[1]+DS.Mass[2]))*R*cos(Tilt[1]) ;
    DS.position[2].z = ((DS.Mass[0]+DS.Mass[1])/(DS.Mass[0]+DS.Mass[1]+DS.Mass[2]))*R*sin(Tilt[1]);
    DS.position[3].x = 0 ;
    DS.position[3].y = R2*cos(Tilt[2]);
    DS.position[3].z = R2*sin(Tilt[2]) ;

    DS.velocity[0].y = 0 , DS.velocity[0].x = -(DS.Mass[2]+DS.Mass[3])/M_total*vout+DS.Mass[1]/(DS.Mass[0]+DS.Mass[1])*vin-DS.Mass[3]/M_total*vout1,DS.velocity[0].z = 0;
    DS.velocity[1].y = 0 , DS.velocity[1].x = -(DS.Mass[2]+DS.Mass[3])/M_total*vout-DS.Mass[0]/(DS.Mass[0]+DS.Mass[1])*vin-DS.Mass[3]/M_total*vout1,DS.velocity[1].z = 0;
    DS.velocity[2].y = 0 , DS.velocity[2].x = ((DS.Mass[0]+DS.Mass[1])*vout-DS.Mass[3]*vout1)/M_total , DS.velocity[2].z = 0;
    DS.velocity[3].y = 0 , DS.velocity[3].x = ((DS.Mass[0]+DS.Mass[1])*vout+(DS.Mass[0]+DS.Mass[1]+DS.Mass[2])*vout1)/M_total, DS.velocity[3].z = 0;
    for(int i=0;i<4;i++)
     {
        DS.acceleration[i].x = GetForce(DS.position[i],DS.position[(i+1)%4],DS.Mass[(i+1)%4],1)\
                             + GetForce(DS.position[i],DS.position[(i+2)%4],DS.Mass[(i+2)%4],1)\
                             + GetForce(DS.position[i],DS.position[(i+3)%4],DS.Mass[(i+3)%4],1);
        DS.acceleration[i].y = GetForce(DS.position[i],DS.position[(i+1)%4],DS.Mass[(i+1)%4],2)\
                             + GetForce(DS.position[i],DS.position[(i+2)%4],DS.Mass[(i+2)%4],2)\
                             + GetForce(DS.position[i],DS.position[(i+3)%4],DS.Mass[(i+3)%4],2);
        DS.acceleration[i].z = GetForce(DS.position[i],DS.position[(i+1)%4],DS.Mass[(i+1)%4],3)\
                             + GetForce(DS.position[i],DS.position[(i+2)%4],DS.Mass[(i+2)%4],3)\
                             + GetForce(DS.position[i],DS.position[(i+3)%4],DS.Mass[(i+3)%4],3);
     }

}

void initParameter_3Orbits(DynamicSystem&DS,const double *Major_Semi_Axis,const double* E,const double*Tilt)
{
    double r = Major_Semi_Axis[0]*(1+E[0]),R=Major_Semi_Axis[1]*(1+E[1]);
    DS.position[0].x = 0 ;
    DS.position[0].y = (DS.Mass[2]/(DS.Mass[0]+DS.Mass[1]+DS.Mass[2]))*R*cos(Tilt[1])+r*DS.Mass[1]/(DS.Mass[0]+DS.Mass[1])*cos(Tilt[0]);
    DS.position[0].z = -(DS.Mass[2]/(DS.Mass[0]+DS.Mass[1]+DS.Mass[2]))*R*sin(Tilt[1])+r*DS.Mass[1]/(DS.Mass[0]+DS.Mass[1])*sin(Tilt[0]);
    DS.position[1].x = 0 ;
    DS.position[1].y = (DS.Mass[2]/(DS.Mass[0]+DS.Mass[1]+DS.Mass[2]))*R*cos(Tilt[1])-r*DS.Mass[0]/(DS.Mass[0]+DS.Mass[1])*cos(Tilt[0]);
    DS.position[1].z = -(DS.Mass[2]/(DS.Mass[0]+DS.Mass[1]+DS.Mass[2]))*R*sin(Tilt[1])-DS.Mass[0]/(DS.Mass[0]+DS.Mass[1])*r*sin(Tilt[0]);
    DS.position[2].x = 0 ;
    DS.position[2].y = -((DS.Mass[0]+DS.Mass[1])/(DS.Mass[0]+DS.Mass[1]+DS.Mass[2]))*R*cos(Tilt[1]);
    DS.position[2].z = ((DS.Mass[0]+DS.Mass[1])/(DS.Mass[0]+DS.Mass[1]+DS.Mass[2]))*R*sin(Tilt[1]);

    DS.velocity[0].y = 0;
    DS.velocity[0].x = -(DS.Mass[2]/(DS.Mass[0]+DS.Mass[1]+DS.Mass[2]))*sqrt((1-E[1])/(1+E[1])*G*(DS.Mass[2]+DS.Mass[1]+DS.Mass[0])/Major_Semi_Axis[1])-DS.Mass[1]/(DS.Mass[0]+DS.Mass[1])*sqrt((1-E[0])/(1+E[0])*G*(DS.Mass[0]+DS.Mass[1])/Major_Semi_Axis[0]);
    DS.velocity[0].z = 0;
    DS.velocity[1].y = 0;
    DS.velocity[1].x = -(DS.Mass[2]/(DS.Mass[0]+DS.Mass[1]+DS.Mass[2]))*sqrt((1-E[1])/(1+E[1])*G*(DS.Mass[2]+DS.Mass[1]+DS.Mass[0])/Major_Semi_Axis[1])+DS.Mass[0]/(DS.Mass[0]+DS.Mass[1])*sqrt((1-E[0])/(1+E[0])*G*(DS.Mass[0]+DS.Mass[1])/Major_Semi_Axis[0]);
    DS.velocity[1].z = 0;
    DS.velocity[2].x = ((DS.Mass[0]+DS.Mass[1])/(DS.Mass[0]+DS.Mass[1]+DS.Mass[2]))*sqrt((1-E[1])/(1+E[1])*G*(DS.Mass[2]+DS.Mass[1]+DS.Mass[0])/Major_Semi_Axis[1]);
    DS.velocity[2].y = 0;
    DS.velocity[2].z = 0;
    for(int i=0;i<3;i++)
     {
        DS.acceleration[i].x = GetForce(DS.position[i],DS.position[(i+1)%3],DS.Mass[(i+1)%3],1)\
                             + GetForce(DS.position[i],DS.position[(i+2)%3],DS.Mass[(i+2)%3],1);
        DS.acceleration[i].y = GetForce(DS.position[i],DS.position[(i+1)%3],DS.Mass[(i+1)%3],2)\
                             + GetForce(DS.position[i],DS.position[(i+2)%3],DS.Mass[(i+2)%3],2);
        DS.acceleration[i].z = GetForce(DS.position[i],DS.position[(i+1)%3],DS.Mass[(i+1)%3],3)\
                             + GetForce(DS.position[i],DS.position[(i+2)%3],DS.Mass[(i+2)%3],3);
     }

}
void ResetParticle(DynamicSystem&DS,const int TDE_Particle)
{
    VEC swap_p,swap_v,swap_a;
    double swap_m;
    swap_p=DS.position[DS.ParticleLabel.seq[TDE_Particle]];
    swap_v=DS.velocity[DS.ParticleLabel.seq[TDE_Particle]];
    swap_a=DS.acceleration[DS.ParticleLabel.seq[TDE_Particle]];
    swap_m=DS.Mass[DS.ParticleLabel.seq[TDE_Particle]];
    int i=0;
    for(i=DS.ParticleLabel.seq[TDE_Particle];i<DS.ParticleNumber-1;i++)
    {
        DS.position[i]=DS.position[i+1];
        DS.velocity[i]=DS.velocity[i+1];
        DS.acceleration[i]=DS.acceleration[i+1];
        DS.Mass[i]=DS.Mass[i+1];
        DS.ParticleLabel.seq[ DS.ParticleLabel.re_seq[i+1] ]=i;
        DS.ParticleLabel.re_seq[i]=DS.ParticleLabel.re_seq[i+1];
    }
    DS.position[i]=swap_p;
    DS.velocity[i]=swap_v;
    DS.acceleration[i]=swap_a;
    DS.Mass[i]=swap_m;
    DS.ParticleLabel.re_seq[i]=TDE_Particle;
    DS.ParticleLabel.seq[TDE_Particle]=DS.ParticleNumber-1;
    DS.ParticleNumber--;
}

void StatusCopy(DynamicSystem&DS_dest,const DynamicSystem&DS_src)
{
memcpy(DS_dest.Mass,DS_src.Mass,sizeof(double)*DS_src.ParticleNumber);
memcpy(DS_dest.position,DS_src.position,sizeof(VEC)*DS_src.ParticleNumber);
memcpy(DS_dest.velocity,DS_src.velocity,sizeof(VEC)*DS_src.ParticleNumber);
memcpy(DS_dest.acceleration,DS_src.acceleration,sizeof(VEC)*DS_src.ParticleNumber);
memcpy(DS_dest.ParticleLabel.seq,DS_src.ParticleLabel.seq,sizeof(int)*DS_src.ParticleNumber);
memcpy(DS_dest.ParticleLabel.re_seq,DS_src.ParticleLabel.re_seq,sizeof(int)*DS_src.ParticleNumber);
DS_dest.ParticleNumber=DS_src.ParticleNumber;
}

void EvolutionTracking(FILE*EvolutionTrackFile,const double &time,const DynamicSystem& DS,const int Number)
{
    fprintf(EvolutionTrackFile,"%lf ",time);
    for(int i = 0 ; i<Number; i++)
    {
        fprintf(EvolutionTrackFile,"%lf %lf %lf %lf %lf %lf ",DS.position[DS.ParticleLabel.seq[i]].x/AU,DS.position[DS.ParticleLabel.seq[i]].y/AU,DS.position[DS.ParticleLabel.seq[i]].z/AU,DS.velocity[DS.ParticleLabel.seq[i]].x,DS.velocity[DS.ParticleLabel.seq[i]].y,DS.velocity[DS.ParticleLabel.seq[i]].z);
    }
    fprintf(EvolutionTrackFile,"\r\n");
}

void LastStatusKeep(FILE*DestinyInformatioFile,const DynamicSystem& DS,const int Number)
{
    for(int i = 0 ; i<Number; i++)
    {
        fprintf(DestinyInformatioFile,"%lf %lf %lf %lf %lf %lf ",DS.position[DS.ParticleLabel.seq[i]].x/AU,DS.position[DS.ParticleLabel.seq[i]].y/AU,DS.position[DS.ParticleLabel.seq[i]].z/AU,DS.velocity[DS.ParticleLabel.seq[i]].x,DS.velocity[DS.ParticleLabel.seq[i]].y,DS.velocity[DS.ParticleLabel.seq[i]].z);
    }
    fprintf(DestinyInformatioFile,"\r\n");
}

void WriteToFileBit(FILE*DataRecord_file,const double& time,VEC*position,VEC*velocity,const int ParticleNumber)
{
    fwrite(&time,sizeof(double),1,DataRecord_file);
    fwrite(position,sizeof(VEC),ParticleNumber,DataRecord_file);
    fwrite(velocity,sizeof(VEC),ParticleNumber,DataRecord_file);
}

void GetCentroid(const VEC&position1,const  VEC&position2,const double *M,VEC&centroidv)
{
    centroidv.x = (position1.x*M[0]+position2.x*M[1])/(M[0]+M[1]);
    centroidv.y = (position1.y*M[0]+position2.y*M[1])/(M[0]+M[1]);
    centroidv.z = (position1.z*M[0]+position2.z*M[1])/(M[0]+M[1]);
}

void DynamicSystemInitial(DynamicSystem& DS,double* M,const int IniParticleNumber)
{
    DS.ParticleNumber = IniParticleNumber;
    DS.position = new VEC[IniParticleNumber];
    DS.velocity = new VEC[IniParticleNumber];
    DS.acceleration = new VEC[IniParticleNumber];
    DS.Mass = new double[IniParticleNumber];
    MatterSet(DS.Mass,M,IniParticleNumber);
    DS.ParticleLabel.seq = new int[IniParticleNumber];
    DS.ParticleLabel.re_seq = new int[IniParticleNumber];
    SetNatureParticleLabel(DS.ParticleLabel,IniParticleNumber);
}

void DestroyDynamicSystem(DynamicSystem& DS)
{
    delete DS.Mass;
    delete DS.ParticleLabel.seq;
    delete DS.ParticleLabel.re_seq;
    delete DS.position;
    delete DS.velocity;
    delete DS.acceleration;
    return;
}
void DestinyFinish(FILE* DestinyInformatioFile,const double& CurrentTime,DynamicSystem&DS,const int IniParticleNumber)
{
    if(DS.ParticleNumber<IniParticleNumber)
    {
        fprintf(DestinyInformatioFile,"%d %lf ",UNEVENTFUL,CurrentTime);
        LastStatusKeep(DestinyInformatioFile,DS,IniParticleNumber);
    }
    else
    {
        fprintf(DestinyInformatioFile,"%d %lf %d %lf ",UNEVENTFUL,CurrentTime,UNEVENTFUL,CurrentTime);
        LastStatusKeep(DestinyInformatioFile,DS,IniParticleNumber);
    }
}
int JudgeDestiny(FILE* DestinyInformatioFile,const double& CurrentTime,DynamicSystem&DS,const int IniParticleNumber)
{
    static double DistanceMeasure;
    static double TEDR1 = GetTEDRadius(3,DS.Mass[2],DS.Mass[0]);
    static double TEDR2 = GetTEDRadius(3,DS.Mass[2],DS.Mass[1]);
    static double MergeLimit = GetMergeLimit(DS.Mass[0],DS.Mass[1]);
    static int destiny = 0;

    if(destiny==UNEVENTFUL)
    {
        DistanceMeasure = GetDistance(DS.position[0],DS.position[1]);
        if(DistanceMeasure < MergeLimit)
        {
            fprintf(DestinyInformatioFile,"%d %lf %d %lf ",MERGER,CurrentTime,UNEVENTFUL,CurrentTime);
            LastStatusKeep(DestinyInformatioFile,DS,IniParticleNumber);
            return END_EVOLUTION;
        }
    }
    if(DS.ParticleLabel.seq[0]<DS.ParticleNumber)
    {
        DistanceMeasure = GetDistance(DS.position[DS.ParticleLabel.seq[0]],DS.position[DS.ParticleLabel.seq[2]]);
        if(DistanceMeasure < TEDR1)
        {
            if(CurrentTime==0)
            {
                fprintf(DestinyInformatioFile,"%d %lf %d %lf ",INITIAL_INVALID,CurrentTime,INITIAL_INVALID,CurrentTime);
                LastStatusKeep(DestinyInformatioFile,DS,IniParticleNumber);
                return END_EVOLUTION;
            }
            else
            {
                destiny+=M0_TDE;
                fprintf(DestinyInformatioFile,"%d %lf ",M0_TDE,CurrentTime);
                ResetParticle(DS,0);
                if(DS.ParticleNumber<=IniParticleNumber-2)
                {
                    LastStatusKeep(DestinyInformatioFile,DS,IniParticleNumber);
                    return END_EVOLUTION;
                }
            }
        }
    }
    if(DS.ParticleLabel.seq[1]<DS.ParticleNumber)
    {
        DistanceMeasure = GetDistance(DS.position[DS.ParticleLabel.seq[1]],DS.position[DS.ParticleLabel.seq[2]]);
        if(DistanceMeasure < TEDR2)
        {
            if(CurrentTime==0)
            {
                fprintf(DestinyInformatioFile,"%d %lf %d %lf ",INITIAL_INVALID,CurrentTime,INITIAL_INVALID,CurrentTime);
                LastStatusKeep(DestinyInformatioFile,DS,IniParticleNumber);
                return END_EVOLUTION;
            }
            else
            {
                destiny+=M1_TDE;
                fprintf(DestinyInformatioFile,"%d %lf ",M1_TDE,CurrentTime);
                ResetParticle(DS,1);
                if(DS.ParticleNumber<=IniParticleNumber-2)
                {
                    LastStatusKeep(DestinyInformatioFile,DS,IniParticleNumber);
                    return END_EVOLUTION;
                }
            }
        }
    }
}

const double ERR=1e-7;
bool errCheck_Orbits(const DynamicSystem&DS_try,const DynamicSystem&DS_htry,const int OrbitNumber)
{
    bool e=0,te;
    double dr;
    double r1,tmp;
    VEC centroid_try,centroid_htry;
    if(OrbitNumber==3)
    {
        dr=GetDistance(DS_try.position[0],DS_try.position[1]);
        r1=GetDistance(DS_htry.position[0],DS_htry.position[1]);
        e= fabs(dr-r1)/dr > ERR;

        dr=GetDistance(DS_try.position[0],DS_try.position[2]);
        r1=GetDistance(DS_htry.position[0],DS_htry.position[2]);
        te= fabs(dr-r1)/dr > ERR;
        e = e||te;
    }
    else if(OrbitNumber==4)
    {
        dr=GetDistance(DS_try.position[DS_try.ParticleLabel.seq[3]],DS_try.position[DS_try.ParticleLabel.seq[2]]);
        r1=GetDistance(DS_htry.position[DS_htry.ParticleLabel.seq[3]],DS_htry.position[DS_htry.ParticleLabel.seq[2]]);
        e = fabs(dr-r1)/dr > ERR;

        dr=GetDistance(DS_try.position[DS_try.ParticleLabel.seq[0]],DS_try.position[DS_try.ParticleLabel.seq[1]]);
        r1=GetDistance(DS_htry.position[DS_htry.ParticleLabel.seq[0]],DS_htry.position[DS_htry.ParticleLabel.seq[1]]);
        tmp=fabs(dr-r1)/dr;
        if(tmp!=0)
        {
            te =  tmp > ERR;
            e=e||te;
        }

        GetCentroid(DS_try.position[DS_try.ParticleLabel.seq[0]],DS_try.position[DS_try.ParticleLabel.seq[1]],DS_try.Mass,centroid_try);
        GetCentroid(DS_htry.position[DS_htry.ParticleLabel.seq[0]],DS_htry.position[DS_htry.ParticleLabel.seq[1]],DS_htry.Mass,centroid_htry);
        dr=GetDistance(centroid_try,DS_try.position[DS_try.ParticleLabel.seq[2]]);
        r1=GetDistance(centroid_htry,DS_htry.position[DS_htry.ParticleLabel.seq[2]]);
        tmp=fabs(dr-r1)/dr;
        if(tmp!=0)
        {
            te =  tmp > ERR;
            e=e||te;
        }
    }
    else
    {
        e=0;
    }
    return e;
}
