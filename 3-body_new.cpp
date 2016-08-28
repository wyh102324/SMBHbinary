#include<time.h>
#include<math.h>
#include<string.h>
#include<stdio.h>
#include<stdlib.h>

#include"astroSimulation.h"
#include"selfAdaptiveMethod.cpp"

#define m1 (2*m_sun)
#define m2 (1*m_sun)
#define m3 (1e6*m_sun)


#define e1 (0.001)


#define a1 (0.05*AU)
#define a2_min (1*AU)
#define a2_max (250*AU)



#define i2 (0.0/180*PI)


void CreateInformationFile(double*Mass,double*Major_Semi_Axis,double*Eccentricity)
{
   time_t GetTime = time(NULL);
   tm* CreatedTime = localtime(&GetTime);
   FILE*Information_file = fopen("ParameterInfo.txt","w");
   fprintf(Information_file,"-------Created in %d %02d %02d %02d:%02d:%02d---------\r\n",CreatedTime->tm_year + 1900,CreatedTime->tm_mon + 1,CreatedTime->tm_mday,CreatedTime->tm_hour,CreatedTime->tm_min,CreatedTime->tm_sec);
   fprintf(Information_file,"Major-Semi-Axis : \r\n a1=%lf AU\r\n a2=[1,250] AU\r\nEccentricity : \r\n e1=%lf \r\n \
e2=[0,1] \r\nMass : \r\n m1=%lf M_sun \r\n m2=%lf M_sun \r\n m3=%lf M_sun \r\n",Major_Semi_Axis[0]/AU,Eccentricity[0],Mass[0]/m_sun,Mass[1]/m_sun,Mass[2]/m_sun);
   fclose(Information_file);
}


void RunSamples(int SampleNumber,double* Mass,double* Major_Semi_Axis,double* Eccentricity,double* Tilt)
{
    double t_inner,t_out, EndTime=0;
    FILE* DestinyInformatioFile=fopen("Destiny.txt","w");

    srand( (unsigned)time( NULL ) );
    for(int i=0;i<SampleNumber;i++)
    {
        Major_Semi_Axis[1]= (a2_min*a2_min+(a2_max*a2_max-a2_min*a2_min)*( (double)rand()/RAND_MAX));
        Major_Semi_Axis[1]=sqrt(Major_Semi_Axis[1]);
        Tilt[0]=( (double)rand()/RAND_MAX )*PI;
        Tilt[1]=( (double)rand()/RAND_MAX )*PI;
        Eccentricity[1]=( (double)rand()/RAND_MAX);
        Eccentricity[1]*=Eccentricity[1];
        Eccentricity[1]=sqrt(1-Eccentricity[1]);
        EndTime=2*PI*sqrt(pow(Major_Semi_Axis[1],3)/G/(Mass[0]+Mass[1]+Mass[2]));

        fprintf(DestinyInformatioFile,"%lf %lf %lf %lf ",Major_Semi_Axis[1]/AU,Eccentricity[1],Tilt[0]*180/PI,Tilt[1]*180/PI);
        printf("%lf %lf %lf %lf \r\n",Major_Semi_Axis[1]/AU,Eccentricity[1],Tilt[0]*180/PI,Tilt[1]*180/PI);
        Dynamics(EVOLUTION_TRACK_OFF,DestinyInformatioFile, Mass, Major_Semi_Axis, Eccentricity, Tilt, EndTime, 3);
    }
}

int main()
{
    double Mass[3]={m1,m2,m3};
    double Major_Semi_Axis[2]={a1,0};
    double Eccentricity[2]={e1,0};
    double Tilt[2]={0,0};

    CreateInformationFile(Mass,Major_Semi_Axis,Eccentricity);
    RunSamples(1000,Mass,Major_Semi_Axis,Eccentricity,Tilt);

    return 0;
}
