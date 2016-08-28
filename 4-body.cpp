#include<time.h>
#include<math.h>
#include<string.h>
#include<stdio.h>
#include<stdlib.h>

#include"astroSimulation.h"
#include"self-adapted.cpp"

#define m1 (2*m_sun)
#define m2 (1*m_sun)
#define m3 (1e6*m_sun)
#define m4 (1e6*m_sun)

#define e1 (0.001)
#define e2 (0.01)
#define e3 (0.01)

#define a1 (0.05*AU)
#define a2_min (1*AU)
#define a2_max (150*AU)
#define a3 (400*AU)


#define i2 (0.0/180*PI)


void CreateInformationFile(double* Mass, double* Major_Semi_Axis, double* Eccentricity)
{
   time_t GetTime = time(NULL);
   tm* CreatedTime = localtime(&GetTime);
   FILE*Information_file = fopen("ParameterInfo.txt", "w");
   fprintf(Information_file, "-------Created in %d %02d %02d %02d:%02d:%02d---------\r\n", CreatedTime->tm_year + 1900, CreatedTime->tm_mon + 1, CreatedTime->tm_mday, CreatedTime->tm_hour, CreatedTime->tm_min, CreatedTime->tm_sec);
   fprintf(Information_file, "Major-Semi-Axis : \r\n a1 = %lf AU \r\n a3 = %lf AU \r\nEccentricity : \r\n e1 = %lf \r\n \
e2 = %lf \r\n e3 = %lf\r\nMass : \r\n m1 = %lf M_sun \r\n m2 = %lf M_sun \r\n m3 = %lf M_sun \r\n m4 = %lf M_sun\r\n" , Major_Semi_Axis[0]/AU, Major_Semi_Axis[2]/AU, Eccentricity[0], Eccentricity[1], Eccentricity[2], Mass[0]/m_sun, Mass[1]/m_sun, Mass[2]/m_sun, Mass[3]/m_sun);
   fclose(Information_file);
}


void RunSamples(int SampleNumber, double* Mass, double* Major_Semi_Axis, double* Eccentricity, double* Tilt)
{
    double t_inner, t_out,  EndTime = 0;
    FILE* DestinyInformatioFile = fopen("Destiny.txt", "w");

    srand( (unsigned)time( NULL ) );
    for(int i = 0;i < SampleNumber;i++)
    {
        Major_Semi_Axis[1] = a2_min + ( (double)rand()/RAND_MAX)*(a2_max - a2_min);
        Tilt[0] = ( (double)rand()/RAND_MAX )*PI;
        Tilt[2] = ( (double)rand()/RAND_MAX )*PI;
        t_inner = (Mass[0] + Mass[1])/Mass[2]*pow((Major_Semi_Axis[1]/Major_Semi_Axis[0]), 3)*pow( (1 - Eccentricity[1]*Eccentricity[1]), 1.5)/sqrt(G*(Mass[0] + Mass[1])/pow(Major_Semi_Axis[0], 3));
        t_out = (Mass[0] + Mass[1] + Mass[2])/Mass[3]*pow((Major_Semi_Axis[2]/Major_Semi_Axis[1]), 3)*pow( (1 - Eccentricity[2]*Eccentricity[2]), 1.5)/sqrt(G*(Mass[0] + Mass[1] + Mass[2] + Mass[3])/pow(Major_Semi_Axis[1], 3));
        EndTime = 0.1*(t_inner + t_out);

        fprintf(DestinyInformatioFile, "%lf %lf %lf ", Major_Semi_Axis[1]/AU, Tilt[0]*180/PI, Tilt[2]*180/PI);
        Dynamics(EVOLUTION_TRACK_OFF, DestinyInformatioFile,  Mass, Major_Semi_Axis, Eccentricity, Tilt, EndTime, 4);
    }
}

int main()
{
    double Mass[4] = {m1, m2, m3, m4};
    double Major_Semi_Axis[3] = {a1, 0, a3};
    double Eccentricity[3] = {e1, e2, e3};
    double Tilt[3] = {0, i2, 0};

    CreateInformationFile(Mass, Major_Semi_Axis, Eccentricity);
    RunSamples(5, Mass, Major_Semi_Axis, Eccentricity, Tilt);

    return 0;
}
