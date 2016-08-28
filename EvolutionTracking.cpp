#include<time.h>
#include<math.h>
#include<string.h>
#include<stdio.h>
#include<stdlib.h>

#include"astroSimulation.h"
#include"self-adapted.cpp"

#define m1 (1*m_sun)
#define m2 (0.25*m_sun)
#define m3 (0.6*m_sun)


#define e1 (0.01)
#define e2 (0.6)


#define a1 (60*AU)
#define a2 (800*AU)


#define i1 (90.02/180*PI)
#define i2 (7.98/180*PI)


void CreateInformationFile(double* Mass, double* Major_Semi_Axis, double* Eccentricity)
{
   time_t GetTime = time(NULL);
   tm* CreatedTime = localtime(&GetTime);
   FILE*Information_file = fopen("ParameterInfo.txt", "w");
   fprintf(Information_file, "-------Created in %d %02d %02d %02d:%02d:%02d---------\r\n", CreatedTime->tm_year + 1900, CreatedTime->tm_mon + 1, CreatedTime->tm_mday, CreatedTime->tm_hour, CreatedTime->tm_min, CreatedTime->tm_sec);
   fprintf(Information_file, "Major-Semi-Axis : \r\n a1 = %lf AU\r\n a2 = %lf AU\r\nEccentricity : \r\n e1 = %lf \r\n \
e2 = %lf \r\nMass : \r\n m1 = %lf M_sun \r\n m2 = %lf M_sun \r\n m3 = %lf M_sun \r\n", Major_Semi_Axis[0]/AU, Major_Semi_Axis[1]/AU, Eccentricity[0], Eccentricity[1], Mass[0]/m_sun, Mass[1]/m_sun, Mass[2]/m_sun);
   fclose(Information_file);
}


void RunSamples(int SampleNumber, double* Mass, double* Major_Semi_Axis, double* Eccentricity, double* Tilt)
{
    double t_inner, t_out, EndTime = 0;
    FILE* DestinyInformatioFile = fopen("Destiny.txt", "w");

    //EndTime = 2*PI*sqrt(pow(Major_Semi_Axis[1], 3)/G/(Mass[0] + Mass[1] + Mass[2]));
    EndTime = 3e6*24*3600*365;
    fprintf(DestinyInformatioFile, "%lf %lf %lf %lf ", Major_Semi_Axis[1]/AU, Eccentricity[1], Tilt[0]*180/PI, Tilt[2]*180/PI);
    printf("%lf %lf %lf %lf \r\n", Major_Semi_Axis[1]/AU, Eccentricity[1], Tilt[0]*180/PI, Tilt[2]*180/PI);
    Dynamics(EVOLUTION_TRACK_ON, DestinyInformatioFile,  Mass, Major_Semi_Axis, Eccentricity, Tilt, EndTime, 3);
}

int main()
{
    double Mass[3] = {m1, m2, m3};
    double Major_Semi_Axis[2] = {a1, a2};
    double Eccentricity[2] = {e1, e2};
    double Tilt[2] = {i1, i2};

    CreateInformationFile(Mass, Major_Semi_Axis, Eccentricity);
    RunSamples(1, Mass, Major_Semi_Axis, Eccentricity, Tilt);

    return 0;
}
