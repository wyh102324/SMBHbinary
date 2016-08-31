/*
                Self-Adaptive Code
                 Created by Y.H.Wang  Sept.7th.2015


*/

#include"interaction.cpp"
#include"functionsLib.cpp"
#include"DynamicSystem.cpp"
#include"Integrator.cpp"

//-------------------code parameters---------------
#define StepLengthMinLimit 1e-6
#define StepLengthMaxLimit 1e9
#define DataPointNumber 1000
//===========================================================================================

void Self_Adaptive_Method(double& CurrentTime, double&StepLength, DynamicSystem&DS)
{
    static DynamicSystem DS_tmp, DS_try, DS_htry;//try: try evolve one step with current steplength; htry: with half of the current steplength;
    if(CurrentTime == 0)//initialize three temporary dynamic systems might be used in adaptive method
    {
        DynamicSystemInitial(DS_tmp, DS.Mass, DS.ParticleNumber);
        DynamicSystemInitial(DS_try, DS.Mass, DS.ParticleNumber);
        DynamicSystemInitial(DS_htry, DS.Mass, DS.ParticleNumber);
    }
    else if(DS.ParticleNumber < DS_tmp.ParticleNumber)//copy the system when particles number changes: the sequence of particle label
    {
         StatusCopy(DS_tmp, DS);
         StatusCopy(DS_try, DS);
         StatusCopy(DS_htry, DS);
    }

    Integrator_Leapfrog(DS_try, StepLength, DS);//evlove one step with current steplength
    Integrator_Leapfrog(DS_tmp, StepLength*0.5, DS);//evlove two steps with a half of current steplength
    Integrator_Leapfrog(DS_htry, StepLength*0.5, DS_tmp);
    if(errCheck_Orbits(DS_try, DS_htry, DS_try.ParticleNumber))//check the error between the two methods
    {//if current err is larger than errLimit
        for(StepLength *= 0.5; errCheck_Orbits(DS_try, DS_htry, DS_try.ParticleNumber) && StepLength > StepLengthMinLimit ;StepLength *= 0.5) //dichotomise steplength until the error between one step and two half step is smaller than errLimit
        {
            StatusCopy(DS_try, DS_htry);
            Integrator_Leapfrog(DS_tmp, StepLength*0.5, DS);
            Integrator_Leapfrog(DS_htry, StepLength*0.5, DS_tmp);
        }
        CurrentTime += StepLength;
        StatusCopy(DS, DS_htry);
    }
    else
    {//if current err is smaller than errLimit
        for(StepLength *= 2; !errCheck_Orbits(DS_try, DS_htry, DS_try.ParticleNumber) && StepLength < StepLengthMaxLimit ;StepLength *= 2) //double steplength until the error between one step and two half step is about to bigger than errLimit
        {
            StatusCopy(DS_tmp, DS_try);
            Integrator_Leapfrog(DS_htry, StepLength*0.5, DS_tmp);
            Integrator_Leapfrog(DS_try, StepLength, DS);
        }
        StepLength *= 0.25;
        CurrentTime += StepLength;
        StatusCopy(DS, DS_tmp);
    }

}

void Dynamics(bool EvolutionTrackSwitch, FILE* DestinyInformatioFile, double* M, const double* Major_Semi_Axis, const double* Eccentricity, const double* Tilt, const double EndTime, const int IniParticleNumber)
{
    FILE*EvolutionTrackFile = 0;//file pointer to record status(v,x,t)
    if(EvolutionTrackSwitch == EVOLUTION_TRACK_ON)
    {
        char filename[100];
        sprintf(filename, "%lf %lf %lf.txt", Major_Semi_Axis[1]/AU, Tilt[0]/PI*180, Tilt[2]/PI*180);
        EvolutionTrackFile = fopen(filename, "w");
    }

    DynamicSystem DS;
    DynamicSystemInitial(DS, M, IniParticleNumber);

    if(IniParticleNumber == 3)//initialize x v a with special orbital parameters
        initParameter_3Orbits(DS, Major_Semi_Axis, Eccentricity, Tilt);
    else if(IniParticleNumber == 4)
        initParameter_4Orbits(DS, Major_Semi_Axis, Eccentricity, Tilt);
    else
        return;

    int  DataInterval = 10000;
    double CurrentTime = 0, StepLength = EndTime*1e-9;

    for(int counter;CurrentTime < EndTime;counter++)
    {
        if(JudgeDestiny(DestinyInformatioFile, CurrentTime, DS, IniParticleNumber) == END_EVOLUTION)//check the status of TDE, merger
        {
            if(EvolutionTrackSwitch == EVOLUTION_TRACK_ON && counter%DataInterval == 0)
                EvolutionTracking(EvolutionTrackFile, CurrentTime, DS, IniParticleNumber);
            return;
        }
        else
        {
            Self_Adaptive_Method(CurrentTime, StepLength, DS);
        }

        if(EvolutionTrackSwitch == EVOLUTION_TRACK_ON)
        {
            if(counter == 0)//set steps that the data outstream write
            {
                DataInterval = (int)EndTime/StepLength;
                DataInterval /= DataPointNumber;
            }
            if(counter%DataInterval == 0)
            {
                EvolutionTracking(EvolutionTrackFile, CurrentTime, DS, IniParticleNumber);//record the status(x,v,t)
            }
        }
    }

    DestinyFinish(DestinyInformatioFile, CurrentTime, DS, IniParticleNumber);//check the status of TDE, merger of the system at the time limit.

    if(EvolutionTrackSwitch == EVOLUTION_TRACK_ON)
    {
        EvolutionTracking(EvolutionTrackFile, CurrentTime, DS, IniParticleNumber);//record the status(x,v,t) at the time limit
        fclose(EvolutionTrackFile);
    }

    DestroyDynamicSystem(DS);
    return;
}

