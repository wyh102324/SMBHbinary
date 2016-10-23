
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
                /*DS.acceleration[i].x += GetForce(DS.position[i], DS.position[j], DS_src.Mass[j], 1);
                DS.acceleration[i].y += GetForce(DS.position[i], DS.position[j], DS_src.Mass[j], 2);
                DS.acceleration[i].z += GetForce(DS.position[i], DS.position[j], DS_src.Mass[j], 3);*/
                DS.acceleration[i].x += GetPNForce(DS.position[i], DS.position[j], DS_src.Mass[i], DS_src.Mass[j], 1);
                DS.acceleration[i].y += GetPNForce(DS.position[i], DS.position[j], DS_src.Mass[i], DS_src.Mass[j], 2);
                DS.acceleration[i].z += GetPNForce(DS.position[i], DS.position[j], DS_src.Mass[i], DS_src.Mass[j], 3);
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

