#ifndef properties_h
#define properties_h

#include "vector3d.h"


//G was so small the precision was not large enough, in some cases the force needed to be replaced by epsilon
//solar mass was a way to scale the mass if necessary
const double solarMass = 1;


//light particles are the fastest
double velocityLightMin = 11;
double velocityLightMax = 15;

double velocityMediumMin = 6;
double velocityMediumMax = 10;

//heavy particles are the slowest
double velocityHeavyMin = 1;
double velocityHeavyMax = 5;

//mass
double massLightMin = 1*solarMass;
double massLightMax = 5*solarMass;

double massMediumMin = 6*solarMass;
double massMediumMax = 10*solarMass;

double massHeavyMin = 11*solarMass;
double massHeavyMax = 15*solarMass;


//colours
//vec3 colourLight = vec3(0,0,1);
//vec3 colourMedium = vec3(0,1,0);
//vec3 colourHeavy = vec3(1,0,0);

#endif