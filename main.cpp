#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <mpi.h>
#include <math.h>

#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"

#define epsilon 0.000000000000000222

int main(int argc, char* argv[]){
    
    bool inputError = false;
    
	if( argc != 10){
		printf("Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrefix\n", argv[0]);
	}
    else if((int)atof(argv[1]) < 0){
        printf("Number of Light Particles can't be negative.");
        inputError = true;
    }
    else if((int) atof(argv[2]) < 0){
        printf("Number of Medium Particles can't be negative.");
        inputError = true;
    }
    else if((int) atof(argv[3]) < 0){
        printf("Number of Heavy Particles can't be negative.");
        inputError = true;
    }
    else if((int) atof(argv[4]) < 0){
        printf("Number of steps can't be negative.");
        inputError = true;
    }
    else if((int) atof(argv[5]) < 0){
        printf("Number of substeps can't be negative.");
        inputError = true;
    }
    else if((double) atof(argv[6]) <= 0){
        printf("Time must be greater than 0 seconds.");
        inputError = true;
    }
    else if((int) atof(argv[7]) <= 0){
        printf("Image width must be greter than 0.");
        inputError = true;
    }
    else if((int) atof(argv[8]) <= 0){
        printf("Image height must be greter than 0.");
        inputError = true;
    }

    
    
	MPI_Init(&argc,&argv);

	int p, my_rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	//variables
	int numParticlesLight = 0;
	int numParticleMedium = 0;
	int numParticleHeavy = 0;

	int numSteps = 0;
	int subSteps = 0;
	double timeSubStep;

	int width, height;

	unsigned char* image;

	//root node stuff goes here
	if(my_rank == 0){
        numParticlesLight = atoi(argv[1]);
        numParticleMedium = atoi(argv[2]);
        numParticleHeavy  = atoi(argv[3]);
        
        numSteps = atoi(argv[4]);
        subSteps = atoi(argv[5]);
        
        timeSubStep = (double) atof(argv[6]);
        
        width = atoi(argv[7]);
        height = atoi(argv[8]);
        
        
        
		//almost done, just save the image
		//saveBMP(argv[9], image, width, height);
	}
	//all other nodes do this
	else{

	}

	free(image);

	MPI_Finalize();
	return 0;
}