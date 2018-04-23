#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <mpi.h>
#include <math.h>
#include <cmath>

#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"

#define epsilon 0.000000000000000222

int main(int argc, char* argv[]){
    
    int inputError = 0;
    
    if( argc != 10){
        printf("Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrefix\n", argv[0]);
        inputError = 1;
    }
    else if((int)atof(argv[1]) < 0){
        printf("Number of Light Particles can't be negative.\n");
        inputError = 1;
    }
    else if((int) atof(argv[2]) < 0){
        printf("Number of Medium Particles can't be negative.\n");
        inputError = 1;
    }
    else if((int) atof(argv[3]) < 0){
        printf("Number of Heavy Particles can't be negative.\n");
        inputError = 1;
    }
    else if((int) atof(argv[4]) < 0){
        printf("Number of steps can't be negative.\n");
        inputError = 1;
    }
    else if((int) atof(argv[5]) < 0){
        printf("Number of substeps can't be negative.\n");
        inputError = 1;
    }
    else if((double) atof(argv[6]) <= 0){
        printf("Time 1 be greater than 0 seconds.\n");
        inputError = 1;
    }
    else if((int) atof(argv[7]) <= 0){
        printf("Image width must be greter than 0.\n");
        inputError = 1;
    }
    else if((int) atof(argv[8]) <= 0){
        printf("Image height must be greter than 0.\n");
        inputError = 1;
    }
    
    
    
//    MPI_Init(&argc,&argv);
//    
//    int p, my_rank;
//    
//    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
//    MPI_Comm_size(MPI_COMM_WORLD, &p);
//    
//    //variables
//    int numParticlesLight = 0;
//    int numParticleMedium = 0;
//    int numParticleHeavy = 0;
//    
//    int numSteps = 0;
//    int subSteps = 0;
//    double timeSubStep;
//    
//    int width, height;
//    
//    
//    unsigned char* image;
    
//      //root node stuff goes here
//      if(my_rank == 0){
//            numParticlesLight = atoi(argv[1]);
//            numParticleMedium = atoi(argv[2]);
//            numParticleHeavy  = atoi(argv[3]);
//    
//            numSteps = atoi(argv[4]);
//            subSteps = atoi(argv[5]);
//    
//            timeSubStep = (double) atof(argv[6]);
//    
//            width = atoi(argv[7]);
//            height = atoi(argv[8]);
//    
//    
//    
//          //almost done, just save the image
//          //saveBMP(argv[9], image, width, height);
//      }
//      //all other nodes do this
//      else{
//    
//      }
    
   
    
    if(inputError == 0){
        
        MPI_Init(&argc,&argv);
        
        int p, my_rank;
        
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &p);
        
        //variables
        int numParticlesLight = 0;
        int numParticleMedium = 0;
        int numParticleHeavy = 0;
        
        int numSteps = 0;
        int numSubSteps = 0;
        double timeSubStep;
        
        int width, height;
        
        
        unsigned char* image;

        
        //root node stuff goes here
        if(my_rank == 0){
            
            //grab arguments
            numParticlesLight = atoi(argv[1]);
            numParticleMedium = atoi(argv[2]);
            numParticleHeavy  = atoi(argv[3]);
            
            numSteps = atoi(argv[4]);
            numSubSteps = atoi(argv[5]);
            
            timeSubStep = (double) atof(argv[6]);
            
            width = atoi(argv[7]);
            height = atoi(argv[8]);
            
            //depth of z plane
            const int zplane=100;
            
            //print args to the console
            printf("numParticlesLight: %d\nnumParticlesMedium: %d\nnumParticlesLarge: %d\nnumSteps: %d\nsubSteps: %d\ntimeSubStep: %f\nWidth: %d\nHeight: %d\n",numParticlesLight,numParticleMedium,numParticleHeavy,numSteps,numSubSteps,timeSubStep, width,height);
            
            
            /* =====================================================================
             
                                            Setup
             
             ===================================================================== */
            
            const int totalParticles=numParticlesLight+numParticleMedium+numParticleHeavy;
            
            int frameSize=width*height;
            
            //allocate space for pixel colours
            image = (unsigned char *)malloc(3*frameSize);
            
            //default all pixels to white/black
            
            memset(image, 0, sizeof(unsigned char) * 3 * frameSize);

            
            
            
            
//            for(int a=0; a<(3*frameSize);a++){//a+=3
//                image[a]=(unsigned char) 0;//0;
//                //image[a+1]=(unsigned char) 0;
//                //image[a+2]=(unsigned char) 0;
//            }
            
            vec3 particles [totalParticles];
            int count=0;
            double velocity=0;
            double mass=0;
            double xVal = 0;
            double yVal = 0;
            double zVal = 0;
            int direction = 0;
            
            /*  === Directions Layout ===
             
             x       y       z     dir (int)
             1       0       0  --> 0
             0       1       0  --> 1
             0       0       1  --> 2
             1       1       0  --> 3
             1       0       1  --> 4
             0       1       1  --> 5
             1       1       1  --> 6
             
             ================ */
            
            
            
            const double G = -0.00000000006673;
            
            //create all three sizes of pixels
            
            for(int i=0; i<numParticlesLight; i++){
                velocity=drand48()*(velocityLightMax-velocityLightMin)+velocityLightMin;
                mass=drand48()*(massLightMax-massLightMin)+massLightMin;
                xVal = (drand48() * (width-1));
                yVal = (drand48() * (height-1));
                zVal = (drand48() * (zplane));
                direction = (int) (drand48() * 6.99);
                //printf("Coordinates X,Y,Z: %d %d %d\n", xVal, yVal, zVal);
                particles[i] = vec3(xVal,yVal,zVal,255,0,0,mass,velocity,direction); //vec3(1,1,1,255,255,255);
                
            }
            count=numParticlesLight;
            
            for(int j=count; j<(numParticleMedium+count); j++){
                velocity=drand48()*(velocityMediumMax-velocityMediumMin)+velocityMediumMin;
                mass=drand48()*(massMediumMax-massMediumMin)+massMediumMin;
                xVal = (drand48() * (width-1));
                yVal = (drand48() * (height-1));
                zVal = (drand48() * (zplane));
                direction = (int) (drand48() * 6.99);
                particles[j] = vec3(xVal,yVal,zVal,0,255,0,mass,velocity,direction);//vec3(10,10,10,255,255,255);
                
            }
            count=numParticlesLight+numParticleMedium;
            
            for(int k=count; k<(numParticleHeavy+count); k++){
                velocity=drand48()*(velocityHeavyMax-velocityHeavyMin)+velocityHeavyMin;
                mass=drand48()*(massHeavyMax-massHeavyMin)+massHeavyMin;
                xVal = (drand48() * (width-1));
                yVal = (drand48() * (height-1));
                zVal = (drand48() * (zplane));
                direction = (int) (drand48() * 6.99);
                particles[k] = vec3(xVal,yVal,zVal,0,0,255,mass,velocity,direction);//vec3(100,100,100,255,255,255);
                
            }

            
            //write pixels onto screen
            for(int g=0;g<totalParticles;g++){
                int xValue=(int) particles[g].getX();
                int yValue=(int) particles[g].getY();
                image[(xValue*3)+(yValue*width*3)]=(unsigned char) particles[g].getR();
                image[(xValue*3)+(yValue*width*3)+1]=(unsigned char) particles[g].getG();
                image[(xValue*3)+(yValue*width*3)+2]=(unsigned char) particles[g].getB();
                
            }
            
            char file[20];
            strcpy(file, argv[9]);
            strcat(file,"_00000.bmp");
            const char* filename = file;
            
            
            //printf("this is arg 9: %s %s %s \n", argv[9],argv[8], argv[7]);
            const unsigned char* result = (image);
            saveBMP (filename, result, width, height);
            
            
            
            
            
            /* =====================================================================
             
             Computing velocitites and positions for next frame
             
             ===================================================================== */
            
            
            
            
            
            double outerDotMagSquared = 0;
            //double forcesArr = (double *)malloc(totalParticles*totalParticles);
            
            double * forcesArr;
            forcesArr = (double *)malloc(totalParticles*totalParticles);
            memset(forcesArr, 0, sizeof(double) * totalParticles*totalParticles);
            double currentForce = 0;
            //forces = totalParticles*totalParticles*sizeof(double);
            
            for(int step = 0; step < numSteps; step++){
                for(int subStep = 0; subStep < numSubSteps; subStep++){
                    for(int outerDot = 0; outerDot < totalParticles; outerDot++){
                        currentForce = 0;
                        outerDotMagSquared = particles[outerDot].MagnitudeSquared();
                        
                        for(int innerDot = outerDot+1; innerDot < totalParticles; innerDot++){
                            
                            currentForce = particles[innerDot].getMass()/(abs(outerDotMagSquared - particles[innerDot].MagnitudeSquared()));
                            forcesArr[innerDot + width*outerDot] = currentForce;
                            forcesArr[outerDot + width*innerDot] = -1 * currentForce;
                            
                        }
                        
                    }
                }
            }
            
            
            double forces = 0;
            double totForce =  0;
            double currentMass = 0;
            for(int index = 0; index < totalParticles; index++){
                currentMass = particles[index].getMass();
                totForce = 0;
                forces = 0;
                for(int length = 0; length < width; length ++){
                    forces += forcesArr[length + index*width];
                }
                
                totForce = G*forces*currentMass;
                //printf("Position before -> X: %f \t Y: %f \t Z: %f --- Direction: %d \t\t", particles[index].getX(), particles[index].getY(), particles[index].getZ(), particles[index].getDirection());
                particles[index].setPosition( particles[index], timeSubStep );
                //printf("Position AFTER -> X: %f \t Y: %f \t Z: %f\n", particles[index].getX(), particles[index].getY(), particles[index].getZ());
//                if(particles[index].getDirection() == 0){
//                    printf("getX(): %f\n", particles[index].getX());
//                }
                
                //particles[index].getVelocity() + ((timeSubStep * totForce) / currentMass)
                printf("force: %f  --- first velocity: %f    \t", totForce, particles[index].getVelocity());
                particles[index].setVelocity( particles[index], timeSubStep, totForce );
                printf("after: %f \n", particles[index].getVelocity());
            }
            
            for(int a=0; a<(3*frameSize);a++){//a+=3
                image[a]=(unsigned char) 0;//0;
                //image[a+1]=(unsigned char) 0;
                //image[a+2]=(unsigned char) 0;
            }
            
            
            //write pixels onto screen
            for(int g=0;g<totalParticles;g++){
                int xValue=(int) particles[g].getX();
                int yValue=(int) particles[g].getY();
                if(!(particles[g].getX() > width || particles[g].getX() < 0 || particles[g].getY() > height || particles[g].getY() < 0)){
                    image[(xValue*3)+(yValue*width*3)]=(unsigned char) particles[g].getR();
                    image[(xValue*3)+(yValue*width*3)+1]=(unsigned char) particles[g].getG();
                    image[(xValue*3)+(yValue*width*3)+2]=(unsigned char) particles[g].getB();
                }
                
            }
            
            char file1[20];
            strcpy(file1, argv[9]);
            strcat(file1,"_00001.bmp");
            const char* filename1 = file1;
            
            
            //printf("this is arg 9: %s %s %s \n", argv[9],argv[8], argv[7]);
            const unsigned char* result1 = (image);
            saveBMP (filename1, result1, width, height);
            //saveBMP      (const char* filename, const unsigned char* image, int width, int height);
            free(forcesArr);
            free(image);

            
            //almost done, just save the image
            //saveBMP(argv[9], image, width, height);
        }
        //all other nodes do this
        else{
            
        }

        
        MPI_Finalize();
        
    } else{
        printf("\n--> Program did not run do to input error. Please view message above and fix input arguments.\n");
    }
    
    
    
    return 0;
}