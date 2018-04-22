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
        int subSteps = 0;
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
            subSteps = atoi(argv[5]);
            
            timeSubStep = (double) atof(argv[6]);
            
            width = atoi(argv[7]);
            height = atoi(argv[8]);
            
            //depth of z plane
            const int zplane=100;
            
            //print args to the console
            printf("numParticlesLight: %d\nnumParticlesMedium: %d\nnumParticlesLarge: %d\nnumSteps: %d\nsubSteps: %d\ntimeSubStep: %f\nWidth: %d\nHeight: %d\n",numParticlesLight,numParticleMedium,numParticleHeavy,numSteps,subSteps,timeSubStep, width,height);
            
            
            /* =====================================================================
             
                                            Setup
             
             ===================================================================== */
            
            const int totalParticles=numParticlesLight+numParticleMedium+numParticleHeavy;
            
            int frameSize=width*height;
            
            //allocate space for pixel colours
            image = (unsigned char *)malloc(3*frameSize);
            
            //default all pixels to white/black
            for(int a=0; a<(3*frameSize);a++){//a+=3
                image[a]=(unsigned char) 0;//0;
                //image[a+1]=(unsigned char) 0;
                //image[a+2]=(unsigned char) 0;
            }
            
            vec3 particles [totalParticles];
            
            //create all three sizes of pixels
            int count=0;
            for(int i=0; i<numParticlesLight; i++){
                double xVal = (drand48() * (width-1));
                double yVal = (drand48() * (height-1));
                double zVal = (drand48() * (zplane));
                //printf("Coordinates X,Y,Z: %d %d %d\n", xVal, yVal, zVal);
                particles[i] = vec3(xVal,yVal,zVal,255,0,0); //vec3(1,1,1,255,255,255);
                
            }
            count=numParticlesLight;
            
            for(int j=count; j<(numParticleMedium+count); j++){
                int xVal = (int) (drand48() * (width-1));
                int yVal = (int) (drand48() * (height-1));
                int zVal = (int) (drand48() * (zplane));
                particles[j] = vec3(xVal,yVal,zVal,0,255,0);//vec3(10,10,10,255,255,255);
                
            }
            count=numParticlesLight+numParticleMedium;
            
            for(int k=count; k<(numParticleHeavy+count); k++){
                int xVal = (int) (drand48() * (width-1));
                int yVal = (int) (drand48() * (height-1));
                int zVal = (int) (drand48() * (zplane));
                particles[k] = vec3(xVal,yVal,zVal,0,0,255);//vec3(100,100,100,255,255,255);
                
            }

            
            //write pixels onto screen
            for(int g=0;g<totalParticles;g++){
                int xValue=(int) particles[g].getX();
                int yValue=(int) particles[g].getY();
                
                image[(xValue*3)+(yValue*width*3)]=(unsigned char) particles[g].getR();
                image[(xValue*3)+(yValue*width*3)+1]=(unsigned char) particles[g].getG();
                image[(xValue*3)+(yValue*width*3)+2]=(unsigned char) particles[g].getB();
                
                
                printf("this is %d col %d %d %d %f %f %f \n",g, particles[g].getR(),particles[g].getG(), particles[g].getB(), particles[g].getX(), particles[g].getY(), particles[g].getZ());
            }
            
            char file[20];
            strcpy(file, argv[9]);
            strcat(file,"_00000.bmp");
            const char* filename = file;
            
            
            //printf("this is arg 9: %s %s %s \n", argv[9],argv[8], argv[7]);
            const unsigned char* result = (image);
            saveBMP (filename, result, width, height);
            //saveBMP      (const char* filename, const unsigned char* image, int width, int height);
            
            //free(image);

            
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