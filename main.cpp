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
    
    
    if( argc != 10){
        printf("Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrefix\n", argv[0]);
    }
    else if((int)atof(argv[1]) < 0){
        printf("Number of Light Particles can't be negative.\n");
    }
    else if((int) atof(argv[2]) < 0){
        printf("Number of Medium Particles can't be negative.\n");
    }
    else if((int) atof(argv[3]) < 0){
        printf("Number of Heavy Particles can't be negative.\n");
    }
    else if((int) atof(argv[4]) < 0){
        printf("Number of steps can't be negative.\n");
    }
    else if((int) atof(argv[5]) < 0){
        printf("Number of substeps can't be negative.\n");
    }
    else if((double) atof(argv[6]) <= 0){
        printf("Time must be greater than 0 seconds.\n");
        printf("Time 1 be greater than 0 seconds.\n");
    }
    else if((int) atof(argv[7]) <= 0){
        printf("Image width must be greter than 0.\n");
    }
    else if((int) atof(argv[8]) <= 0){
        printf("Image height must be greter than 0.\n");
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
        
<<<<<<< HEAD

    const int totalParticles=numParticlesLight+numParticleMedium+numParticleHeavy;
    
    int frameSize=width*height;
    
    //allocate space for pixel colours
    image = (unsigned char *)malloc(3*frameSize);

    //default all pixels to black/white 
    for(int a=0; a<(3*frameSize);a++){//a+=3
        image[a]=(unsigned char) 0;//0;
        //image[a+1]=(unsigned char) 0;
        //image[a+2]=(unsigned char) 0;
    }   

    vec3 particles [totalParticles];
    
    //create all three sizes of pixels
    int count=0;
    double velocity=0;
    double mass=0;
    double X=0;
    double Y=0;
    double Z=0;
    
    for(int i=0; i<numParticlesLight; i++){
        velocity=drand48()*(velocityLightMax-velocityLightMin)+velocityLightMin;
        mass=drand48()*(massLightMax-massLightMin)+massLightMin;
        X=drand48()*(width);
        Y=drand48()*(height);
        Z=drand48()*(zplane);

        //printf("this is the x %f, this is the y %f, this is the z %f \n",X,Y,Z);
        particles[i] = vec3(X,Y,Z,255,255,255,mass,velocity); //vec3(1,1,1,255,255,255,mass,velocity);
=======
        int width, height;
>>>>>>> dd9de02f3543ecadd5a588ce92a0db8715f94856
        
        
<<<<<<< HEAD
    for(int j=count; j<(numParticleMedium+count); j++){
        velocity=drand48()*(velocityMediumMax-velocityMediumMin)+velocityMediumMin;
        mass=drand48()*(massMediumMax-massMediumMin)+massMediumMin;
        X=drand48()*(width);
        Y=drand48()*(height);
        Z=drand48()*(zplane);

        //printf("this is the x %f, this is the y %f, this is the z %f \n",X,Y,Z);
        particles[j] = vec3(X,Y,Z,255,255,255,mass,velocity);//vec3(10,10,10,255,255,255,mass,velocity);
=======
        unsigned char* image;

>>>>>>> dd9de02f3543ecadd5a588ce92a0db8715f94856
        
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
                int xVal = (int) (drand48() * (width-1));
                int yVal = (int) (drand48() * (height-1));
                int zVal = (int) (drand48() * (zplane));
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

<<<<<<< HEAD
    for(int k=count; k<(numParticleHeavy+count); k++){
        velocity=drand48()*(velocityHeavyMax-velocityHeavyMin)+velocityHeavyMin;
        mass=drand48()*(massHeavyMax-massHeavyMin)+massHeavyMin;
        X=drand48()*(width);
        Y=drand48()*(height);
        Z=drand48()*(zplane);

        //printf("this is the x %f, this is the y %f, this is the z %f \n",X,Y,Z);
        particles[k] = vec3(X,Y,Z,255,255,255,mass,velocity);//vec3(100,100,100,255,255,255,mass,velocity);
        
    }
     
    


    //printf("testing \n");
    //write pixels onto screen
    for(int g=0;g<totalParticles;g++){
        int xValue=(int) particles[g].getX();
        int yValue=(int) particles[g].getY();
        
        //printf("this is the x %d, this is the y %d \n",xValue,yValue);        

        image[(xValue*3)+(yValue*width*3)]=(unsigned char) particles[g].getR();
        image[(xValue*3)+(yValue*width*3)+1]=(unsigned char) particles[g].getG();
        image[(xValue*3)+(yValue*width*3)+2]=(unsigned char) particles[g].getB();


        //printf("this is %d col %d %d %d  \n",g, particles[g].getR(),particles[g].getG(), particles[g].getB());
=======
        
        MPI_Finalize();
        
    } else{
        printf("\n--> Program did not run do to input error. Please view message above and fix input arguments.\n");
>>>>>>> dd9de02f3543ecadd5a588ce92a0db8715f94856
    }
    
    
    
<<<<<<< HEAD

    //Equations of motion - performing steps and substeps
    //at each substep calculate, at each step take a photo








        free(image);
    }


    //MPI_Finalize();
=======
>>>>>>> dd9de02f3543ecadd5a588ce92a0db8715f94856
    return 0;
}