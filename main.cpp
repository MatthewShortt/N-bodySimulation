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
   
    
    if(inputError == 0){
        double start=0;
        double end=0;
        double timeRequired=0;
        double subStepTimes [atoi(argv[4])*atoi(argv[5])];
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
            image = (unsigned char *)malloc(3*frameSize*sizeof(unsigned char));
            
            //default all pixels to white/black
            
            //memset(image, 0, sizeof(unsigned char) * 3 * frameSize);

            for(int a=0; a<(3*frameSize);a++){
                image[a]=(unsigned char) 0;
 
            }
            
            vec3 particles [totalParticles];
            int count=0;
            double velocity=0;
            double velocityX=0;
            double velocityY=0;
            double velocityZ=0;
            double mass=0;
            double xVal = 0;
            double yVal = 0;
            double zVal = 0;
            int direction = 0;  
            
            const double G = 0.00000000006673;
            
            //create all three sizes of pixels
            
            for(int i=0; i<numParticlesLight; i++){
                //determining direction by computing x y and z components of velocity
                velocity=drand48()*(velocityLightMax-velocityLightMin)+velocityLightMin;

                direction = (int) (drand48() * 1.99);
                velocityX=drand48()*(velocity)*pow(-1,direction);
                
                
                direction = (int) (drand48() * 1.99);
                velocityY=drand48()*(sqrt(pow(velocity,2)-pow(velocityX,2)))*pow(-1,direction);

                direction = (int) (drand48() * 1.99);
                velocityZ=sqrt(pow(velocity,2)-pow(velocityX,2)-pow(velocityX,2))*pow(-1,direction);

                mass=drand48()*(massLightMax-massLightMin)+massLightMin;
                xVal = (drand48() * (width-1));
                yVal = (drand48() * (height-1));
                zVal = (drand48() * (zplane));
       
                particles[i] = vec3(xVal,yVal,zVal,255,0,0,mass,velocityX, velocityY, velocityZ); //vec3(1,1,1,255,255,255);
                
            }
            count=numParticlesLight;
            
            for(int j=count; j<(numParticleMedium+count); j++){
                velocity=drand48()*(velocityMediumMax-velocityMediumMin)+velocityMediumMin;

                direction = (int) (drand48() * 1.99);
                velocityX=drand48()*(velocity)*pow(-1,direction);
                
                
                direction = (int) (drand48() * 1.99);
                velocityY=drand48()*(sqrt(pow(velocity,2)-pow(velocityX,2)))*pow(-1,direction);

                direction = (int) (drand48() * 1.99);
                velocityZ=sqrt(pow(velocity,2)-pow(velocityX,2)-pow(velocityX,2))*pow(-1,direction);
                
                mass=drand48()*(massMediumMax-massMediumMin)+massMediumMin;
                xVal = (drand48() * (width-1));
                yVal = (drand48() * (height-1));
                zVal = (drand48() * (zplane));
            
                particles[j] = vec3(xVal,yVal,zVal,0,255,0,mass,velocityX, velocityY, velocityZ);//vec3(10,10,10,255,255,255);
                
            }
            count=numParticlesLight+numParticleMedium;
            
            for(int k=count; k<(numParticleHeavy+count); k++){
                velocity=drand48()*(velocityHeavyMax-velocityHeavyMin)+velocityHeavyMin;

                direction = (int) (drand48() * 1.99);
                velocityX=drand48()*(velocity)*pow(-1,direction);
                
                
                direction = (int) (drand48() * 1.99);
                velocityY=drand48()*(sqrt(pow(velocity,2)-pow(velocityX,2)))*pow(-1,direction);

                direction = (int) (drand48() * 1.99);
                velocityZ=sqrt(pow(velocity,2)-pow(velocityX,2)-pow(velocityX,2))*pow(-1,direction);

                mass=drand48()*(massHeavyMax-massHeavyMin)+massHeavyMin;
                xVal = (drand48() * (width-1));
                yVal = (drand48() * (height-1));
                zVal = (drand48() * (zplane));
             
                particles[k] = vec3(xVal,yVal,zVal,0,0,255,mass,velocityX, velocityY, velocityZ);//vec3(100,100,100,255,255,255);
                
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
            
       
            const unsigned char* result = (image);
            saveBMP (filename, result, width, height);
            
            
            
            
            
            /* =====================================================================
             
             Computing velocitites and positions for next frame
             
             ===================================================================== */
            
            
            
            
            
            //double outerDotMag = 0;
            //double forcesArr = (double *)malloc(totalParticles*totalParticles);
            
            double * forcesArr;
            forcesArr = (double *)malloc(3*totalParticles*totalParticles*sizeof(double));
            
            
            //memset(forcesArr, 0, sizeof(double) * totalParticles*totalParticles);
            for(int i = 0; i < (3*totalParticles*totalParticles) ; i++){
                forcesArr[i] = 0;
            }

            
            double currentForceX = 0;
            double currentForceY = 0;
            double currentForceZ = 0;
            double magX=0;
            double magY=0;
            double magZ=0;

            double forcesX = 0;
            double forcesY = 0;
            double forcesZ = 0;
            double totForceX =  0;
            double totForceY =  0;
            double totForceZ =  0;
            double currentMass = 0;
            
            // STEP
            for(int step = 0; step < numSteps; step++){
                
                //SUBSTEP
                for(int subStep = 0; subStep < numSubSteps; subStep++){

                    start=MPI_Wtime();

                    //CALCULATING FORCES AT EACH SUBSTEP FOR EACH PARTICLE
                    for(int outerDot = 0; outerDot < totalParticles; outerDot++){
                        currentForceX = 0;
                        currentForceY = 0;
                        currentForceZ = 0;
                        magX=0;
                        magY=0;
                        magZ=0;
                        
                        
                        //SUMMING THE FORCES OF INNERDOT PARTICLES ACTING UPON OUTERDOT PARTICLES
                        for(int innerDot = outerDot+1; innerDot < totalParticles; innerDot++){
                            magX=pow(particles[outerDot].getX() - particles[innerDot].getX(),2);
                            magY=pow(particles[outerDot].getY() - particles[innerDot].getY(),2);
                            magZ=pow(particles[outerDot].getZ() - particles[innerDot].getZ(),2);
                            
                            currentForceX = particles[innerDot].getMass()*(particles[outerDot].getX()-particles[innerDot].getX())/(pow(sqrt(magX+magY+magZ),3));
                            currentForceY = particles[innerDot].getMass()*(particles[outerDot].getY()-particles[innerDot].getY())/(pow(sqrt(magX+magY+magZ),3));
                            currentForceZ = particles[innerDot].getMass()*(particles[outerDot].getZ()-particles[innerDot].getZ())/(pow(sqrt(magX+magY+magZ),3));
                            
                            if(!isnan(currentForceX)){
                              
                                forcesArr[innerDot*3 + totalParticles*outerDot*3] = currentForceX;
                                forcesArr[outerDot*3 + totalParticles*innerDot*3] = -1 * currentForceX;
                            }
                            if(!isnan(currentForceY)){
                                
                                forcesArr[innerDot*3 + totalParticles*outerDot*3 + 1] = currentForceY;
                                forcesArr[outerDot*3 + totalParticles*innerDot*3 + 1] = -1 * currentForceY;
                            }
                            if(!isnan(currentForceZ)){
                                
                                forcesArr[innerDot*3 + totalParticles*outerDot*3 + 2] = currentForceZ;
                                forcesArr[outerDot*3 + totalParticles*innerDot*3 + 2] = -1 * currentForceZ;
                            }
                        }
                        
                    }
                  
                    //RECALCULATING POSITION AND VELOCITY FOR EACH PARTICLE AT EACH SUBSTEP
                    forcesX = 0;
                    forcesY = 0;
                    forcesZ = 0;
                    totForceX =  0;
                    totForceY =  0;
                    totForceZ =  0;
                    currentMass = 0;
                    for(int index = 0; index < totalParticles; index++){
                        currentMass = particles[index].getMass();
                        totForceX =  0;
                        totForceY =  0;
                        totForceZ =  0;
                        forcesX = 0;
                        forcesY = 0;
                        forcesZ = 0;
                        for(int length = 0; length < totalParticles; length ++){
                           
                            forcesX += forcesArr[length*3 + index*totalParticles*3];
                            forcesY += forcesArr[length*3 + index*totalParticles*3 + 1];
                            forcesZ += forcesArr[length*3 + index*totalParticles*3 + 2];
                 
                        }
                        
                        totForceX = (-1)*G*forcesX*currentMass;
                        totForceY = (-1)*G*forcesY*currentMass;
                        totForceZ = (-1)*G*forcesZ*currentMass;
                        
                        if(totForceX < epsilon && totForceX>0){
                            totForceX = epsilon;
                        }
                        if(totForceX > epsilon && totForceX<0){
                            totForceX = (-1)*epsilon;
                        }
                        if(totForceY < epsilon && totForceY>0){
                            totForceY = epsilon;
                        }
                        if(totForceY > epsilon && totForceY<0){
                            totForceY = (-1)*epsilon;
                        }
                        if(totForceZ < epsilon && totForceZ>0){
                            totForceZ = epsilon;
                        }
                        if(totForceZ > epsilon && totForceZ<0){
                            totForceZ = (-1)*epsilon;
                        }
                        
                        
                        particles[index].setPosition( particles[index], timeSubStep );
                        
                        particles[index].setVelocityX( particles[index], timeSubStep, totForceX );
                        particles[index].setVelocityY( particles[index], timeSubStep, totForceY );
                        particles[index].setVelocityZ( particles[index], timeSubStep, totForceZ );
                    }
//////////
                    end=MPI_Wtime();
                    timeRequired=end-start;
                    subStepTimes[step*numSubSteps+subStep]=timeRequired;
                    start=0;
                    end=0;
                    timeRequired=0;
                }
                
                //ALL PIXELS BLACK
                for(int a=0; a<(3*frameSize);a++){
                    image[a]=(unsigned char) 0; //all pixels return to black
                }
                
                //WRITE PIXELS TO IMAGE FOR THIS FRAME
                for(int g=0;g<totalParticles;g++){
                    int xValue=(int) particles[g].getX();
                    int yValue=(int) particles[g].getY();
                    if(!(particles[g].getX() > width || particles[g].getX() < 0 || particles[g].getY() > height || particles[g].getY() < 0)){
                        image[(xValue*3)+(yValue*width*3)]=(unsigned char) particles[g].getR();
                        image[(xValue*3)+(yValue*width*3)+1]=(unsigned char) particles[g].getG();
                        image[(xValue*3)+(yValue*width*3)+2]=(unsigned char) particles[g].getB();
                    }
                    
                }
                
                
                //CORRECTLY FORMAT NAME OF FILE (number of zeros)
                int numDigits = 0;
                int zeroChecker = step+1;
                while(zeroChecker >= 1){
                    zeroChecker /= 10;
                    numDigits++;
                }

                std::string two = "";
                switch(numDigits){
                    case 5:
                        two += "_";
                        break;
                    case 4:
                        two += "_0";
                        break;
                    case 3:
                        two += "_00";
                        break;
                    case 2:
                        two += "_000";
                        break;
                    case 1:
                        two += "_0000";
                        break;
                    default:
                        two += "_00000";
                        break;
                }
              
                //GENERATE FRAME NAME FOR IMAGE
                std::string one=argv[9];
                std::string three= std::to_string(step+1);
                std::string four=".bmp";
                
                std::string resultat= one+two+three+four;
                
                const char* filename1 = resultat.c_str();
                
                const unsigned char* result1 = (image);
                saveBMP (filename1, result1, width, height);
              
                
            }

   
            free(forcesArr);
            free(image);
        
            double minTime=1000;
            double maxTime=0;
            double countTime=0;
            for(int i=0;i<(numSteps*numSubSteps);i++){
                countTime+=subStepTimes[i];
                
                if(subStepTimes[i]<minTime){
                    minTime=subStepTimes[i];
                }
                if(subStepTimes[i]>maxTime){
                    maxTime=subStepTimes[i];
                }
            }

            countTime=countTime/(numSteps*numSubSteps);
            
            printf("min time of substeps: %f \n",minTime);
            printf("max time of substeps: %f \n",maxTime);
            printf("average time of substeps: %f \n", countTime);

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