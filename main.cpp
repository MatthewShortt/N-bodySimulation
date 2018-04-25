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
        
        
        int glob_steps = atoi(argv[4]) * atoi(argv[5]);;
        
        
        
        int numSteps = 0;
        int numSubSteps = 0;
        double timeSubStep;
        
        int width, height;
        
        int totalParticles = (atoi(argv[1]) + atoi(argv[2]) + atoi(argv[3]));
        
        unsigned char* image;
        
        
        double * forcesArr = (double *)malloc(3*totalParticles*totalParticles*sizeof(double));
        
        
        //int length_forcesArr = 3*totalParticles*totalParticles;
        //double buff[length_forcesArr];
        
        for(int i = 0; i < (3*totalParticles*totalParticles) ; i++){
            forcesArr[i] = 0;
        }
        
        vec3 particles [totalParticles];
        
        vec3 tempParticles [totalParticles];
        
        MPI_Status status;
        
        int blockSizeP, startIndexForP;
        
        double glob_timeSubStep = (double) atof(argv[6]);
        
        int num_members = 10;
        int lengths [num_members];
        MPI_Datatype types [num_members];
        
        MPI_Aint offsets [num_members];
        offsets [0]=offsetof (vec3, x);
        offsets [1]=offsetof (vec3, y);
        offsets [2]=offsetof (vec3, z);
        offsets [3]=offsetof (vec3, r);
        offsets [4]=offsetof (vec3, g);
        offsets [5]=offsetof (vec3, b);
        offsets [6]=offsetof (vec3, mass);
        offsets [7]=offsetof (vec3, velocityX);
        offsets [8]=offsetof (vec3, velocityY);
        offsets [9]=offsetof (vec3, velocityZ);
        
        for(int i = 0; i < num_members; i++){
            lengths[i] = 1;
            if(i == 3 || i == 4 || i == 5){
                types[i] = MPI_INT;
            } else{
                types[i] = MPI_DOUBLE;
            }
        }
        //int lengths [num_members] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        //= {offsetof (vec3, x), offsetof(vec3, y), offsetof (vec3, z)};
        MPI_Datatype MPI_VEC3;
        MPI_Type_create_struct (num_members, lengths, offsets, types, &MPI_VEC3);
        MPI_Type_commit(&MPI_VEC3);
        
        const double G = 0.00000000006673;
        
        //root node stuff goes here
        if(my_rank == 0){
            
            //grab arguments
            numParticlesLight = atoi(argv[1]);
            numParticleMedium = atoi(argv[2]);
            numParticleHeavy  = atoi(argv[3]);
            
            numSteps = atoi(argv[4]);
            numSubSteps = atoi(argv[5]);
            
            timeSubStep = (double) atof(argv[6]);
            
            timeSubStep = timeSubStep + 1 - 1;//for testing
            
            width = atoi(argv[7]);
            height = atoi(argv[8]);
            
            //depth of z plane
            const int zplane=100;
            
            
            /* =====================================================================
             
             Setup
             
             ===================================================================== */
            
            totalParticles=numParticlesLight+numParticleMedium+numParticleHeavy;
            
            int frameSize=width*height;
            
            //allocate space for pixel colours
            image = (unsigned char *)malloc(3*frameSize*sizeof(unsigned char));
            
            //default all pixels to white/black
            for(int a=0; a<(3*frameSize);a++){
                image[a]=(unsigned char) 0;
                
            }
            
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
            
//            double currentForceX = 0;
//            double currentForceY = 0;
//            double currentForceZ = 0;
//            double magX=0;
//            double magY=0;
//            double magZ=0;
//            
//            double forcesX = 0;
//            double forcesY = 0;
//            double forcesZ = 0;
//            double totForceX =  0;
//            double totForceY =  0;
//            double totForceZ =  0;
//            double currentMass = 0;
            
            
            //Matts vars
            int numRows        = (totalParticles)/(p-1);
            int numExtraRows   = (totalParticles)%(p-1);
            int startIndexForP = 0;
            blockSizeP     = 0;
            int startIncrement = 0;
            int rowSize        = 3;
        
            int tempStartIndex1 = 0;
            int tempNumParticles1 = 0;
            
            // STEP
            for(int step = 0; step < numSteps; step++){
                
                //SUBSTEP
                for(int subStep = 0; subStep < numSubSteps; subStep++){
                    

                    start=MPI_Wtime();
                    startIncrement = 0;

                    memcpy(tempParticles, particles, totalParticles*sizeof(vec3));
                    
                    for(int q = 1; q < p; q++){
                        startIndexForP = startIncrement;
                        
                        if(q <= numExtraRows){
                            blockSizeP = (numRows + 1)*rowSize;
                        } else{
                            blockSizeP = (numRows    )*rowSize;
                        }

                        startIncrement += blockSizeP;
                        MPI_Send(&startIndexForP, 1             , MPI_INT , q, 1, MPI_COMM_WORLD);
                        MPI_Send(&blockSizeP    , 1             , MPI_INT , q, 1, MPI_COMM_WORLD);
                        MPI_Send(&tempParticles     , totalParticles, MPI_VEC3, q, 1, MPI_COMM_WORLD);
                        
                    }
                    
                    for(int src = 1; src < p; src++){
                        
                        MPI_Recv(&startIndexForP, 1, MPI_INT, src, 2, MPI_COMM_WORLD, &status);
                        MPI_Recv(&blockSizeP    , 1, MPI_INT, src, 2, MPI_COMM_WORLD, &status);
                        MPI_Recv(&tempParticles    , totalParticles, MPI_VEC3, src, 2, MPI_COMM_WORLD, &status);
                        
                        tempNumParticles1 = blockSizeP / 3;
                        tempStartIndex1 = startIndexForP / 3;
                        
                        for(int i = tempStartIndex1; i < (tempNumParticles1+tempStartIndex1); i++){
                            particles[i] = tempParticles[i-tempStartIndex1];
                        }
  
                    }
                    
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
            
            int tempNumParticles = 0;
            int tempStartIndex = 0;
            
            int counter2 = 0;
            
            
            double posX, posY, posZ, velX, velY, velZ;
            
            
            for(int cnt = 0; cnt < glob_steps; cnt++){
                MPI_Recv(&startIndexForP, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(&blockSizeP    , 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(&particles     , totalParticles, MPI_VEC3, 0, 1, MPI_COMM_WORLD, &status);

                tempNumParticles = blockSizeP / 3;
                tempStartIndex = startIndexForP / 3;
                
                double * tempForces = (double *)malloc(blockSizeP*totalParticles*sizeof(double));

                
                //CALCULATING FORCES AT EACH SUBSTEP FOR EACH PARTICLE
                for(int outerDot = tempStartIndex; outerDot < (tempStartIndex + tempNumParticles); outerDot++){
                    currentForceX = 0;
                    currentForceY = 0;
                    currentForceZ = 0;
                    magX=0;
                    magY=0;
                    magZ=0;
                    
                    
                    //SUMMING THE FORCES OF INNERDOT PARTICLES ACTING UPON OUTERDOT PARTICLES
                    for(int innerDot = 0; innerDot < totalParticles; innerDot++){
                        
                        if(innerDot != outerDot){
                            magX=pow(particles[outerDot].getX() - particles[innerDot].getX(),2);
                            magY=pow(particles[outerDot].getY() - particles[innerDot].getY(),2);
                            magZ=pow(particles[outerDot].getZ() - particles[innerDot].getZ(),2);
                            
                            currentForceX = particles[innerDot].getMass()*(particles[outerDot].getX()-particles[innerDot].getX())/(pow(sqrt(magX+magY+magZ),3));
                            currentForceY = particles[innerDot].getMass()*(particles[outerDot].getY()-particles[innerDot].getY())/(pow(sqrt(magX+magY+magZ),3));
                            currentForceZ = particles[innerDot].getMass()*(particles[outerDot].getZ()-particles[innerDot].getZ())/(pow(sqrt(magX+magY+magZ),3));
                        }

                       
                        if(!isnan(currentForceX)){
                            
                            tempForces[(innerDot)*3 + totalParticles*(outerDot-tempStartIndex)*3] = currentForceX;
                        }
                        if(!isnan(currentForceY)){
                            
                            tempForces[(innerDot)*3 + totalParticles*(outerDot-tempStartIndex)*3 + 1] = currentForceY;
                        }
                        if(!isnan(currentForceZ)){
                            
                            tempForces[(innerDot)*3 + totalParticles*(outerDot-tempStartIndex)*3 + 2] = currentForceZ;
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
                
                vec3 tempVec;
                
                
                for(int index = tempStartIndex; index < (tempNumParticles+tempStartIndex); index++){
                    
                    currentMass = particles[index].getMass();
                    totForceX =  0;
                    totForceY =  0;
                    totForceZ =  0;
                    forcesX = 0;
                    forcesY = 0;
                    forcesZ = 0;
                    for(int length = 0; length < totalParticles; length ++){
                        
                        forcesX += tempForces[length*3 + (index-tempStartIndex)*totalParticles*3];
                        forcesY += tempForces[length*3 + (index-tempStartIndex)*totalParticles*3 + 1];
                        forcesZ += tempForces[length*3 + (index-tempStartIndex)*totalParticles*3 + 2];
                        
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
                    
                    
                    posX = particles[index].getX() + glob_timeSubStep * particles[index].getVelocityX();
                    posY = particles[index].getY() + glob_timeSubStep * particles[index].getVelocityY();
                    posZ = particles[index].getZ() + glob_timeSubStep * particles[index].getVelocityZ();

                    velX = particles[index].getVelocityX() + ((glob_timeSubStep*totForceX)/ particles[index].getMass());
                    velY = particles[index].getVelocityY() + ((glob_timeSubStep*totForceY)/ particles[index].getMass());
                    velZ = particles[index].getVelocityZ() + ((glob_timeSubStep*totForceZ)/ particles[index].getMass());
                    
                    tempVec = vec3(posX,posY,posZ,particles[index].getR(),particles[index].getG(),particles[index].getB(),particles[index].getMass(),velX, velY, velZ);
                    particles[index] = tempVec;
                    
                    counter2 ++;
                }
                
                free(tempForces);
                
                
                
                MPI_Send(&startIndexForP, 1             , MPI_INT , 0, 2, MPI_COMM_WORLD);
                MPI_Send(&blockSizeP    , 1             , MPI_INT , 0, 2, MPI_COMM_WORLD);
                MPI_Send(&particles     , totalParticles, MPI_VEC3, 0, 2, MPI_COMM_WORLD);

            }
            
        }
        
        
        MPI_Finalize();
        
    } else{
        printf("\n--> Program did not run do to input error. Please view message above and fix input arguments.\n");
    }
    
    
    
    return 0;
}