#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <mpi.h>
#include <math.h>
#include <cmath>

#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"

#include <cstdint>


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
        //double start=0;
        //double end=0;
        //double timeRequired=0;
        //double subStepTimes [atoi(argv[4])*atoi(argv[5])];
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
        
        
        unsigned char image[3*atoi(argv[7])*atoi(argv[8])];
        
        MPI_Status status;
        
        
        //variables that need to be declared for MPI
        int frameSize, blockSizeP, startIndexForP;
        
        
        //unsigned char * buff[128];
        
        
        
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
            //const int zplane=100;
            
            //print args to the console
            printf("numParticlesLight: %d\nnumParticlesMedium: %d\nnumParticlesLarge: %d\nnumSteps: %d\nsubSteps: %d\ntimeSubStep: %f\nWidth: %d\nHeight: %d\n",numParticlesLight,numParticleMedium,numParticleHeavy,numSteps,numSubSteps,timeSubStep, width,height);
            
            
            /* =====================================================================
             
             Setup
             
             ===================================================================== */
            
            //const int totalParticles=numParticlesLight+numParticleMedium+numParticleHeavy;
            
            frameSize=width*height;
            
            //allocate space for pixel colours
            //image = (unsigned char *)malloc(3*frameSize*sizeof(unsigned char));
            
            //            for(int a=0; a<(3*frameSize);a++){//a+=3
            //                image[a]=(unsigned char) 123;//0;
            //                //image[a+1]=(unsigned char) 0;
            //                //image[a+2]=(unsigned char) 0;
            //            }
            
            
            //printf("address: %p\n", (void *) &image);
            //default all pixels to white/black
            
            //memset(image, 0, sizeof(unsigned char) * 3 * frameSize);
            
            
            
            /*   Practice MPI start  */
            
            int msg_type = 1;
            blockSizeP = (3 * width * height) / ((p-1)*2); //block size, ie the number of values each processor will calculate
            printf("Practice block size: %d\n\n", blockSizeP);


            
                
                
                startIndexForP = 0;
                msg_type = 1; //from processor 0 to another processor
                for(int q = 1; q < p; q++){
                    
                    printf("Made it inside for...\n");
                    
                    startIndexForP = (q-1)*blockSizeP + i*(p-1)*blockSizeP; //calculating the start index for each processor
                    
                    /*       What ur sending, size of that ,  type   , p, from who, over channel   */
                    MPI_Send(&startIndexForP       , 1         , MPI_INT          , q, msg_type, MPI_COMM_WORLD);
                    MPI_Send(&blockSizeP           , 1         , MPI_INT          , q, msg_type, MPI_COMM_WORLD);
                    MPI_Send(&frameSize            , 1         , MPI_INT          , q, msg_type, MPI_COMM_WORLD);
                    MPI_Send(&image[startIndexForP], blockSizeP, MPI_UNSIGNED_CHAR, q, msg_type, MPI_COMM_WORLD);
                    
                    
                    
                }
                
                printf("Outside for!\n");
                
                /** ==============================
                 -- Receive Results form Processors
                 ============================== **/
                
                
                
                for(int src = 1; src < p; src++){
                    printf("Made it inside for #2...\n");
                    MPI_Recv(&startIndexForP   , 1         , MPI_INT          , src, 2, MPI_COMM_WORLD, &status);
                    MPI_Recv(&blockSizeP       , 1         , MPI_INT          , src, 2, MPI_COMM_WORLD, &status);
                    MPI_Recv(&frameSize        , 1         , MPI_INT          , src, 2, MPI_COMM_WORLD, &status);
                    MPI_Recv(&image[startIndexForP], blockSizeP, MPI_UNSIGNED_CHAR, src, 2, MPI_COMM_WORLD, &status);
                    printf("just after recvs...\n");
                    //                for(int i = 0; i < (3*frameSize); i++){
                    //                    printf("Value at %d: %d\n", i, (int)image[i]);
                    //                }
                    //                printf("\n");
                    
                }
                




                
                printf("DONE!!!!\n");
                
            

            
            
            for(int i = 0; i < (3*frameSize); i++){
                printf("Value at %d: %d\n", i, (int)image[i]);
            }
            

            
            
            
            
            
            
            

            
            char file[20];
            strcpy(file, argv[9]);
            strcat(file,"_00000.bmp");
            const char* filename = file;
            
            
            //printf("this is arg 9: %s %s %s \n", argv[9],argv[8], argv[7]);
            const unsigned char* result = (image);
            saveBMP (filename, result, width, height);
            
            /*   Practice MPI end  */
            
            /* Repaste code here start */
            
            /* Repaste code here end */
            
            //free(image);
            
        }
        //all other nodes do this
        else{
            
            int msg_type = 1;
            //for(int i = 0;  i < 2; i++){
            msg_type = 1; //from process 0
            
            MPI_Recv(&startIndexForP, 1         , MPI_INT          , 0, msg_type, MPI_COMM_WORLD, &status);
            MPI_Recv(&blockSizeP    , 1         , MPI_INT          , 0, msg_type, MPI_COMM_WORLD, &status);
            MPI_Recv(&frameSize     , 1         , MPI_INT          , 0, msg_type, MPI_COMM_WORLD, &status);
            MPI_Recv(&image         , blockSizeP, MPI_UNSIGNED_CHAR, 0, msg_type, MPI_COMM_WORLD, &status);
            
            //printf("I am rank: %d with start index: %d and blockSize: %d  and startIndexForP + blockSizeP: %d \t\t", my_rank, startIndexForP, blockSizeP, startIndexForP + blockSizeP);
            
            unsigned char my_colour = (unsigned char) (my_rank * 25);
            
            
            
            int temper = (int) my_colour;
            //printf("my colour: %d\n", temper);
            
            for(int j=0; j < 2; j++){
                
                for(int indexP = 0; indexP < blockSizeP; indexP++){
                    //printf("Assigning value: %d to index %d\n", temper, indexP);
                    image[indexP] = (unsigned char) temper;
                }
                
                msg_type = 2; //from proessor that isn't 0
                
                MPI_Send(&startIndexForP, 1         , MPI_INT          , 0, msg_type, MPI_COMM_WORLD);
                MPI_Send(&blockSizeP    , 1         , MPI_INT          , 0, msg_type, MPI_COMM_WORLD);
                MPI_Send(&frameSize     , 1         , MPI_INT          , 0, msg_type, MPI_COMM_WORLD);
                MPI_Send(&image, blockSizeP, MPI_UNSIGNED_CHAR, 0, msg_type, MPI_COMM_WORLD);
            }
            
            //}
            
            
            
        }
        
        
        MPI_Finalize();
        
    } else{
        printf("\n--> Program did not run do to input error. Please view message above and fix input arguments.\n");
    }
    
    
    
    return 0;
}





//            for(int a=0; a<(3*frameSize);a++){//a+=3
//                image[a]=(unsigned char) 0;//0;
//                //image[a+1]=(unsigned char) 0;
//                //image[a+2]=(unsigned char) 0;
//            }
//
//            vec3 particles [totalParticles];
//            int count=0;
//            double velocity=0;
//            double mass=0;
//            double xVal = 0;
//            double yVal = 0;
//            double zVal = 0;
//            int direction = 0;
//
//            /*  === Directions Layout ===
//
//             x       y       z     dir (int)
//             1       0       0  --> 0
//             0       1       0  --> 1
//             0       0       1  --> 2
//             1       1       0  --> 3
//             1       0       1  --> 4
//             0       1       1  --> 5
//             1       1       1  --> 6
//
//             ================ */
//
//
//
//            //const double G = -0.00000000006673;
//
//            //create all three sizes of pixels
//
//            for(int i=0; i<numParticlesLight; i++){
//                velocity=drand48()*(velocityLightMax-velocityLightMin)+velocityLightMin;
//                mass=drand48()*(massLightMax-massLightMin)+massLightMin;
//                xVal = (drand48() * (width-1));
//                yVal = (drand48() * (height-1));
//                zVal = (drand48() * (zplane));
//                direction = (int) (drand48() * 6.99);
//                //printf("Coordinates X,Y,Z: %d %d %d\n", xVal, yVal, zVal);
//                particles[i] = vec3(xVal,yVal,zVal,255,0,0,mass,velocity,direction); //vec3(1,1,1,255,255,255);
//
//            }
//            count=numParticlesLight;
//
//            for(int j=count; j<(numParticleMedium+count); j++){
//                velocity=drand48()*(velocityMediumMax-velocityMediumMin)+velocityMediumMin;
//                mass=drand48()*(massMediumMax-massMediumMin)+massMediumMin;
//                xVal = (drand48() * (width-1));
//                yVal = (drand48() * (height-1));
//                zVal = (drand48() * (zplane));
//                direction = (int) (drand48() * 6.99);
//                particles[j] = vec3(xVal,yVal,zVal,0,255,0,mass,velocity,direction);//vec3(10,10,10,255,255,255);
//
//            }
//            count=numParticlesLight+numParticleMedium;
//
//            for(int k=count; k<(numParticleHeavy+count); k++){
//                velocity=drand48()*(velocityHeavyMax-velocityHeavyMin)+velocityHeavyMin;
//                mass=drand48()*(massHeavyMax-massHeavyMin)+massHeavyMin;
//                xVal = (drand48() * (width-1));
//                yVal = (drand48() * (height-1));
//                zVal = (drand48() * (zplane));
//                direction = (int) (drand48() * 6.99);
//                particles[k] = vec3(xVal,yVal,zVal,0,0,255,mass,velocity,direction);//vec3(100,100,100,255,255,255);
//
//            }
//
//
//            //write pixels onto screen
//            for(int g=0;g<totalParticles;g++){
//                int xValue=(int) particles[g].getX();
//                int yValue=(int) particles[g].getY();
//                image[(xValue*3)+(yValue*width*3)]=(unsigned char) particles[g].getR();
//                image[(xValue*3)+(yValue*width*3)+1]=(unsigned char) particles[g].getG();
//                image[(xValue*3)+(yValue*width*3)+2]=(unsigned char) particles[g].getB();
//
//            }
//
//            char file[20];
//            strcpy(file, argv[9]);
//            strcat(file,"_00000.bmp");
//            const char* filename = file;
//
//
//            //printf("this is arg 9: %s %s %s \n", argv[9],argv[8], argv[7]);
//            const unsigned char* result = (image);
//            saveBMP (filename, result, width, height);
//
//
//
//
//
//            /* =====================================================================
//
//             Computing velocitites and positions for next frame
//
//             ===================================================================== */
//
//
//
//
//
//            double outerDotMagSquared = 0;
//            //double forcesArr = (double *)malloc(totalParticles*totalParticles);
//
//            double * forcesArr;
//            forcesArr = (double *)malloc(totalParticles*totalParticles*sizeof(double));
//
//
//            //memset(forcesArr, 0, sizeof(double) * totalParticles*totalParticles);
//            for(int i = 0; i < (totalParticles*totalParticles) ; i++){
//                forcesArr[i] = 0;
//            }
//
//            double currentForce = 0;
//
//            double forces = 0;
//            double totForce =  0;
//            double currentMass = 0;
//
//            // STEP
//            for(int step = 0; step < numSteps; step++){
//
//                //SUBSTEP
//                for(int subStep = 0; subStep < numSubSteps; subStep++){
//
//                    start=MPI_Wtime();
//
//                    //CALCULATING FORCES AT EACH SUBSTEP FOR EACH PARTICLE
//                    for(int outerDot = 0; outerDot < totalParticles; outerDot++){
//                        currentForce = 0;
//                        outerDotMagSquared = particles[outerDot].MagnitudeSquared();
//
//                        //SUMMING THE FORCES OF INNERDOT PARTICLES ACTING UPON OUTERDOT PARTICLES
//                        for(int innerDot = outerDot+1; innerDot < totalParticles; innerDot++){
//
//                            currentForce = particles[innerDot].getMass()/(abs(outerDotMagSquared - particles[innerDot].MagnitudeSquared()));
//
//                            forcesArr[innerDot + totalParticles*outerDot] = currentForce;
//                            forcesArr[outerDot + totalParticles*innerDot] = -1 * currentForce;
//
//                        }
//
//                    }
//
//                    //RECALCULATING POSITION AND VELOCITY FOR EACH PARTICLE AT EACH SUBSTEP
//                    forces = 0;
//                    totForce =  0;
//                    currentMass = 0;
//                    for(int index = 0; index < totalParticles; index++){
//                        currentMass = particles[index].getMass();
//                        totForce = 0;
//                        forces = 0;
//                        for(int length = 0; length < totalParticles; length ++){
//                            forces += forcesArr[length + index*totalParticles];
//                        }
//
//                        totForce = forces*currentMass;
//                        if(totForce < epsilon){
//                            totForce = epsilon;
//                        }
//
//                        particles[index].setPosition( particles[index], timeSubStep );
//
//                        particles[index].setVelocity( particles[index], timeSubStep, totForce );
//                    }
//
//                    end=MPI_Wtime();
//                    timeRequired=end-start;
//                    subStepTimes[step*numSubSteps+subStep]=timeRequired;
//                    start=0;
//                    end=0;
//                    timeRequired=0;
//                }
//
//                //ALL PIXELS BLACK
//                for(int a=0; a<(3*frameSize);a++){
//                    image[a]=(unsigned char) 0; //all pixels return to black
//                }
//
//                //WRITE PIXELS TO IMAGE FOR THIS FRAME
//                for(int g=0;g<totalParticles;g++){
//                    int xValue=(int) particles[g].getX();
//                    int yValue=(int) particles[g].getY();
//                    if(!(particles[g].getX() > width || particles[g].getX() < 0 || particles[g].getY() > height || particles[g].getY() < 0)){
//                        image[(xValue*3)+(yValue*width*3)]=(unsigned char) particles[g].getR();
//                        image[(xValue*3)+(yValue*width*3)+1]=(unsigned char) particles[g].getG();
//                        image[(xValue*3)+(yValue*width*3)+2]=(unsigned char) particles[g].getB();
//                    }
//
//                }
//
//
//                //CORRECTLY FORMAT NAME OF FILE (number of zeros)
//                int numDigits = 0;
//                int zeroChecker = step+1;
//                while(zeroChecker >= 1){
//                    zeroChecker /= 10;
//                    numDigits++;
//                }
//
//                std::string two = "";
//                switch(numDigits){
//                    case 5:
//                        two += "_";
//                        break;
//                    case 4:
//                        two += "_0";
//                        break;
//                    case 3:
//                        two += "_00";
//                        break;
//                    case 2:
//                        two += "_000";
//                        break;
//                    case 1:
//                        two += "_0000";
//                        break;
//                    default:
//                        two += "_00000";
//                        break;
//                }
//
//                //GENERATE FRAME NAME FOR IMAGE
//                std::string one=argv[9];
//                std::string three= std::to_string(step+1);
//                std::string four=".bmp";
//
//                std::string resultat= one+two+three+four;
//
//                const char* filename1 = resultat.c_str();
//
//                const unsigned char* result1 = (image);
//                saveBMP (filename1, result1, width, height);
//
//
//            }
//
//            free(forcesArr);
//            free(image);
//
//        double minTime=1000;
//        double maxTime=0;
//        double countTime=0;
//        for(int i=0;i<(numSteps*numSubSteps);i++){
//            countTime+=subStepTimes[i];
//            printf("this is time: %f \n",subStepTimes[i]);
//            if(subStepTimes[i]<minTime){
//                minTime=subStepTimes[i];
//            }
//            if(subStepTimes[i]>maxTime){
//                maxTime=subStepTimes[i];
//            }
//        }
//
//        countTime=countTime/(numSteps*numSubSteps);
//
//        printf("min time of substeps: %f \n",minTime);
//        printf("max time of substeps: %f \n",maxTime);
//        printf("average time of substeps: %f \n", countTime);
