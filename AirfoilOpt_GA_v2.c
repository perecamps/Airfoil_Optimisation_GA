// FINAL MASTER THESIS: Master's degree in Aerospace Engineering
// TITLE:               Genetic Algorithm for Airfoil Optimisation
// CONTACT:             Pere Camps Castellanos - pere.camps.c@gmail.com - linkedin.com/in/pere-camps-castellanos-35678b143/
// UNIVERSITY:          UPC, Escola Superior d'Enginyeries Industrial, Aeroespacial i Audiovisual de Terrassa (ESEIAAT)



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ABSTRACT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this Master Thesis, an airfoil optimisation using a Genetic Algorithm is developed.  This project has been proposed by myself and done with the guidance and supervision of professor Manel Soria.

The main idea of the project is to develop from scratch an algorithm capable of finding the optimal airfoil for specific flow conditions, such as the angle of attack, the Reynolds number, and the Mach number. 
The objective is to create a useful tool for aerospace engineering students so they can use it on their projects and designs during the college years.\\

The work has a first theoretical part about Genetic Algorithms, in which the basic concepts needed to understand the current project are explained. Then, the implementation of the algorithm is fully expounded 
and all the intern processes of the genetic algorithm can be consulted. Several validations of the code have also been made.\\

The Genetic Algorithm created uses crossovers and mutations. The airfoil parametrisation used has been the PARSEC parametrisation and the computation of the aerodynamic coefficients is done with XFOIL. 
The whole code is written in C language and the results analysis and graphs are done with MATLAB and XFLR5.\\

Finally, the algorithm is tested with two real design cases, an airfoil for a heavy lifter aircraft that participated in the Air Cargo Challenge 2017 in Stuttgart, and an airfoil for a glider that flew in the 
Paper Air Challenge 2015 in ESEIAAT, Terrasa. The results and improvements offered by the algorithm are compared with the results that the designers of these aircraft obtained manually during the design process.

*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

/* CONSTANTS AND DEFINITIONS */

// General Constants - DO NOT MODIFY
#define PI 3.14159265359             // PI
#define N 100                        // airfoil discretization  

// Genetic Algorithm Parameters
#define Generation_NOT_IMPROVE 8     // Maximum number of generations without improvement allowed before stop
#define MAX_INFEASIBLE         9     // Maximum number of infeasible generations in a row allowed before stop
#define MAX_GENERATIONS        100   // Maximum number of generations allowed
#define MUTATION               1     // Mutation type [1, 2 or 3].
#define MUTATION_PROB          0.7   // Mutation Probability
#define individuals            160   // Initial Population

// Fitness function
#define FITNESS_TYPE 1       // 1 - Max Efficiency // 2 - Objective Cl and min Cd (define ClT) // 3 - Max Cl // 4 Min Cd
#define ClT 1.8              // Ovjecticve Cl used with FITNESS_TYPE 2

// CONSTRAINTS
#define MAX_THICKNESS 1      //  1 for no Constraint
#define MIN_THICKNESS 0      //  0 for no Constraint
#define TE_THIKNESS   0.000  //  0 for no Constrint
#define MIN_TE_ANGLE  1      //  1 for no Constraint
#define MIN_CM       -1      // -1 for no Constraint
#define EFF_LIM       200    // Efficiency limit to discard ridiculous solutions


// Flow Conditions - Modifyable
#define ALPHA     2            // Angle of attack [º]
#define REYNOLDS  1e6          // Re = rho*V*c/mu [-]
#define MACH      0.05         // v/a = v/sqrt(gamma*R*T) [-]

/* PARSEC LIMIT PARAMETRES - It is recommended NOT to modify them */
// R upper
#define rupi  0.001     
#define rups  0.2
// R lower
#define rloi  0.001     
#define rlos  0.15
// X upper
#define xupi  0.1       
#define xups  0.65
// Z upper
#define zupi  0.001      
#define zups  0.4
// X lower
#define xloi  0.1     
#define xlos  0.55
// Z lower
#define zloi -0.3       
#define zlos -0.01
// ZX upper
#define zxui  -1.8//-01.8     
#define zxus  -0.1 //-0.400
// ZX lower
#define zxli  0.1 //0.04      
#define zxls  1.2 //0.9
// Z trailing edge
#define ztei  0         
#define ztes  0.0
// dZ trailing edge
#define dzti  TE_THIKNESS         
#define dzts  TE_THIKNESS
// alpha
#define alpi  -32*PI/180    
#define alps  10*PI/180
// beta
#define beti  MIN_TE_ANGLE*PI/180      
#define bets  25*PI/180 


double limits_inferiors[12] = {rupi,rloi,xupi,zupi,xloi,zloi,zxui,zxli,ztei,dzti,alpi,beti};
double limits_superiors[12] = {rups,rlos,xups,zups,xlos,zlos,zxus,zxls,ztes,dzts,alps,bets};


/* FUNCTIONS */

// NOTE: Each function is explained below the main

// GENETIC
int      Genetic_Algorithm(double *PAR_dbl,double *newPAR_dbl, unsigned short *PAR_ush, unsigned short *newPAR_ush, double **sol_PAR_dbl, double *fits, double * Fits2Prob, double alpha, double Re, double Ma, double *generation_TOP_fit);
void     Fits_To_Probability(double *fits, double *Fits2Prob);
void     Natural_Selection (double *Fits2Prob, double *PAR_dbl, unsigned short *PAR_ush, double **C1, double **C2, unsigned short **C1ush, unsigned short **C2ush);
void     Crossover(double *newPAR_dbl, unsigned short *newPAR_ush, unsigned long *NewInd, double *C1, double *C2, unsigned short *C1ush, unsigned short *C2ush, int d);
void     OnePointCrossoverMutation(unsigned short p1,  unsigned short p2,unsigned short *f1, unsigned short *f2, float prob);

// PARSEC - XFOIL
void     Generate_Random_parameters(double *PAR_dbl, unsigned short *GA_2);
void     GenerateAirfoil(double *PAR_dbl, char a_name[], unsigned short *ALARM);    
double   callXfoil(double alpha, double Re, double Ma, double *PAR_dbl, double *Cl, double *Cd, double *Cm);
int      fexists(const char * filename);
unsigned long CountLines(char *filename, int header_lines); 

// RANDOM
double   Random_UD(double i, double j);
float    ran1(long *idum);
float    uniform(void);
void     randomize(void);
unsigned char random_bit(void);
unsigned short USHRTran(void);

// ERROR AND INFORMATION
void     ExitError(const char *miss, int errcode); 
void     PrintUserInfo(int i, int max);

// SYSTEM SOLVER
double   determinant(double a[6][6], double k);                                                                         
void     cofactor(double num[6][6], double f, double inv[6][6]);                                                        
void     transpose(double num[6][6], double fac[6][6], double r, double inverse[6][6]);                                 
void     multiplyMatrices(double first[6][6], double second[6][1], double mult[6][1], int r1, int c1, int r2, int c2);



/* MAIN */
int main(int argc, char *argv[]){

    randomize();

    double alpha = ALPHA, Re = REYNOLDS, Ma = MACH;
    unsigned short ga;

    double *fits = NULL, *Fits2Prob = NULL, *sol_PAR_dbl = NULL, *generation_TOP_fit = NULL;
    double *PAR_dbl = malloc(12*individuals*sizeof(double)), *newPAR_dbl = malloc(12*individuals*sizeof(double));
    unsigned short *PAR_ush = malloc(12*individuals*sizeof(int)),*newPAR_ush = malloc(12*individuals*sizeof(int));

    fits      = (double *) malloc(individuals*sizeof(double));
    Fits2Prob = (double *) malloc(individuals*sizeof(double));
    generation_TOP_fit = (double *) malloc(MAX_GENERATIONS*sizeof(double));  

    char *command = "start cmd @cmd /k \"kill_XFOIL\" ";   
    system(command);                                        // Execute the parallel code that kills XFOIL in case of stuck

    ga = Genetic_Algorithm(PAR_dbl,newPAR_dbl,PAR_ush,newPAR_ush, &sol_PAR_dbl, fits, Fits2Prob, alpha, Re, Ma, generation_TOP_fit);
    if(ga == 1){
        printf("Problem Solved\n");
    }
    printf("GENERATION TOP FITS \n");
    for(int i=0; i<20; i++){  
        printf("%f ",generation_TOP_fit[i]);
    }
    //Show Results
    //Save Results
    system("taskkill /F /IM kill_XFOIL.exe");
    printf("\n~~~~END~~~~");
}




/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */




/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ GA Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// These functions are related to the Genetic Algorithm. Besides the main genetic algorithm function, there are the functions 
// related to the Natural Selection, the Crossover and Mutation.

int Genetic_Algorithm(double *PAR_dbl, double *newPAR_dbl, unsigned short *PAR_ush, unsigned short *newPAR_ush, double **sol_PAR_dbl, double *fits, double * Fits2Prob, double alpha, double Re, double Ma, double *generation_TOP_fit){

    /*
    DESCRIPTION:
     It is the main function of the genetic algorithm. The inputs are basically the data structures used to save all the information, and the output is a 1 or a 0, meaning if the algorithm has worked or not.
     First, the initial population is created randomly. Then, the fitness is computed for each individual, and the best of them is saved into a variable called best_fitness. Then, while the solution is not 
     found or certain conditions are not met, the population is reproducing and evolving. First, a whole new generation is created via crossover and mutation. Then, the fitness is evaluated again and the best 
     of them is saved into new_best_fitness. If the new fitness is better than the previous one, best_fitness = new_best_fitness and a new reproduction is performed. However, if the new fitness is worst than 
     the one from the previous generation, a counter is activated and the reproduction continues. If this counter arrives at a certain value set by the programmer, the algorithm stops because it is considered 
     that the solution has been found and it can not be more optimised. In addition, if the whole population is infeasible, and it happens a certain number of consecutive times, it is considered that the problem 
     has ended without a solution and the algorithm also stops. It is also possible to stop the algorithm within a certain number of generations. 
    
    INPUTS:
     - double *PAR_dbl: pointer to a vector with the population parameters stored in double
     - double *newPAR_dbl: pointer to a vector with the population parameters after crossover stored in double
     - unsigned short *PAR_ush: same as *PAR_dbl but stored in unsigned short
     - unsigned short *newPAR_ush: same as *newPAR_dbl but stored in unsigned short
     - double **sol_PAR_dbl: pointer to a pointer that points where the best individual is
     - double *fits *Fits2Prob: pointers to a vector where the fits and the probabilities are stored
     - double alpha, double Re, double Ma: the flow conditions
     - double *generation_TOP_fit: pointer to a vector to store the top fit of each generation

     OUTPUTS:
     - int GA: an int with 0 of 1 depending on if the solution has been found or not
    */

    double *best_individual, *new_best_individual, *C1, *C2, *auxPAR_dbl, best_fit =__DBL_MAX__, new_best_fit =__DBL_MAX__, best_cl, best_cd, best_cm, Cl, Cd, Cm, *Cls = NULL, *Cds = NULL, *Cms = NULL; //cl_maxfit, cd_maxfit
    unsigned long Solution = 0, generation = 0, NewInd, counter_infeasible = 0, counter_No_Improve=0;
    int d;
    unsigned short *auxPAR_ush, *C1ush, *C2ush;

    Cls = (double *) malloc(individuals*sizeof(double));
    Cds = (double *) malloc(individuals*sizeof(double));
    Cms = (double *) malloc(individuals*sizeof(double));
    


    printf("Creating first Generation...\n");
    int ind = individuals;
    int pos;

    for(int i=0; i<individuals; i++){ 
        
        Generate_Random_parameters(PAR_dbl + i*12, PAR_ush + i*12);    
        fits[i] = callXfoil(alpha, Re, Ma, PAR_dbl+i*12, &Cl, &Cd, &Cm);
        Cls[i]  = Cl;
        Cds[i]  = Cd;
        Cms[i]  = Cm;

        if(fits[i]<best_fit){
            best_fit = fits[i];
            best_cl  = Cls[i];
            best_cd  = Cds[i];
            best_cm  = Cms[i];
            best_individual = PAR_dbl+i*12;
            pos = i;
        }
        PrintUserInfo(i,ind);
    }

    //Save results 
    FILE * Results;
        Results = fopen("Results.txt","a");
        if (Results == NULL) ExitError("the data file does not exist or cannot be opened",01);
        fprintf(Results,"%d %f %f %f %f %f \n",0, best_cl, best_cd, best_cl/best_cd ,best_fit, best_cm);
    if (fclose(Results)!= 0) ExitError("the data file does not exist or cannot be closed",02);


    printf("\nBEST FIT: %f \n",best_fit);
    printf("BEST Eff: %f \n",best_cl/best_cd);
    printf("BEST Cl: %f \n",best_cl);
    printf("BEST Cd: %f \n",best_cd);
    generation_TOP_fit[0] = best_fit;

    char temp[10];
    unsigned short buf=0;
    sprintf(temp, "airfoil_gen%d.dat",generation);
    GenerateAirfoil(best_individual,temp,&buf);
    printf("file %s created\n",temp);

    printf("-------------------------------------------------------\n");

    // CREATE NEW INDIVIDUALS
    while (Solution == 0){
        NewInd=0;            
        generation++;

        if(generation == MAX_GENERATIONS){
            printf("\nMaximum number of generations reached\n");
            Solution = 1;
        }

        printf("\n\nGeneration  %d   \n",generation);

        // Transform the vector of fits to vector of probailities
        Fits_To_Probability(fits, Fits2Prob);

        printf("   -> Reproducing\n");
        for(int i=0;i<individuals/2;i++){
            // select two individuals from old generation for mating (biassed in favour of the fitter ones)
            Natural_Selection (Fits2Prob, PAR_dbl, PAR_ush, &C1, &C2, &C1ush, &C2ush);

            // crossover of parents (Cs1 and Cs2) in bit 'd' (from right to left). The result is added in newCij
            Crossover(newPAR_dbl, newPAR_ush, &NewInd, C1, C2, C1ush, C2ush, d);
            PrintUserInfo(i,individuals/2);
        }
            
        // Swap pointers to continue working with Cij, not newCij

        //double
        auxPAR_dbl = PAR_dbl;
           PAR_dbl = newPAR_dbl;
        newPAR_dbl = auxPAR_dbl;
        
        //unsigned short
        auxPAR_ush = PAR_ush;
           PAR_ush = newPAR_ush;
        newPAR_ush = auxPAR_ush;

        //compute fitness of new generation
        new_best_fit = __DBL_MAX__;

        printf("  -> Computing fitness...\n");

        //system(command);
        for(int i=0; i<individuals; i++){

            //printf("%d ",i);
            fits[i]= callXfoil(alpha, Re, Ma, PAR_dbl+i*12, &Cl, &Cd, &Cm);
            Cls[i] = Cl;
            Cds[i] = Cd;
            Cms[i] = Cm;

            if(fits[i]<new_best_fit){
                new_best_fit = fits[i];
                best_cl = Cls[i];
                best_cd = Cds[i];
                best_cm = Cms[i];
                new_best_individual = PAR_dbl+i*12;
                pos = i;
            }
            PrintUserInfo(i,ind);
        }


        // Stopping algorithm
        if(new_best_fit == __DBL_MAX__){ // check if 10 generations do not have feasible individuals
            counter_infeasible++;
            printf("INFEASIBLE (%d/%d) \n",counter_infeasible,MAX_INFEASIBLE);
            if(counter_infeasible>MAX_INFEASIBLE){
                printf("Best fitness: infinite \n");
                break;
            }
        }else{ // check if population has converged (not fitter individuals created in 4 generations)
            counter_infeasible = 0;
            if(new_best_fit >= best_fit){

                counter_No_Improve++;
                printf("\n  ~~ NOT IMPROVNG! (%d/%d) ~~\n",counter_No_Improve,Generation_NOT_IMPROVE);
                char temp[10];
                sprintf(temp, "airfoil_gen%d.dat",generation);
                GenerateAirfoil(best_individual,temp,&buf);
                printf("file %s created\n",temp);
                printf("-------------------------------------------------------\n");

                if(counter_No_Improve == Generation_NOT_IMPROVE) {
                    printf("Not improve > %d \n",Generation_NOT_IMPROVE);
                    *sol_PAR_dbl = best_individual;
                     Solution    = 1;  //BOOL TRUE
                }
            }else{ // the most common case: new generation has at least one fitter individual
                printf("\n  ~~ IMPROVING! ~~\n");
                best_fit        = new_best_fit;
                best_individual = new_best_individual;
                counter_No_Improve=0;
                printf("\nBEST FIT: %f \n",best_fit);
                printf("BEST Eff: %f \n",best_cl/best_cd);
                printf("BEST Cl: %f \n",best_cl);
                printf("BEST Cd: %f \n",best_cd);

                char temp[10];
                sprintf(temp, "airfoil_gen%d.dat",generation);
                GenerateAirfoil(best_individual,temp,&buf);
                printf("file %s created\n",temp);
                printf("-------------------------------------------------------\n");
            }
        }


        //Save results 
        FILE * Results;
            Results = fopen("Results.txt","a");
            if (Results == NULL) ExitError("the data file does not exist or cannot be opened",01);
            fprintf(Results,"%d %f %f %f %f %f \n",generation, best_cl, best_cd, best_cl/best_cd ,best_fit, best_cm);
        if (fclose(Results)!= 0) ExitError("the data file does not exist or cannot be closed",02);

        generation_TOP_fit[generation] = best_fit;

    }
    if (Solution == 0){
        printf("Error ssss");
    }
    
    GenerateAirfoil(best_individual,"best_airfoil.dat",&buf);
    return Solution;
} //END GA

void Fits_To_Probability(double *fits, double *Fits2Prob){
    
    /*
    DESCRIPTION: 
     This function is the one that transforms the vector of fits into a vector of probabilities. 
     In order to do so, it distributes, in a vector from 0 to 1, a certain space for each individual proportional to its fitness. 
     The higher is the fitness, the lower is the space, having 0 space for infinite fitnesses. Thus, this function generates a vector 
     of doubles Fits2Prob which is a vector from 0 to 1 divided into as elements as individuals in the population, and the individuals 
     with a better fitness have more space assigned in this vector.

    INPUTS: 
    - double *fits: pointer to a vector of doubles to store the fits
    - double *Fits2Prob: pointer to a vector of doubles to store the probabilities

    OUPUTS: (void)
    */

    double tot=0;
    for(int i=0;i<individuals;i++){
        Fits2Prob[i]=1e10/fits[i];
        tot += Fits2Prob[i];
    }
    Fits2Prob[0]=Fits2Prob[0]/tot;
    for(int i=1;i<individuals;i++){
        Fits2Prob[i]=Fits2Prob[i-1]+Fits2Prob[i]/tot;
    }

}

void Natural_Selection (double *Fits2Prob, double *PAR_dbl, unsigned short *PAR_ush, double **C1, double **C2, unsigned short **C1ush, unsigned short **C2ush){

    /*
    DESCRIPTION:
     This function selects two individuals for reproduction. It is inside of a loop of half the individuals. In each iteration, it generates a random number 
     between 0 and 1 and looks into the Fits2Prob vector to see to which individual correspond this random number and selects it as the parent 1. 
     Then repeats this process for parent 2. As a consequence, individuals with a better fit to the environment will be more frequently selected, but also, 
     individuals with a low fit can have the opportunity to be selected and transfer its genetic code to the next generation. On the other hand, infeasible individuals 
     or individuals with a very high fitness value will hardly transfer their genetic code and will, eventually, disappear completely. 
     This function generates two pointers to the two selected parents information.

    INPUTS:
     - double *Fits2Prob: pointer to a vector of doubles where the probabilities are stored
     - double *PAR_dbl, unsigned short *PAR_ush: pointers to vectors where the airfoil parameters are stored (in double and unsigned short)
     - double **C1, **C2: pointer to a pointer that will point the position of the selected airfoils to reproduce (parent 1 and 2)
     - unsigned short **C1ush, **C2ush: the same as **C1 and **C2 but with unsigned short

    OUPUTS: (void)
    */

    double prob;
    prob = Random_UD(0,1);
    for(int i=0;i<individuals;i++){
        if(Fits2Prob[i]>prob){
            *C1 = PAR_dbl + i*12;  // C1 es el puntero a la ubi de Cij + tants cops la dimensio de la matriu individu, on "tants cops" es la i (l'individu en questio)
            *C1ush = PAR_ush + i*12;
            break;
        }
    }
    prob = Random_UD(0,1);
    for(int i=0;i<individuals;i++){
       if(Fits2Prob[i]>prob){
            *C2 = PAR_dbl + i*12;  // C1 es el puntero a la ubi de Cij + tants cops la dimensio de la matriu individu, on "tants cops" es la i (l'individu en questio)
            *C2ush = PAR_ush + i*12;
           break;
       }
    }
    
}


void Crossover(double *newPAR_dbl, unsigned short *newPAR_ush, unsigned long *NewInd, double *C1, double *C2, unsigned short *C1ush, unsigned short *C2ush, int d){
    /*
    DESCRIPTION:
    The main objective of this function is to generate a completely new population. In order to do so, it takes the two parents chromosomes, 
    and gene by gene, it calls the OnePointCrossoverMutation to generate the complete genome of the two children. The function saves the information 
    of the two new children in vectors newPAR_dbl and newPAR_ush. 

    INPUTS:
    - double *newPAR_dbl unsigned short *newPAR_ush: pointer to the vectors where the new parameters will be stored
    - unsigned long *NewInd: variable that is used as a reference to store the new individuals in the vectors newPAR
    - double *C1, *C2: not used in this fuction in this version of the code
    - unsigned short *C1ush, *C2ush: pointer to a pointer that will point the position of the selected airfoils to reproduce (parent 1 and 2)
    - int d:  not used in this fuction in this version of the code

    OUTPUTS: (void)
    */

    unsigned short f1,f2;

    for(int i=0;i<12;i++){

        OnePointCrossoverMutation(C1ush[i], C2ush[i], &f1, &f2, MUTATION_PROB);

        int check1 = 0, check2 = 0;

        // Calculem els doubles associats als ints
        double v1,v2;
        v1 = limits_inferiors[i]+f1*(limits_superiors[i]-limits_inferiors[i])/(pow(2,16)-1);
        v2 = limits_inferiors[i]+f2*(limits_superiors[i]-limits_inferiors[i])/(pow(2,16)-1);

        if(f1>65535 || f1 < 0){
            printf("OUT OF LIMITS\n");
        }
        if(f2>65535 || f2 < 0){
            printf("OUT OF LIMITS\n");
        }

        //copiem els enters que toquin
        newPAR_ush[(*NewInd)*12+i]     = f1;
        newPAR_ush[((*NewInd)+1)*12+i] = f2;

        //copiem els dobles que toquin
        newPAR_dbl[(*NewInd)*12+i]     = v1;
        newPAR_dbl[((*NewInd)+1)*12+i] = v2;

    }
    *NewInd = (*NewInd)+2; // Hem afegit 2 individus
}

#define USHRT_WIDTH 16
//#define ULONG_WIDTH __WORDSIZE

void OnePointCrossoverMutation(unsigned short p1,  unsigned short p2, unsigned short *f1, unsigned short *f2, float prob){ 

    /* 
    DESCRIPTION:
     This function takes as an input two genes and performs the crossover and mutation to generate the two genes of the offspring. 
     The explanation of the implementation of crossover and mutation is widely explained in the report: section  4.4 Crossover And Mutation IMPLEMENTATION
    
    INPUTS:
     - unsigned short p1 p2: value of the airfoil parameter to perform the crossover (parent information, a gene).
     - unsigned short *f1 *f2: pointer to a vector of unsigned shorts where the result of each crossover is stored (children information, the genes)
     - float prob: probability of the mutation to happen

    OUPUTS: (void)
    */

    unsigned short oneU = 1U, mask =  0xFFFFU << (1 + (USHRTran() % (USHRT_WIDTH - 1))), mask_mut1, mask_mut2;

    // OPTION 1 - 1 bit mutation //
    if(MUTATION == 1){
        mask_mut1 = (oneU << (USHRTran() % USHRT_WIDTH)); //^ canvia el bit de la posició en la que hi hagi un 1, la part de la dreta son tot 0 i un 1
        mask_mut2 = (oneU << (USHRTran() % USHRT_WIDTH)); 
    }
        
    // OPTION 2 - 2 bit mutation //
    if(MUTATION == 2){

        unsigned short n1 = 1, n2 = 1, pos1, pos2;
        while (n2 == n1){
            n1 = Random_UD(0,15);
            n2 = Random_UD(0,15);
        }
        if(n1>n2){
            pos1 = n1;
            pos2 = n2;
        }else{
            pos1 = n2;
            pos2 = n1;
        }
        mask_mut1 = ( ( oneU << (pos1-pos2) ) + 1 ) << pos2;
        n1 = 1; n2 = 1;
        while (n2 == n1){
            n1 = Random_UD(0,15);
            n2 = Random_UD(0,15);
        }
        if(n1>n2){
            pos1 = n1;
            pos2 = n2;
        }else{
            pos1 = n2;
            pos2 = n1;
        }
        mask_mut2 = ( ( oneU << (pos1-pos2) ) + 1 ) << pos2;

    }
    
    // OPTION 3 - 0101010101010101 //
    if(MUTATION == 3){
        mask_mut1 = 0x5555U;
        mask_mut2 = 0x5555U;
    }


    // CROSSOVER + MUTATION
    *f1 = (p1 & mask) | (p2 & ~mask);
    if(uniform() < prob){
        *f1 = (*f1) ^ mask_mut1;
        //printf("Mutation1\n");
    } 
    *f2 = (p2 & mask) | (p1 & ~mask);
    if(uniform() < prob){
        *f2 = (*f2) ^ mask_mut2;
        //printf("Mutation2\n");
    } 

}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ AIRFOIL PARAMETRISATION AND AERODYNAMIC COMPUTATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

//These functions are related with the computation of aerodynamic coefficients.

void Generate_Random_parameters(double *PAR_dbl, unsigned short *PAR_ush){
    /*
    DESCRIPTION:
     This function takes as an input two pointers to places in the memory where it has to allocate the data that will generate. 
     Its objective is to generate the 12 random PARSEC parameters for an individual. 
     It saves these parameters as double and also linearised as unsigned short for its use the genetic algorithm. 
     To do so, it generates a random number between two predefined limits for each PARSEC parameter. 
     The data is saved in the the vectors PAR_dbl and PAR_ush, which contain the information of each individual of the population,
     the first stored as double values and the second as unsigned short.

    INPUTS:
     - double *PAR_dbl: pointer to the vector of airfoil parameters stored as double
     - unsigned short *PAR_ush: pointer to the vector of airfoil parameters stored as unsigned short

    OUTPUTS: (void)

    */

    for( int i = 0; i<12;i++){
        PAR_dbl[i] = (double) Random_UD(limits_inferiors[i],limits_superiors[i]);
        PAR_ush[i] = (PAR_dbl[i]-limits_inferiors[i])/(limits_superiors[i]-limits_inferiors[i])*(pow(2,16)-1);
    }
}

void GenerateAirfoil(double *PAR_dbl, char a_name[], unsigned short *ALARM){
    /*
    DESCRIPTION
     This function is the one that solves the system of equations of the PARSEC parametrisation and generates 
     a file with the current airfoil coordinates. To solve the system, it uses the functions at the end of this code. 
     The generated file will be used by XFOIL to compute the aerodynamic coefficients.
    
    INPUTS:
     - double *PAR_dbl: pointer to the vector of airfoil parameters stored as double
     - char a_name[]: the name of the file with the airfoil data to be generated
     - unsigned short *ALARM: a variable to control if the data is generated correctly

    OUTPUTS: (void)
    */

    double rleup, rlelo, Xup, Zup, Xlo, Zlo, Zxxup, Zxxlo, Z_te, dZ_te, alpha_te, beta_te;

    rleup    = PAR_dbl[0]; //0.02;          // Radius leading edge Upper Surface
    rlelo    = PAR_dbl[1];//0.005;          // Radius leading edge Lower Surface
    Xup      = PAR_dbl[2]; //0.43;          // X Coordinate of Maximum thicnkess Upper Surface
    Zup      = PAR_dbl[3]; //0.12;          // Z Coordinate of Maximum thickess Upper Surface
    Xlo      = PAR_dbl[4]; //0.23;          // X Coordinate of Maximum thicnkess Lower Surface
    Zlo      = PAR_dbl[5]; //-0.018;        // Z Coordinate of Maximum thickess Lower Surface
    Zxxup    = PAR_dbl[6]; //-0.8;          // Curvature Uper Surface
    Zxxlo    = PAR_dbl[7]; //0.4;           // Curvature Lower Surface
    Z_te     = PAR_dbl[8]; //-0.001;        // Z Coordinate trailing edge
    dZ_te    = PAR_dbl[9]; //0;             // Thickness trailing edge. Set as 0
    alpha_te = PAR_dbl[10]; //-5*PI/180;    // alpha trailing edge
    beta_te  = PAR_dbl[11]; //25*PI/180;     // Beta trailing edge
    
    


    // A upper
    double Aup[6][6] = {  
        {1,                      1,                     1,                     1,                     1,                     1                    },
        {pow(Xup,1./2.),         pow(Xup,3./2.),        pow(Xup,5./2.),        pow(Xup,7./2.),        pow(Xup,9./2.),        pow(Xup,11./2.)      },
        {1./2.,                  3./2.,                 5./2.,                 7./2.,                 9./2.,                 11./2.               },
        {1./2.*pow(Xup,-1./2.),  3./2.*pow(Xup,1./2.),  5./2.*pow(Xup,3./2.),  7./2.*pow(Xup,5./2.),  9./2.*pow(Xup,7./2.),  11./2.*pow(Xup,9./2.)},
        {-1./4.*pow(Xup,-3./2.), 3./4.*pow(Xup,-1./2.), 15./4.*pow(Xup,1./2.), 35./4.*pow(Xup,3./2.), 63./4.*pow(Xup,5./2.), 99./4.*pow(Xup,7./2.)},
        {1,                      0,                     0,                     0,                     0,                     0                    } 
    };

    double RHSup[6][1]  = { 
        {Z_te+1./2.*dZ_te},
        {Zup},
        {tan(alpha_te-beta_te/2.)},
        {0}, 
        {Zxxup}, 
        {sqrt(rleup)} 
    };
    // A lower
    double Alo[6][6] = {  
        {1,                      1,                     1,                     1,                     1,                     1                    },
        {pow(Xlo,1./2.),         pow(Xlo,3./2.),        pow(Xlo,5./2.),        pow(Xlo,7./2.),        pow(Xlo,9./2.),        pow(Xlo,11./2.)      },
        {1./2.,                  3./2.,                 5./2.,                 7./2.,                 9./2.,                 11./2.               },
        {1./2.*pow(Xlo,-1./2.),  3./2.*pow(Xlo,1./2.),  5./2.*pow(Xlo,3./2.),  7./2.*pow(Xlo,5./2.),  9./2.*pow(Xlo,7./2.),  11./2.*pow(Xlo,9./2.)},
        {-1./4.*pow(Xlo,-3./2.), 3./4.*pow(Xlo,-1./2.), 15./4.*pow(Xlo,1./2.), 35./4.*pow(Xlo,3./2.), 63./4.*pow(Xlo,5./2.), 99./4.*pow(Xlo,7./2.)},
        {1,                      0,                     0,                     0,                     0,                     0                    } 
    };
    double RHSlo[6][1]  = { 
        {Z_te-1./2.*dZ_te},
        {Zlo},
        {tan(alpha_te+beta_te/2.)},
        {0}, 
        {Zxxlo}, 
        {-sqrt(rlelo)}
    };
    

    //SOLVING THE SYSTEM Aup*coefup=RHSup and Alo*coeflo=RHSlo
    double k = 6;
    double inv_up[6][6], inv_lo[6][6], cup[6][1], clo[6][1];  //coef up and coef lo

    if (determinant(Aup, k) == 0 || determinant(Alo, k)==0)
        ExitError("Inverse of the matrix is not possible. Determinant = 0.",03);
    else{
        cofactor(Aup, k, inv_up);  // Here the inverse matrix is generated (inv_up)
        cofactor(Alo, k, inv_lo);
    }
    multiplyMatrices(inv_up, RHSup, cup, 6,6,6,1); //inv_up*RHS_up = cup;
    multiplyMatrices(inv_lo, RHSlo, clo, 6,6,6,1);

    //We have cup i clo, tge coefficients. Now we have to generate the airfoil shape
    double x[N], xaux[N], zu[N], zl[N], flipzu[N], airfoil_x[2*N-1], airfoil_y[2*N-1];


    // Cosine distribution for x coordinates (from 0 to 1)
    for(int i=0;i<N;i++){
        x[i] = 0.5*(1-cos( PI*i/(N-1) ));
        xaux[i] = 1- x[i];
        zu[i] = 0; zl[i] = 0; //initializing
    }

    //Z upper and Z lower coordinates  (PARSEC FORMULATION)
    for(int i=1;i<=6;i++){
        for(int j=1;j<=N;j++){
            zu[j-1] += cup[i-1][0]*pow(x[j-1],i-0.5);
            zl[j-1] += clo[i-1][0]*pow(x[j-1],i-0.5);
        }
    }

    for(int j=1;j<N;j++){
        if (zl[j-1]>zu[j-1]){
                *ALARM = 1;
            } 
    }

    // FLIP Z upper 
    unsigned long count = 0;
    for(int i=N-1;i>=0;i--){
        flipzu[count] = zu[i];
        count++;
    }

    // Generate airfoil format
    for(int i=0;i<N;i++){
        airfoil_x[i] = xaux[i];
        airfoil_y[i] = flipzu[i];
    }
    count = 0;
    for(int i=N-1;i<(2*N-1);i++){
        airfoil_x[i] = x[count];
        airfoil_y[i] = zl[count];
        count++;
    }


    /* Generate file for the airfoil */
    FILE * airfoil;
        airfoil = fopen(a_name,"wt");
        if (airfoil == NULL) ExitError("the data file does not exist or cannot be opened",01);
        for(int i=0;i<2*N-1;i++){
            fprintf(airfoil,"%f %f\n",airfoil_x[i],airfoil_y[i]);
        }
    if (fclose(airfoil)!= 0) ExitError("the data file does not exist or cannot be closed",02);

}



double callXfoil(double alpha, double Re, double Ma, double *PAR_dbl, double *Cl, double *Cd, double *Cm){
    /*
    DESCRIPTION:
     This function calls XFOIL for execution. In order to do so, it generates an input file in which the specifications of the XFOIL analysis are set. 
     These specifications are the fluid conditions such as viscosity, Reynolds number, Mach number, and also the angle of attack to do the analysis. 
     The file created with GenerateAirfoil is used. Then the function calls XFOIL with a system command: system("xfoil.exe < xfoil.inp > xfoil.out"). 
     When XFOIL is called and the computation completed, it generates a series of output files. The function checks if all the files are generated 
     and if so, it reads the DATA.dat file, in wich the results of the analysis are saved. This file contains a list with the angle or angles of 
     attack analysed, and the Cl, Cd, Cm, and more data for these angles. Te function reads the data and saves it in memory. Finally, it deletes all 
     the generated files and computes the fitness function, which is the output of this function. Sometimes the XFOIL gets stuck during the execution 
     because the input airfoil has a strange shape and XFOIL calculus can not converge. If so, there is a backup code running in parallel that cheks 
     the elapsed time since the execution of XFOIL and, if it lasts more than 3 seconds to generate the output files, it means that XFOIL has cracked, 
     so XFOIL task is killed, and the fitness function is computed as infinite.

    INPUTS:
     - double alpha, Re, Ma: The flow conditions
     - double *PAR_dbl: Pointer to a vector with the airfoil parameters
     - double *Cl, *Cd, *Cm: Pointers to a vector to store the airfoil coefficents

    OUTPUTS:
     - doule fits: the value of the fitness function for an airfoil
    */

    char airfoil_name[] = "a_to_XFOIL.dat";
    unsigned short ALARM = 0;

    GenerateAirfoil(PAR_dbl, airfoil_name, &ALARM); // GENERATES AN AIRFOIL FILE WICHI WILL BE CALLED BY XFOIL
    
    double fits, Ef;

    //double alpha_i = 0, alpha_f = 2, delta_alpha = 0.5;

    double alpha_i = alpha, alpha_f = alpha, delta_alpha = 0.5;
    //double alpha_i = alpha, alpha_f = alpha+0.5, delta_alpha = 0.5;

    if(ALARM == 0){

        int n_alphas = (alpha_f-alpha_i)/delta_alpha+1;

        // Generate file to call Xfoil
        FILE * xfoil;
            xfoil = fopen("xfoil.inp","wt");
            if (xfoil == NULL) ExitError("the data file does not exist or cannot be opened",01);
            fprintf(xfoil,"plop\n");
            fprintf(xfoil,"G\n");
            fprintf(xfoil,"\n");
            fprintf(xfoil,"load %s\n",airfoil_name);
            fprintf(xfoil,"nombreairfoil\n");
            fprintf(xfoil,"OPER\n");
            fprintf(xfoil,"ITER 100\n");
            fprintf(xfoil,"\n\n");
            fprintf(xfoil,"OPER\n");
            fprintf(xfoil,"RE %f\n",Re);
            fprintf(xfoil,"MACH %f\n",Ma);
            fprintf(xfoil,"VISC\n");
            fprintf(xfoil,"PACC\n");
            fprintf(xfoil,"\n\n");
            fprintf(xfoil,"as %f %f %f\n",alpha_i,alpha_f,delta_alpha);
            fprintf(xfoil,"dump xfoil_dump.dat\n");
            fprintf(xfoil,"cpwr xfoil_cpwr.dat\n");
            fprintf(xfoil,"pwrt\n");
            fprintf(xfoil,"DATA.dat\n");
            fprintf(xfoil,"plis\n");
            fprintf(xfoil,"\n");
            fprintf(xfoil,"QUIT");

        if (fclose(xfoil)!= 0) ExitError("the data file does not exist or cannot be closed",02);

        clock_t t_start, t_end;
        int tmax = 3, check = 0, check2=0, header_lines = 12;
        double elapsed_time = 0;
        double alphas[n_alphas], Cls[n_alphas], Cds[n_alphas], Cms[n_alphas], aux1[n_alphas], aux2[n_alphas], aux3[n_alphas], aux4[n_alphas];

        // Inicialitzem a 0 per si no trobés res, tenir resultats com a -10
        for(int i = 0;i< n_alphas; i++){
            alphas[i] = 0;
            Cls[i]  = 0;
            Cds[i]  = 0;
            Cms[i]  = 0;
            //Efs[i] = 0;
        }


        t_start = clock();

        // CALL XFOIL
        system("xfoil.exe < xfoil.inp > xfoil.out");  // Sovint es clava aqui directament...

        // CHECK ELAPSED TIME AND KILL TASK IF NECESSARY
        int comptador = 0;
        while (elapsed_time < tmax){
            t_end = clock();
            elapsed_time = ((double) (t_end - t_start)) / CLOCKS_PER_SEC;
            
        
            // IF FILES EXIST, READ RESULTS FROM DATA.dat
            if(fexists("DATA.dat") && fexists("xfoil.inp") && fexists("xfoil.out") && fexists("xfoil_dump.dat") && fexists("xfoil_cpwr.dat")){

                // Compute number of lines of the file to see if there are results
                unsigned long lines = CountLines("DATA.dat",header_lines) -1 ;   //THERE IS AN EXTRA \n in DATA.dat, so we decrease the number of lines by 1
                
                if(lines > 0){ // If there's some result, save it into the vector of alphas, Cl and Cd
                    FILE *results;
                        if ((results = fopen ("DATA.dat", "r")) == NULL){
                            ExitError("the data file does not exist or cannot be opened", 01);
                        }
                        char buffer[100];
                        int i = 0, header_counter = 0;

                        // Save results skipping header lines
                        while(!feof(results)){
                            if(header_counter>=header_lines){
                                fscanf(results, "%lf %lf %lf %lf %lf %lf %lf", &alphas[i], &Cls[i], &Cds[i], &aux1[i], &Cms[i], &aux3[i], &aux4[i] );
                                i++;
                            }else{
                                fgets(buffer, 100, results);
                                header_counter ++;
                            }
                        }
                    fclose(results);
                    // EXIT THE WHILE ONCE THE RESULTS ARE SAVED
                    elapsed_time = tmax;
                    check = 1;

                }else{
                    elapsed_time = tmax;  // If lines = 0, skip the while and go, the results will be 0
                    check = 2;
                } 
            }
        }

        if( check == 1){
            remove("xfoil.inp");
            remove("xfoil.out");
            remove("xfoil_dump.dat");
            remove("xfoil_cpwr.dat");
            remove("DATA.dat");
            
            *Cl = Cls[0]; 
            *Cd = Cds[0];
            *Cm = Cms[0];

            int CONDITIONS = 1;
            double thickness = PAR_dbl[3]-PAR_dbl[5];

            if(*Cd<0.001)       CONDITIONS = 0;
            if(*Cl<0)           CONDITIONS = 0;
            if((*Cl)/(*Cd)>EFF_LIM) CONDITIONS = 0;
            if(*Cm < MIN_CM)    CONDITIONS = 0;
            if(thickness > MAX_THICKNESS || thickness < MIN_THICKNESS) CONDITIONS = 0;
            

            if(CONDITIONS){
                
                Ef = (*Cl)/(*Cd);

                if (FITNESS_TYPE == 1){
                    fits = 1/Ef;
                }
                if (FITNESS_TYPE == 2){
                    fits = (*Cl-ClT)*(*Cl-ClT) + *Cd;
                }
                if (FITNESS_TYPE == 3){
                    fits = 1 / *Cl;
                }
                if (FITNESS_TYPE == 4){
                    fits = *Cd;
                }

            }else{
                *Cl  = 0;
                *Cd  = 0;
                *Cm  = 0;
                Ef = 0;
                fits = __DBL_MAX__;
            } 

        }else if( check == 2){
            remove("xfoil.inp");
            remove("xfoil.out");
            remove("xfoil_dump.dat");
            remove("xfoil_cpwr.dat");
            remove("DATA.dat");
            *Cl  = 0;
            *Cd  = 0;
            *Cm  = 0;
            Ef  = 0;
            fits = __DBL_MAX__;


        }else{
            remove("xfoil.inp");
            remove("xfoil.out");
            printf("XFOIL task killed\n");
            *Cl  = 0;
            *Cd  = 0;
            Ef   = 0;
            fits = __DBL_MAX__;
        }

        if(fexists(airfoil_name)) remove(airfoil_name);

    }else{
        *Cl  = 0;
        *Cd  = 0;
        Ef  = 0;
        fits = __DBL_MAX__;
    }

    return fits;
                            
}


int fexists(const char * filename){
    /*
    DESCRIPTION:
     This function is used in callXfoil to verify if some files exist in the working directory. 
     It takes the file name as an input and the output is a 1 or a 0.
    
    INPUTS:
     - const char *filename: a pointer to the filename to check
    
    OUPUTS:
     - int: 1 or 0
    */

    FILE *file;
    if (file = fopen(filename, "r")){
        fclose(file);
        return 1;
    }
    return 0;
}

unsigned long CountLines(char *filename, int header_lines){
    /*
    DESCRIPTION:
     This function is used in callXfoil to count the lines of the output DATA.dat file to check if the computation has been successful. 
     Sometimes, even though XFOIL does not crack, the DATA.dat file is empty. In this case, it is also needed to compute the fitness as infinite.
    
    INPUTS:
     - char *filename: a pointer to the filename to analyse
     - int header_lines: the number of headlines to skip when counting the lines

    OUTPUTS:
     - unsigned long lines: the number of lines in the file minus the header lines
    */

    FILE *fpointer;
        if ((fpointer = fopen (filename, "r")) == NULL){
            ExitError("the data file does not exist or cannot be opened", 01);
        }
        char buffer[100];
        unsigned long lines = 0;
        int header_counter = 0;
        while(!feof(fpointer)){
            if(header_counter>=header_lines){
                fgets(buffer, 100, fpointer);
                lines++;
            }else{
                fgets(buffer, 100, fpointer);
                header_counter ++;
            }
        }   
    fclose(fpointer);
    return lines;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ RANDOM PARAMETERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// These functions are dedicated to generating reliable random numbers. The random processes in the genetic algorithm 
// are extremely important, and the use of \texttt{rand()} function from C is not an option since it is a pseudo-random
// generator in which the generated numbers are dependant of a seed \cite{WikipediaPseudorandomGenerator}. 
// The use of this function could generate biased answers, so this is why the random generator, ran1(), 
// has been chosen from the Numerical Recipes in C \cite{Press1988NumericalComputing}.

double Random_UD(double i, double j){
    /*
    DESCRIPTION:
    This function returns a random double between the two inputs, i and j. In order to do so, 
    it calls the function uniform(), which returns a random number between 0 and 1 with a uniform distribution
    and linearise the output to the input limits.
    INPUTS:
     - double i, j: The two limits to set the random number
    OUTPUTS:
     - double RD: a random nomber between the two limits
    */

    return i + uniform() * (j-i);
}

/*********************************
 * ran1 from "Numerical Recipes" *
**********************************/
#define IA 16807
#define IM 2147483647L
#define AM (1.0/IM)
#define IQ 127773L
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum){
    /*
    It is a function from the Numerical Recipes in C. This function generates a random number using the Minimal Standard 
    for its purpose between 0.0 ans 1.0 excluding the endpoint value. According to references, there are no known statistical 
    tests that ran1 has failed, except for some specific cases with specific sequences of calls \cite{Press1988NumericalComputing}. 
    The input is an initialiser obtained from randomize(). 

    INPUTS:
    - long *idum: an intitaliser defined in randomise ()

    OUTPUTS:
    - float: a random number
    */

    int j; long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1; else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=(int) (iy/(long)NDIV);
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
} /* ran1 */
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/* Utility functions to simplify the use of ran1 and idum */
long idum;

void randomize(void) {
    /*
    DESCRIPTION:
     It is used to initialise the random process at the beginning of the code. 
     It generates an initialiser taking data from the computer time.
    
    INPUTS: (none)
    OUTPUTS: (void)
    */

    idum = time(NULL); 
    if(idum > 0) idum = -idum;
}

float uniform(void) {
    /*
    DESCRIPTION:
     As explained, this function returns a random number in the interval (0,1) with a uniform distribution.
     In order to do so, it calls the function ran1 with input idum, a negative integer to 
     initialise the random generator. The initialisation is performed in the function randomize().
    
    INPUTS: (none)
    OUTPUTS: 
     - float RD: a random number between 0 and 1
    */

    return ran1(&idum);
} 

unsigned char random_bit(void){ 
    /*
    DESCRIPTION:
     This function returns randomly a bit (1 or 0) using the function ran1().
     Based on von-neuman observation ; rather inefficient.

    INPUTS: (none)
    OUTPUTS:
     - unsigned char f: a random bit (1 or 0)
    */

    unsigned char f, s;
	do { 
        f = 2*ran1(&idum); 
        s = 2*ran1(&idum); 
    } while (f == s);

	return f;
}

unsigned short USHRTran(void){ 
    /*
    DESCRIPTION:
     This function returns a random unsigned short, which means a number between 0 and 65535. 
     It is used in the crossover and mutation, when dealing with binary operators. 
     It uses the function random\_bit() for its purpose.
    
    INPUTS: (none)
    OUTPUTS:
     - unsigned short: a random unsigned short
    */

    register unsigned char i;
	unsigned short oneU = 1U, base = 0U;
	
	for(i=0; i < USHRT_WIDTH; i++) if(random_bit()) { base = oneU; break; }
    for(i++; i < USHRT_WIDTH; i++){ 
        base = base << 1;
		if(random_bit()) base = base | oneU;
	}
	return base;
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ USER INFORMATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
// These functions are dedicated to give information of the code execution to the user.

void ExitError(const char *miss, int errcode) {
    /* 
    DESCRIPTION:
     This function stops the execution if there is an unexpected error 
     and sends a message to the user indicating the type of error
    
    INPUTS: 
      - const char *miss: A pointer to the error message
      - int errorcode: The error code
    OUTPUTS: (void)
    */

    fprintf (stderr, "\nERROR: %s.\nStopping...\n\n", miss); exit(errcode);
} 

void PrintUserInfo(int i, int ind){
    /*
    DESCRIPTION:
     This function prints the evolution of a loop execution. Since in this algorithm several
     loops are used, it is useful to know the evolution of them and check if the execution is
     stuck at some point.

    INPUTS:
     - int i: The current index
     - int d: the number of individuals
    
    OUTPUT: (void)
    */

    double p=(double)i/(double)ind*100;
    if(p == 5){
        printf("  05%%...\n");
    }
    if(p == 10){
        printf("  10%%...\n");
    }
    if(p == 20){
        printf("  20%%...\n");
    }
    if(p == 30){
        printf("  30%%...\n");
    }
    if(p == 40){
        printf("  40%%...\n");
    }
    if(p == 50){
        printf("  50%%...\n");
    }
    if(p == 60){
        printf("  60%%...\n");
    }
    if(p == 70){
        printf("  70%%...\n");
    }
    if(p == 80){
        printf("  80%%...\n");
    }
    if(p == 90){
        printf("  90%%...\n");
    }
    if(i == ind-1){
        printf("  100%%\n");
    }
}


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SYSTEM SOLVER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// These functions are dedicated to solve the PARSEC system of equations to generate the airfoil coordinates from the 12 PARSEC parameters. As a remember, the system of the PARSEC parametrisation
// to solve is a system with the following form: A·x=RHS
// Where A is a 6x6 matrix and x and RHS a 6x1 vector. Thus, to solve such a system the following recombination can be done: x = A^{-1}·RHS
// All in all, it is needed to make the inverse of A function to solve the system, so the following functions have been created.
//
// These functions generate the inverse by dividing the adjugate matrix of A and the determinant of A, and then solve the system by multiplying the inverse of A by the Right Hand Side of the system.


double determinant(double a[6][6], double k){


    double s = 1, det = 0, b[6][6];
    int i, j, m, n, c;
    if (k == 1){
      return (a[0][0]);
    }
    else{
        det = 0;
        for (c = 0; c < k; c++){
            m = 0;
            n = 0;
            for (i = 0;i < k; i++){
                for (j = 0 ;j < k; j++){
                    b[i][j] = 0;
                    if (i != 0 && j != c){
                        b[m][n] = a[i][j];
                        if (n < (k - 2)) n++;
                        else{
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            det = det + s * (a[0][c] * determinant(b, k - 1));
            s = -1 * s;
        }
    }
 
    return (det);
}
 
void cofactor(double num[6][6], double f, double inv[6][6]){

    double b[6][6], fac[6][6];
    int p, q, m, n, i, j;
    for (q = 0;q < f; q++){
        for (p = 0;p < f; p++){
            m = 0;
            n = 0;
            for (i = 0;i < f; i++){
                for (j = 0;j < f; j++){
                    if (i != q && j != p){
                        b[m][n] = num[i][j];
                        if (n < (f - 2)) n++;
                        else{
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
        }
    } 
    transpose(num, fac, f, inv);
}

/*Finding transpose of matrix*/ 
void transpose(double num[6][6], double fac[6][6], double r, double inverse[6][6]){
    int i, j;
    double b[6][6], d;
 
    for (i = 0;i < r; i++){
        for (j = 0;j < r; j++){
            b[i][j] = fac[j][i];
        }
    }
    d = determinant(num, r);
    for (i = 0;i < r; i++){
        for (j = 0;j < r; j++){
            inverse[i][j] = b[i][j] / d;
        }
    }
    
}

void multiplyMatrices(double first[6][6], double second[6][1], double mult[6][1], int r1, int c1, int r2, int c2) {

    // Initializing elements of matrix mult to 0.
    for (int i = 0; i < r1; ++i) {
        for (int j = 0; j < c2; ++j) {
            mult[i][j] = 0;
        }
    }

    // Multiplying first and second matrices and storing in mult.
    for (int i = 0; i < r1; ++i) {
        for (int j = 0; j < c2; ++j) {
            for (int k = 0; k < c1; ++k) {
                mult[i][j] += first[i][k] * second[k][j];
            }
        }
    }

}
