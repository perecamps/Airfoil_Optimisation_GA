# AirfoilOpt_GeneticAlg
Airfoil Optimisation with a Genetic Algorithm. Final Master Thesis of Master's Degree in Aerospace Engineering for the UPC.
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
In this Master Thesis, an airfoil optimisation using a Genetic Algorithm is developed.  This project has been proposed by myself and done with the guidance and supervision of professor Manel Soria.\\

The main idea of the project is to develop from scratch an algorithm capable of finding the optimal airfoil for specific flow conditions, such as the angle of attack, the Reynolds number, and the Mach number. The objective is to create a useful tool for aerospace engineering students so they can use it on their projects and designs during the college years.\\

The work has a first theoretical part about Genetic Algorithms, in which the basic concepts needed to understand the current project are explained. Then, the implementation of the algorithm is fully expounded and all the intern processes of the genetic algorithm can be consulted. Several validations of the code have also been made.\\

The Genetic Algorithm created uses crossovers and mutations. The airfoil parametrisation used has been the PARSEC parametrisation and the computation of the aerodynamic coefficients is done with XFOIL. The whole code is written in C language and the results analysis and graphs are done with MATLAB and XFLR5.\\

Finally, the algorithm is tested with two real design cases, an airfoil for a heavy lifter aircraft that participated in the Air Cargo Challenge 2017 in Stuttgart, and an airfoil for a glider that flew in the Paper Air Challenge 2015 in ESEIAAT, Terrasa. The results and improvements offered by the algorithm are compared with the results that the designers of these aircraft obtained manually during the design process.\\

USER GUIDE
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

In order to execute the code, for now, it is needed a computer with a gcc compiler, which is included in the Minimalist GNU for Windows (MinGW). In addition, the software XFOIL.exe must be in the same folder as the executable. To modify the parameters, the .c file must be opened. In the first lines of the code, the modifiable parameters can be easily identified. These are the following:

· Genetic Algorithm Parameters

    - Generations_NOT_IMPROVE - Set as 8
    - MAX_INFEASIBLE - Set as 9.
    - MAX_GENERATIONS - Set as 100.
    - MUTATION - Set as type 1
    - MUTATION_PROB - Set as 0.7 (70%).
    - individuals - Set as 160 (more than 200 is not recommended)

· Fitness function

    - FITNESS_TYPE - To choose between 1, 2, 3 and 4. Recomended 1 and 2.
    - ClT - Objective Cl. Useful only if fitness type 2 is selected.
    
· Constraints

    - MAX_THICKNESS - To choose between 0 and 1. Set 1 for no constraint.
    - MIN_THICKNESS - To choose between 0 and 1. Set 0 for no constraint.
    - TE_THICKNESS - Set 0 for no constraint
    - MIN_TE_ANGLE - To choose between 1 and 25. Set 1 for no constraint
    - MIN_CM - To choose between -1 and 0. Set -1 for no constraint
    - EFF_LIM - Sometimes an error occurs in XFOIL execution and it returns huge values of Eff (because it computes very low values of Cd). To avoid this error, an Efficiency         limit related to the problem dimension is set. If the computed efficiency overcomes this value, the fit is set as infinite and the airfoil is not considered. It is               recommended to set it as 500 or so, and then modify it after a first execution, taking into account the efficiency values computed.

The execution generates several files, including the best airfoil of each generation (.dat files) and a .txt file containing the aerodynamic coefficients and fit of the best airfoil of each generation. These files can be analysed with a MATLAB code. Thus, it is needed to have the software installed to use it.

All in all, the needed files to execute the code and analyse the results are:

    - AirfoilOot_GA_v2.c
    - AirfoilOot_GA_v2.exe
    - kill_XFOIL.c
    - kill_XFOIL.exe
    - XFOIL.exe
    - plot_airfoils.m
    
The MATLAB file generates a plot in which all the airfoils and its evolution can be seen in the Appendices. It also shows the C_l and C_d evolution, the C_m evolution, the Efficiency evolution and the fitness evolution.
