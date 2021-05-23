# AirfoilOpt_GeneticAlg
Airfoil Optimisation with a Genetic Algorithm. Final Master Thesis of Master's Degree in Aerospace Engineering for the UPC.
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
In this Master Thesis, an airfoil optimisation using a Genetic Algorithm is developed.  This project has been proposed by myself and done with the guidance and supervision of professor Manel Soria.\\

The main idea of the project is to develop from scratch an algorithm capable of finding the optimal airfoil for specific flow conditions, such as the angle of attack, the Reynolds number, and the Mach number. The objective is to create a useful tool for aerospace engineering students so they can use it on their projects and designs during the college years.\\

The work has a first theoretical part about Genetic Algorithms, in which the basic concepts needed to understand the current project are explained. Then, the implementation of the algorithm is fully expounded and all the intern processes of the genetic algorithm can be consulted. Several validations of the code have also been made.\\

The Genetic Algorithm created uses crossovers and mutations. The airfoil parametrisation used has been the PARSEC parametrisation and the computation of the aerodynamic coefficients is done with XFOIL. The whole code is written in C language and the results analysis and graphs are done with MATLAB and XFLR5.\\

Finally, the algorithm is tested with two real design cases, an airfoil for a heavy lifter aircraft that participated in the Air Cargo Challenge 2017 in Stuttgart, and an airfoil for a glider that flew in the Paper Air Challenge 2015 in ESEIAAT, Terrasa. The results and improvements offered by the algorithm are compared with the results that the designers of these aircraft obtained manually during the design process.\\
