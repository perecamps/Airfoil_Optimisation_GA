// PROGRAM TO KILL THE XFOIL PROCESS IF IT GETS STUCK WHILE RUNNING THE GENETIC ALGORITHM. 
// TIME: 3s

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define tmax 3

int fexists(const char * filename){
    /* try to open file to read */
    FILE *file;
    if (file = fopen(filename, "r")){
        fclose(file);
        return 1;
    }
    return 0;
}

clock_t t_start, t_end;
double elapsed_time;
int condicio;

int main(){
    int a = 0;
    while (a == 0){
        condicio = 0;
        t_start = clock();
        while (condicio == 0){
            t_end = clock();
            elapsed_time = ((double) (t_end - t_start)) / CLOCKS_PER_SEC;
            if(fexists("DATA.dat") && fexists("xfoil.inp") && fexists("xfoil.out") && fexists("xfoil_dump.dat") && fexists("xfoil_cpwr.dat")){
                //printf("1\n");
                condicio = 1;
            }else{
                if (elapsed_time>tmax){
                    printf("killing\n");
                    system("taskkill /F /IM xfoil.exe");
                    condicio = 1;
                }

            }

        }
    }

}
