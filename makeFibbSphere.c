/******************
 * makeFibbSphere *
 *****************/
/* input: a radius (r) and number of particles (N)
 * output xyz file containing centers of particles
 * on a sphere of radius r using Fibonacci sphere
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "stdlib.h"
#define PI 3.1415926

int main(int argc,char* argv[])
{
    int i,j,k;
    int numPart;
    double x,y,z;
    double radius,th,phi;
    double rx,ry,rz,rrad,v1,v2,fac,rsq,norm; 

    double **point;
    FILE *op;
    char outfile[156];

    if(argc!=3){
        fprintf(stderr,"Usage: %s <radius of big particle> <number of particles on surface>.\n",argv[0]);
        exit(1);
    }
    radius=atof(argv[1]);
    numPart=atoi(argv[2]);

    point = malloc(numPart*sizeof(double *));
    for(i=0;i<numPart;i++){
        point[i] = malloc(3*sizeof(double));
    }
    srand(time(0));
    phi=PI*(3.0-pow(5.0,0.5));
    for(i=0;i<numPart;i++){
        ry=1.0-((float)(i)/(numPart-1.0))*2.0;
        rrad=pow(1.0-ry*ry,0.5);
        th=phi*i;
        rx=cos(th)*rrad;
        rz=sin(th)*rrad;
        x=point[i][0]=radius*rx;
        y=point[i][1]=radius*ry;
        z=point[i][2]=radius*rz;

//        fprintf(stderr,"C %.5f %.5f %.5f\n",x,y,z);
    }

    // output as an xyz file
    sprintf(outfile,"np.xyz");
    op = fopen(outfile,"w");
    fprintf(op,"%d\nparticles on surface of sphere\n",numPart+1);
    fprintf(op,"C %f %f %f\n",0.0000,0.0000,0.0000);
    for(i=0;i<numPart;i++){
        fprintf(op,"C %f %f %f\n",point[i][0],point[i][1],point[i][2]);
    }
    fclose(op);
}

