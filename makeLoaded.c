/**************
 * makeLoaded *
 *************/
/* 
 * Make initial data file for particle-loaded bead simulation
 * Some parts adapted from chain.f in the LAMMPS distribution
 * input: number of monomers (beads) in chain and bond length.
 * output xyz file of chain 
 *
 */

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "stdlib.h"
#define PI 3.1415926

#define PARTCENTERTYPE 1
#define PARTTYPE 2
#define VINYLBEADTYPE 3
#define SILICONEBEADTYPE 4

#define atomTypes 4
#define bondTypes 2
#define angleTypes 1

int main(int argc,char* argv[])
{
    int i,j,k;
    int regen;
    int numStrands,strandBeads,numPart,partBeads;
    double partRad,partVol,polyVol,cutoff,rhostar_i,rhostar_f,rbond,partSphereRad,ri0i2;
    double packingFrac,halfEdge_f,onePartVol;
    int numBeads,bead,mol,bType,aType;
    double x[5];
    int y[2];
    char ele[10];
    double x0,x1,x2,y0,y1,y2,z0,z1,z2,volume_i,volume_f,boxEdge_f;
    float ran0(),pbc();
    double xprd,yprd,zprd,r,dr,rsq;
    double xinner,yinner,zinner;
    double dx,dy,dz,bpx,bpy,bpz;
    double xsurf,ysurf,zsurf;
    double xboundlo,xboundhi,
           yboundlo,yboundhi,
           zboundlo,zboundhi; 
    double **point;
    double edge,origin;
    FILE *fp,*op;
    char inFile[80],outFile[156],buf[360],shell[180];
    bool noFiller = false;

    if(argc!=3){
        fprintf(stderr,"Usage: %s <input file name> <dummy>.\n",argv[0]);
        exit(1);
    }

    // read in input data file 
    sprintf(inFile,"%s",argv[1]);
    fprintf(stderr,"opening input file %s...\n",inFile);
    fp = fopen(inFile,"r");
    fgets(buf,360,fp); // header
    fgets(buf,360,fp); // blank line
    fgets(buf,360,fp); // 
    sscanf(buf,"%d",&y[0]);
    strandBeads=y[0];
    fgets(buf,360,fp); // 
    sscanf(buf,"%lf",&x[0]);
    ri0i2=x[0];
    fgets(buf,360,fp); // 
    sscanf(buf,"%lf",&x[0]);
    rbond=x[0];
    fgets(buf,360,fp); // 
    sscanf(buf,"%lf",&x[0]);
    rhostar_i=x[0];
    fgets(buf,360,fp); // 
    sscanf(buf,"%lf",&x[0]);
    rhostar_f=x[0];
    fgets(buf,360,fp); // blank line
    fgets(buf,360,fp); // 
    sscanf(buf,"%lf",&x[0]);
    packingFrac=x[0];
/*    fgets(buf,360,fp); // 
    sscanf(buf,"%d",&y[0]);
    // numPart is dummy since using volume basis
    numPart=y[0]; */
    fgets(buf,360,fp); // 
    sscanf(buf,"%lf",&x[0]);
    partRad=x[0];
/*
    fgets(buf,360,fp); // 
    sscanf(buf,"%d",&y[0]);
    partBeads=y[0]; */
    // re-assign since using volume basis 
    partBeads=(int)(4*PI*pow(partRad,2.0));
    fgets(buf,360,fp); // 
    sscanf(buf,"%lf",&x[0]);
    partSphereRad=x[0];
    fgets(buf,360,fp); // 
    sscanf(buf,"%lf",&x[0]);
    cutoff=x[0];
    fgets(buf,360,fp); // 
    sscanf(buf,"%lf",&x[0]);
    boxEdge_f=x[0];

    fclose(fp);

    if(packingFrac==0.0){
        noFiller = true;
    }
    
    if(noFiller){
        numPart = 0;
        partVol = 0.0;
        polyVol = pow(boxEdge_f,3.0); // new for this set of runs// jjk //strandBeads<50?pow(25.0*2.0,3.0)*rhostar_f:pow((25.0+((strandBeads-50)/3.0))*2.0,3.0)*rhostar_f;
    }
    else{
    //partVol = numPart*(4./3.)*PI*(pow(partSphereRad+partRad,3.));
    onePartVol = (4./3.)*PI*(pow(partSphereRad+partRad,3.));
    numPart = (int)((packingFrac*pow(boxEdge_f,3.0))/onePartVol);
    partVol = numPart*onePartVol;
    polyVol = pow(boxEdge_f,3.)-partVol; //partVol/packingFrac*(1-packingFrac);
    }
    numStrands = (int)(polyVol*rhostar_f/strandBeads);

    volume_i = partVol+polyVol*rhostar_f/rhostar_i;
    xprd = pow(volume_i,(1.0/3.0));
    yprd = xprd;
    zprd = xprd;

    xboundlo = -xprd/2.0;
    xboundhi = -xboundlo;
    yboundlo = xboundlo;
    yboundhi = xboundhi;
    zboundlo = xboundlo;
    zboundhi = xboundhi;

    if(noFiller){
        volume_f = pow(boxEdge_f,3.); //polyVol;
    }
    else{
        volume_f = pow(boxEdge_f,3.); //partVol/packingFrac;
    }
    halfEdge_f = pow(volume_f,(1./3.))/2.;

    numBeads = numStrands*strandBeads + numPart*partBeads;

    fprintf(stdout,"rhostar_i = %f\n",rhostar_i);
    fprintf(stdout,"rhostar_f = %f\n",rhostar_f);
    fprintf(stdout,"pack frac = %f\n",packingFrac);
    fprintf(stdout,"cutoff = %f\n",cutoff);
    fprintf(stdout,"num filler particles = %d\n",numPart);
    fprintf(stdout,"num particles on NP = %d\n",partBeads);
    fprintf(stdout,"half edge, initial = %f\n",xboundhi);
    fprintf(stdout,"half edge, final = %f\n",halfEdge_f);
    fprintf(stdout,"num strands = %d, num beads = %d\n",numStrands, numBeads);
    
    sprintf(outFile,"README.txt");
    op=fopen(outFile,"w");
    fprintf(op,"polymer density, initial, rhostar_i = %f\n",rhostar_i);
    fprintf(op,"polymer density, final, rhostar_f = %f\n",rhostar_f);
    fprintf(op,"NP packing fraction = %f\n",packingFrac);
    fprintf(op,"cutoff distance (used in setup only) = %f\n",cutoff);
    fprintf(op,"num filler nanoparticles(NP) = %d\n",numPart);
    fprintf(op,"num particles on NP = %d\n",partBeads);
    fprintf(op,"half edge, initial = %f\n",xboundhi);
    fprintf(op,"half edge, final = %f\n",halfEdge_f);
    fprintf(op,"number of polymer strands = %d\n",numStrands);
    fprintf(op,"beads per polymer strand = %d\n",strandBeads);
    fprintf(op,"total beads = %d\n",numBeads);
    fclose(op);
    
    point = malloc(numBeads*sizeof(double *));
    for(i=0;i<numBeads;i++){
        point[i] = malloc(8*sizeof(double));
    }
    // update LAMMPS input file
    sprintf(shell,"sed -e 's/HALFSIM/%f/g' in.template > in.run",halfEdge_f);
    system(shell);

    if(!noFiller){
    // make one particle
    sprintf(shell,"./sphere %5.2f %d",partRad,partBeads-1);
    system(shell);

    // PACKMOL some spheres into a relevant box
    edge=pow(volume_i,1./3.)-1.;
    origin=-(edge/2.)+1.;
    sprintf(shell,"sed -e 's/ORIGIN/%f/g' makePartBox.inp | sed -e 's/S_EDGE/%f/' | sed -e 's/NUMPART/%d/' > tmpINP",origin,edge,numPart);
//    fprintf(stderr,"%s",shell);
    system(shell);
    // PACKMOL $PATH HERE
    system("~/packmol/packmol < tmpINP");

    // read in particle xyz's 
    // & assign mol & type
    sprintf(inFile,"np_box.xyz");
    fprintf(stderr,"opening file %s...\n",inFile);
    fp = fopen(inFile,"r");
    fgets(buf,360,fp);
    sscanf(buf,"%d",&y[0]);
    if(numPart*partBeads != y[0]){
        fprintf(stderr,"Mismatch, number of sphere particles.Expected %d, found %d.\n",numPart*partBeads,y[0]);
        exit(0);
    }
    fgets(buf,360,fp); // comment line
    fprintf(stderr,"%s\n",buf);
    for(i=0;i<numPart;i++){
        for(j=0;j<partBeads;j++){
            bead = i*partBeads+j;
            fgets(buf,360,fp);
            //            fprintf(stderr,"%s",buf);
            sscanf(buf,"%s %lf %lf %lf",ele,&x[0],&x[1],&x[2]);
            point[bead][0]=i+1; // mol
            if(j==0){ // bead center
                point[bead][1]=PARTCENTERTYPE;
            }
            else{
                point[bead][1]=PARTTYPE;
            }
            point[bead][2]=x[0];
            point[bead][3]=x[1];
            point[bead][4]=x[2];
            point[bead][5] = point[bead][6] = point[bead][7] = 0.0; //iflags
        }
    }
    fclose(fp);
    fprintf(stderr,"closed %s\n",inFile);
    }
    // now generate and store the chains
MAKE_POLY: i=0;
    for(i=0;i<numStrands;i++){
        mol = numPart+1+i; 

RAND_START:j=0;
           bead = numPart*partBeads+i*strandBeads+j;

           x1 = 0.0;
           y1 = 0.0;
           z1 = 0.0;
           x2 = xboundlo + ran0()*xprd; // uniform random, 0 to 1
           y2 = yboundlo + ran0()*yprd;
           z2 = zboundlo + ran0()*zprd;
           //check distances to particle centers
           if(!noFiller){
           for(k=0;k<numPart;k++){ //don't think we need PBC check here
               bpx = x2-point[k*partBeads][2];
               bpy = y2-point[k*partBeads][3];
               bpz = z2-point[k*partBeads][4];
               if(sqrt(bpx*bpx+bpy*bpy+bpz*bpz) < partRad+cutoff){
                   goto RAND_START;
               }
           }
           }
           //           fprintf(stderr,"start is ok\n");

           //generate random point
           point[bead][0]=mol;
           point[bead][1]=VINYLBEADTYPE; //vinyl end
           point[bead][2]=x2;
           point[bead][3]=y2;
           point[bead][4]=z2;
           point[bead][5]=0;
           point[bead][6]=0;
           point[bead][7]=0;

           for(j=1;j<strandBeads;j++){
               bead++;
               x0 = x1;
               y0 = y1;
               z0 = z1;
               x1 = x2;
               y1 = y2;
               z1 = z2;
               regen=0;
RAND_PART: rsq=99.0;
           regen++;
           xinner = 2.0*ran0() - 1.0;
           yinner = 2.0*ran0() - 1.0;
           zinner = 2.0*ran0() - 1.0;
           rsq = xinner*xinner + yinner*yinner + zinner*zinner;
           if(rsq>1.0) goto RAND_PART;

           //project onto surface of unit radius sphere
           r = pow(rsq,0.5);
           xsurf = xinner/r;
           ysurf = yinner/r;
           zsurf = zinner/r;
           // scale to bond length
           x2 = x1 + xsurf*rbond;
           y2 = y1 + ysurf*rbond;
           z2 = z1 + zsurf*rbond;

           // check distance constraint
           dx = x2 - x0;
           dy = y2 - y0;
           dz = z2 - z0;
           r = sqrt(dx*dx + dy*dy + dz*dz);

           if (j>2 && r<=ri0i2) goto RAND_PART;

           if(!noFiller){
           //check distances to particle centers
           for(k=0;k<numPart;k++){ //don't think we need PBC check here
               bpx = x2-point[k*partBeads][2];
               bpy = y2-point[k*partBeads][3];
               bpz = z2-point[k*partBeads][4];
               if(sqrt(bpx*bpx+bpy*bpy+bpz*bpz) < partRad+cutoff){
                   if(regen > 9999){
                       fprintf(stderr,"not finding it... lets's start over and try again...\n");
                       goto MAKE_POLY;    
                   }
                   goto RAND_PART;
               }
           }
           }
           // move new bead within PBC
           if (x2 < xboundlo) x2 += xprd;
           if (x2 >= xboundhi) x2 -= xprd;
           if (y2 < yboundlo) y2 += yprd;
           if (y2 >= yboundhi) y2 -= yprd;
           if (z2 < zboundlo) z2 += zprd;
           if (z2 >= zboundhi) z2 -= zprd;

           // store new point
           point[bead][0] = mol;
           
           if(j==strandBeads-1){
               point[bead][1] = VINYLBEADTYPE;
           }
           else{
               point[bead][1] = SILICONEBEADTYPE;
           }
           point[bead][2] = x2;
           point[bead][3] = y2;
           point[bead][4] = z2;

           // check and store image flags
           for(k=0;k<3;k++){
               dr = point[bead][2+k]-point[bead-1][2+k];
               if(fabs(dr) < 2.0*rbond){
                   point[bead][5+k]=point[bead-1][5+k];
               }
               else if(dr < 0.0){
                   point[bead][5+k]=point[bead-1][5+k]+1;
               }
               else if(dr > 0.0){
                   point[bead][5+k]=point[bead-1][5+k]-1;
               }
               else{
                   fprintf(stderr,"logic mixup! exiting...\n");
                   exit(0);
               }
           }
           }
    }

    // now print the LAMMPS data file
    sprintf(outFile,"loaded.data");
    op=fopen(outFile,"w");
    fprintf(op,"Bead-spring with particles\n\n");
    fprintf(op,"%d atoms\n",numBeads);
    fprintf(op,"%d bonds\n",numStrands*(strandBeads-1));
    fprintf(op,"%d angles\n",numStrands*(strandBeads-2));
    fprintf(op,"0 dihedrals\n");
    fprintf(op,"0 impropers\n\n");
    fprintf(op,"%d atom types\n",atomTypes);
    fprintf(op,"%d bond types\n",bondTypes);
    fprintf(op,"%d angle types\n",angleTypes);
    fprintf(op,"0 dihedral types\n");
    fprintf(op,"0 improper types\n\n");
    fprintf(op,"%f %f xlo xhi\n",xboundlo,xboundhi);
    fprintf(op,"%f %f ylo yhi\n",yboundlo,yboundhi);
    fprintf(op,"%f %f zlo zhi\n\n",zboundlo,zboundhi);

    fprintf(op,"Masses\n\n"); // massess all =1, for now
    for(i=0;i<atomTypes;i++){
        fprintf(op,"%d %f\n",i+1,1.0);
    }
    fprintf(op,"\nAtoms\n\n");
    for(i=0;i<numBeads;i++){
        fprintf(op,"%d %d %d 0.000 %f %f %f %d %d %d\n",i+1,(int)point[i][0],(int)point[i][1],point[i][2],point[i][3],point[i][4],(int)point[i][5],(int)point[i][6],(int)point[i][7]);
    }
    fprintf(op,"\nBonds\n\n");
    for(i=0;i<numStrands;i++){
        for(j=0;j<strandBeads-1;j++){
            bead=numPart*partBeads+i*strandBeads+j;
            if(j==0||j==strandBeads-2){
                bType=1;
            }
            else{
                bType=2;
            }
            fprintf(op,"%d %d %d %d\n",i*(strandBeads-1)+j,bType,bead+1,bead+2);
        }
    }
    fprintf(op,"\nAngles\n\n");
    for(i=0;i<numStrands;i++){
        for(j=0;j<strandBeads-2;j++){
            bead=numPart*partBeads+i*strandBeads+j;
            if(j==0||j==strandBeads-3){
                aType=1;
            }
            else{
                aType=1;
            }
            fprintf(op,"%d %d %d %d %d\n",i*(strandBeads-2)+j,aType,bead+1,bead+2,bead+3);
        }
    }
    fclose(op);
}

float ran0(){
    return (double)rand() / (double)RAND_MAX;
}
//}
