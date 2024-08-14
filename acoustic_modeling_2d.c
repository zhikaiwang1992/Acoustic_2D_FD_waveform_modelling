// 本程序用于实现均匀介质中二维声波方程规则网格高阶差分数值模拟.
// 方程形式为：d2u/d2t=v*v*[d2u/d2x+d2u/d2z];
// 本程序最高可以实现空间16阶，时间2阶;


#include "su.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "FDTD_2D_for.c"


/*********************** self documentation **********************/
char *sdoc[] = {NULL};

//main function;
int main (int argc, char **argv)
{
    int   i=0, j=0, k=0, L=0, m=0, M=0, frw=0;
    int   nx, nz, nt, flag, sx, sz, rz, abc, dsnap, nsnap;
    float dx, dz, dt, f;

    char * srcpath = "";
    char * velpath = "";
    char * outpath = "";
    char file[100];
	
	float pi=3.1415926;
	float c[8][8]={0.0};
	float *S, **v;
	float **seismic, ***snapshot;
	float tdelay;
	FILE  *fidout, *fidin;

	printf("****** start ******\n");
	 
    initargs(argc, argv);

 	//requestdoc(1);
 	
    if (!getparint("nx",&nx))      nx  = 301;
    if (!getparint("nz",&nz))      nz  = 301;
    if (!getparint("nt",&nt))      nt  = 301;
    if (!getparint("sx",&sx))      sx  = nx/2;
    if (!getparint("sz",&sz))      sz  = nz/2;
    if (!getparint("rz",&rz))      rz  = 0;
    if (!getparint("abc",&abc))    abc = 50;
    if (!getparint("flag",&flag))  flag = 6;
    if (!getparint("dsnap",&dsnap)) dsnap = 5;

    if (!getparfloat("dx",&dx))    dx = 10.0;
    if (!getparfloat("dz",&dz))    dz = 10.0;
    if (!getparfloat("dt",&dt))    dt = 0.001;
    if (!getparfloat("f",&f))      f  = 5.0;

    if (!getparstring("out",&outpath));   //err("You must input the outpath");
    if (!getparstring("vel",&velpath));     //err("must input the velocity!\n");
    if (!getparstring("src",&srcpath));

    nsnap = (int)(nt/dsnap);

    printf("nx=%d, nz=%d  nt=%d\n",nx,nz,nt);
    printf("dx=%f, dz=%f, dt=%f\n",dx,dz,dt);
    printf("sx=%d, sz=%d, f =%f\n",sx,sz,f);
    printf("abc=%d, fd_order=%d\n",abc,flag);
        
	S = alloc1float(nt);
	v = alloc2float(nz, nx);
	seismic  = alloc2float(nt, nx);   // array for seismic  data;
	snapshot = alloc3float(nz, nx, nsnap);   // array for snapshot data;
	
    // velocity model
    fidin=fopen(velpath,"rb");
    frw=fread(v[0], sizeof(float), nz*nx, fidin);
    fclose(fidin);
	
  
	tdelay = 1.0/f;
	for(k=0;k<nt;k++)   
	{
	   S[k]=(1.0-2.0*(pi*f*(k*dt-tdelay))*(pi*f*(k*dt-tdelay)))*exp(-1.0*(pi*f*(k*dt-tdelay))*(pi*f*(k*dt-tdelay)));
	}

	fdtd2d(v, S, seismic, snapshot, nx, nz, nt, dx, dz, dt, sx, sz, rz, flag, abc, nsnap, dsnap);
	
	printf("finish modelling\n");
	
	// write seismic and snapshot to seismic.dat and snapshot.dat respectivelity;  
	/*
	for(k=0;k<nsnap;k++)
	{
	    sprintf(file,"%s/snapshot_%dms.bin",outpath,k*dsnap);      
	    fidout=fopen(file,"wb");   //snapshot data;
	    for(j=0;j<nx;j++)
	    for(i=0;i<nz;i++)
	    {
	    	fwrite(&snapshot[k][j][i], sizeof(float), 1, fidout);
	    }
	    fclose(fidout);
        }
        */
	sprintf(file,"%s/shotgather.bin",outpath);      
	fidout=fopen(file, "wb");   //seismic data;
	for(i=0; i<nx; i++)
	    fwrite(seismic[i], sizeof(float), nt, fidout);
	    
	fclose(fidout);
	free1float(S);
	free2float(v);
    free2float(seismic); 
    free3float(snapshot);

	printf("****** finish ******\n"); 
	return 0;
}

