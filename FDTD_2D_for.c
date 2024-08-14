//introduction:
/* 
   fdtd2d: < 2 order acoustic equation finite difference model > <finite difference order : 2-10>  
                  < doundary number: 30 ,but you can change the value of pml in this function>                                
*/
//----------------------Function Prototypes:---------------------------//
void fdtd2d(float **C, float *source, float **seis, float ***snap, int nx, int nz, int nt, float dx, float dz, float dt,
                  int isx, int isz, int rz, int flag, int pml, int nsnap, int dsnap);


//----------------------function defination:---------------------
//<--------------------------------------------------------------------------------->//
void fdtd2d(float **C, float *source, float **seis, float ***snap, int nx, int nz, int nt, float dx, float dz, float dt,
                 int isx, int isz, int rz, int flag, int pml, int nsnap, int dsnap)
{
    int   nx2=0, nz2=0, i=0, ix=0, iz=0, it=0, ixx=0, izz=0, mx=0, mz=0;
    float pi=3.1415926;
    float **P0=NULL, **P1=NULL, **P2=NULL, **vel=NULL, *wx=NULL, *wz=NULL, **w=NULL, **P31=NULL, **P32=NULL;
    float c[8][8], c0[8];
    float px=0.0, pz=0.0;
    FILE  *fid;
    char file[100];
 
    nx2=nx+2*pml;
    nz2=nz+2*pml;  
    wx =alloc1float(nx2);        wz =alloc1float(nz2);        w  =alloc2float(nz2,nx2); 
    P0 =alloc2float(nz2,nx2);    P1 =alloc2float(nz2,nx2);    P2 =alloc2float(nz2,nx2);
    P31=alloc2float(nz2,nx2);    P32=alloc2float(nz2,nx2);    vel=alloc2float(nz2,nx2);
    //>>>>>>>>>>>>>>>>>> velocity model >>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
    for(ix=0;ix<nx;ix++)
        for(iz=0;iz<nz;iz++)
        {
            vel[ix+pml][iz+pml]=C[ix][iz];
        }
    //the velocity of the left and right boundary;
    for(ix=0;ix<pml;ix++)
        for(iz=0;iz<nz2;iz++)
        {
            vel[ix][iz]=vel[pml][iz];
            vel[ix+nx+pml][iz]=vel[ix+nx+pml-1][iz];
        }
    //the velocity of the top and bottom boundary;
    for(iz=0;iz<pml;iz++)
        for(ix=0;ix<nx2;ix++)
        {
            vel[ix][iz]=vel[ix][pml];
            vel[ix][iz+nz+pml]=vel[ix][iz+pml+nz-1];
        }
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
    for(ix=0;ix<nx2;ix++)
        for(iz=0;iz<nz2;iz++)
        {
            P0 [ix][iz]=0.0;      P1 [ix][iz]=0.0;           P2 [ix][iz]=0.0;  
            P31[ix][iz]=0.0;      P32[ix][iz]=0.0;    
        }
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
    
    for(ix=0;ix<nx2;ix++)    wx[ix]=1.0;
    for(iz=0;iz<nz2;iz++)    wz[iz]=1.0;
    for(ix=1;ix<=pml;ix++)
    {  
        wx[ix-1]=cos(pi/2*(ix-pml-1)/pml);
        wx[nx2-ix]=wx[ix-1];
    }
    for(iz=1;iz<=pml;iz++)
    {
        wz[iz-1]=cos(pi/2*(iz-pml-1)/pml);
        //!wz(iz)=(iz-1)/Naz
        wz[nz2-iz]=wz[iz-1];
    }
    for(ix=1;ix<=nx2;ix++)
        for(iz=1;iz<=nz2;iz++)
        {
            w[ix-1][iz-1]=wx[ix-1]*wz[iz-1];
        }
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>// 
    c0[0]=-2.0;
    c0[1]=-5/2.0;
    c0[2]=-49/18.0;
    c0[3]=-205/72.0;
    c0[4]=-5269/1800.0;
    c0[5]=-5369/1800.0;
    c0[6]=-266681/88200.0;
    c0[7]=-1077749/352800.0; 

    c[0][0]=1.0;//2jie
    c[0][1]=4/3.0;c[1][1]=-1/12.0;//4jie
    c[0][2]=3/2.0;c[1][2]=-3/20.0;c[2][2]=1/90.0;//6jie
    c[0][3]=8/5.0;c[1][3]=-1/5.0 ;c[2][3]=8/315.0;c[3][3]=-1/560.0;//8jie
    c[0][4]=5/3.0;c[1][4]=-5/21.0;c[2][4]=5/126.0;c[3][4]=-5/1008.0;c[4][4]=1/3150.0;//10jie
    c[0][5]=12/7.0;c[1][5]=-15/56.0;c[2][5]=10/189.0;c[3][5]=-1/112.0;c[4][5]=2/1925.0;c[5][5]=-1/16632.0;//12jie
    c[0][6]=7/4.0;c[1][6]=-7/24.0;c[2][6]=7/108.0;c[3][6]=-7/528.0;c[4][6]=7/3300.0;c[5][6]=-7/30888.0;
                c[6][6]=1/84084.0;//14jie
    c[0][7]=16/9.0;c[1][7]=-14/45.0;c[2][7]=112/1485.0;c[3][7]=-7/396.0;c[4][7]=112/32175.0;c[5][7]=-2/3861.0;
                c[6][7]=16/315315.0;c[7][7]=-1/411840.0;//16jie
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
    mx=flag/2; mz=flag/2;
    for(it=1;it<=nt;it++)
    {
        //if(fmod(it,500)==0)            printf("it=%d\n",it);
        for(ix=2;ix<=nx2-1;ix++)
        for(iz=2;iz<=nz2-1;iz++)
        {    
            if(iz-mz<1)    mz=iz-1;       if(iz+mz>nz2)    mz=nz2-iz;
            if(ix-mx<1)    mx=ix-1;       if(ix+mx>nx2)    mx=nx2-ix;
            px=c0[mx-1]*P1[ix-1][iz-1];
            pz=c0[mz-1]*P1[ix-1][iz-1];
            for(ixx=1;ixx<=mx;ixx++)  px=px+c[ixx-1][mx-1]*(P1[ix+ixx-1][iz-1]+P1[ix-ixx-1][iz-1]);
            for(izz=1;izz<=mz;izz++)  pz=pz+c[izz-1][mz-1]*(P1[ix-1][iz+izz-1]+P1[ix-1][iz-izz-1]);
            P31[ix-1][iz-1]=vel[ix-1][iz-1]*dt/dx*vel[ix-1][iz-1]*dt/dx*px+vel[ix-1][iz-1]*dt/dz*vel[ix-1][iz-1]*dt/dz*pz
                            +2.0*P1[ix-1][iz-1]-P0[ix-1][iz-1];
            //if((ix==isx) && (iz==isz)) P31[ix-1][iz-1]=P31[ix-1][iz-1]+source[it-1]; //wave(it)*V(iz,ix)*V(iz,ix)*dt*dt
            //
            mx=flag/2;
            mz=flag/2;
            px=0.0;
            pz=0.0;
        }
        //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> boundary >>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
        for(iz=1;iz<=nz2;iz++)
        {
            for(ix=pml;ix>0;ix--)
            {
                P32[ix-1][iz-1]=(2.0-2.0*vel[ix-1][iz-1]*dt/dx-vel[ix-1][iz-1]*dt/dx*vel[ix-1][iz-1]*dt/dx)*P1[ix-1][iz-1]
                                +2.0*vel[ix-1][iz-1]*dt/dx*(1.0+vel[ix-1][iz-1]*dt/dx)*P1[ix][iz-1]-vel[ix-1][iz-1]*dt/dx*vel[ix-1][iz-1]*dt/dx*P1[ix+1][iz-1]
                                +(2.0*vel[ix-1][iz-1]*dt/dx-1.0)*P0[ix-1][iz-1]-2.0*vel[ix-1][iz-1]*dt/dx*P0[ix][iz-1];
            }
            for(ix=nx2-pml+1;ix<=nx2;ix++)
            {
                P32[ix-1][iz-1]=(2.0-2.0*vel[ix-1][iz-1]*dt/dx-vel[ix-1][iz-1]*dt/dx*vel[ix-1][iz-1]*dt/dx)*P1[ix-1][iz-1]
                                +2.0*vel[ix-1][iz-1]*dt/dx*(1.0+vel[ix-1][iz-1]*dt/dx)*P1[ix-2][iz-1]-vel[ix-1][iz-1]*dt/dx*vel[ix-1][iz-1]*dt/dx*P1[ix-3][iz-1]
                                +(2.0*vel[ix-1][iz-1]*dt/dx-1.0)*P0[ix-1][iz-1]-2.0*vel[ix-1][iz-1]*dt/dx*P0[ix-2][iz-1];
            }
        }
        for(ix=1;ix<=nx2;ix++)
        {
            for(iz=pml;iz>0;iz--)
            {
                P32[ix-1][iz-1]=(2.0-2.0*vel[ix-1][iz-1]*dt/dz-vel[ix-1][iz-1]*dt/dz*vel[ix-1][iz-1]*dt/dz)*P1[ix-1][iz-1]
                                +2.0*vel[ix-1][iz-1]*dt/dz*(1.0+vel[ix-1][iz-1]*dt/dz)*P1[ix-1][iz]-vel[ix-1][iz-1]*dt/dz*vel[ix-1][iz-1]*dt/dz*P1[ix-1][iz+1]
                                +(2.0*vel[ix-1][iz-1]*dt/dz-1.0)*P0[ix-1][iz-1]-2.0*vel[ix-1][iz-1]*dt/dz*P0[ix-1][iz];
            }
            for(iz=nz2-pml+1;iz<=nz2;iz++)
            {
                P32[ix-1][iz-1]=(2.0-2.0*vel[ix-1][iz-1]*dt/dz-vel[ix-1][iz-1]*dt/dz*vel[ix-1][iz-1]*dt/dz)*P1[ix-1][iz-1]
                                +2.0*vel[ix-1][iz-1]*dt/dz*(1.0+vel[ix-1][iz-1]*dt/dz)*P1[ix-1][iz-2]-vel[ix-1][iz-1]*dt/dz*vel[ix-1][iz-1]*dt/dz*P1[ix-1][iz-3]
                                +(2.0*vel[ix-1][iz-1]*dt/dz-1.0)*P0[ix-1][iz-1]-2.0*vel[ix-1][iz-1]*dt/dz*P0[ix-1][iz-2];
            }
        }
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
        for(ix=1;ix<=nx2;ix++)
            for(iz=1;iz<=nz2;iz++)
            {
                P2[ix-1][iz-1]=w[ix-1][iz-1]*P31[ix-1][iz-1]+(1.0-w[ix-1][iz-1])*P32[ix-1][iz-1];
            }
        //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
        // add source;
        P2[isx+pml][isz+pml]+=source[it-1]*vel[isx+pml][isz+pml]*vel[isx+pml][isz+pml]*dt*dt/dx/dx;
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<//
       
        for(ix=0;ix<nx;ix++)
        {
            seis[ix][it]=P2[ix+pml][rz+pml];
        }
        for(ix=0;ix<nx2;ix++)
            for(iz=0;iz<nz2;iz++)
            {
                P0[ix][iz] = P1[ix][iz];
                P1[ix][iz] = P2[ix][iz];            
            }
            
        if(it%dsnap == 0)
        {


            sprintf(file,"snapshot_%dms.bin",it);      
	    fid=fopen(file,"wb");   //snapshot data;
            for(ix=0;ix<nx;ix++)
            for(iz=0;iz<nz;iz++)
	    	fwrite(&P1[ix+pml][iz+pml], sizeof(float), 1, fid);

	    fclose(fid);
            //for(ix=0;ix<nx;ix++)
            //for(iz=0;iz<nz;iz++)
            //{
            //    snap[it/dsnap-1][ix][iz] = P2[ix+pml][iz+pml];
            //}
        }

            
    }

    
    //------------forward modeling end-----------------//
    free1float(wx);   free1float(wz);
    free2float(P0);   free2float(P1);    free2float(P2);    free2float(vel);   
    free2float(P31);  free2float(P32);   free2float(w);     
}
