#include <iostream>
#include <cmath>
#include <vector>

#define LX  0.2
#define LY  0.1
#define LZ  .02

#define DT .001;
#define NX  80
#define NY  40
#define NZ  16

#define U_TOPWALL 1
#define U_BOTTOMWALL -1
#define GRAV    -9.81
#define BETA    .003
#define RHO     1
#define CP  1000
#define MU  1.8e-5
#define KCOND .02
#define T_TOPWALL    400
#define T_BOTTOMWALL 300
#define DPDX    -1




using namespace std;

vector<vector<vector<vector<double>>>> nodeLocs(NX+1,vector<vector<vector<double>>>(NY+1,vector<vector<double>>(NZ+1,vector<double>(3))));
vector<vector<vector<vector<double>>>> centroidLocs(NX,vector<vector<vector<double>>>(NY,vector<vector<double>>(NZ,vector<double>(3))));
vector<vector<vector<double>>> P(NX,vector<vector<double>>(NY,vector<double>(NZ)));
vector<vector<vector<double>>> Pprime(NX,vector<vector<double>>(NY,vector<double>(NZ)));
vector<vector<vector<double>>> massbal(NX,vector<vector<double>>(NY,vector<double>(NZ)));
vector<vector<vector<double>>> T(NX,vector<vector<double>>(NY,vector<double>(NZ)));
vector<vector<vector<double>>> Tnew(NX,vector<vector<double>>(NY,vector<double>(NZ)));
vector<vector<vector<double>>> U(NX,vector<vector<double>>(NY,vector<double>(NZ)));
vector<vector<vector<double>>> V(NX,vector<vector<double>>(NY,vector<double>(NZ)));
vector<vector<vector<double>>> W(NX,vector<vector<double>>(NY,vector<double>(NZ)));
vector<vector<vector<double>>> Ustar(NX,vector<vector<double>>(NY,vector<double>(NZ)));
vector<vector<vector<double>>> Vstar(NX,vector<vector<double>>(NY,vector<double>(NZ)));
vector<vector<vector<double>>> Wstar(NX,vector<vector<double>>(NY,vector<double>(NZ)));
vector<vector<vector<double>>> Unew(NX,vector<vector<double>>(NY,vector<double>(NZ)));
vector<vector<vector<double>>> Vnew(NX,vector<vector<double>>(NY,vector<double>(NZ)));
vector<vector<vector<double>>> Wnew(NX,vector<vector<double>>(NY,vector<double>(NZ)));
vector<vector<vector<double>>> Hnx(NX,vector<vector<double>>(NY,vector<double>(NZ)));
vector<vector<vector<double>>> Hny(NX,vector<vector<double>>(NY,vector<double>(NZ)));
vector<vector<vector<double>>> Hnz(NX,vector<vector<double>>(NY,vector<double>(NZ)));
vector<vector<vector<double>>> Holdx(NX,vector<vector<double>>(NY,vector<double>(NZ)));
vector<vector<vector<double>>> Holdy(NX,vector<vector<double>>(NY,vector<double>(NZ)));
vector<vector<vector<double>>> Holdz(NX,vector<vector<double>>(NY,vector<double>(NZ)));
vector<double> uBarProfile(NZ+2);
vector<double> hxBarProfile(NZ);
static double dx = LX/NX;
static double dy = LY/NY;
static double dz = LZ/NZ;
static double Ax = dy*dz;
static double Ay = dx*dz;
static double Az = dx*dy;


double calcHx(int i, int j, int k);
double calcHy(int i, int j, int k);


int main()
{


cout << "Total Elements = " << NX*NY*NZ << " Total Volume = " << LX*LY*LZ << "  < DX, DY, DZ > = <"<<dx<<", "<<dy<<", "<<dz<<">\n";
cout << "Maximum Aspect Ratio = " << max(dx/dy,max(dx/dz,max(dy/dx,max(dz/dx,max(dy/dz,dz/dy))))) << endl;

double x=0;
double y=0;
double z=0;
int nNodes=0;
for (int i=0;i<NX+1;i++)// initialize Node Locations
{
    y=0;
    for (int j=0;j<NY+1;j++)
    {
        z=0;
        for (int k=0; k<NZ+1;k++)
        {
            nodeLocs[i][j][k][0]=x;
            nodeLocs[i][j][k][1]=y;
            nodeLocs[i][j][k][2]=z;
            //cout << nNodes<<" < X, Y, Z > = <"<<nodeLocs[i][j][k][0]<<", "<<nodeLocs[i][j][k][1]<<", "<<nodeLocs[i][j][k][2]<<">\n";
            z+=dz;
            nNodes++;
        }
        y+=dy;
    }
    x+=dx;
}
   
int nCells=0;//initialize Scalar Locations
x=0.5*dx;
for (int i=0;i<NX;i++)
{
    y=0.5*dy;
    for(int j=0;j<NY;j++)
    {
        z=0.5*dz;
        for(int k=0;k<NZ;k++)
        {
            centroidLocs[i][j][k][0]=x;
            centroidLocs[i][j][k][1]=y;
            centroidLocs[i][j][k][2]=z;
            P[i][j][k]=0+DPDX*centroidLocs[i][j][k][0];
            Pprime[i][j][k]=0;
            massbal[i][j][k]=0;
            T[i][j][k]=T_BOTTOMWALL+(T_TOPWALL-T_BOTTOMWALL)*centroidLocs[i][j][k][2]/LZ;;
            U[i][j][k]=0;
            V[i][j][k]=0;
            W[i][j][k]=0;
            Ustar[i][j][k]=U_BOTTOMWALL+(U_TOPWALL-U_BOTTOMWALL)*centroidLocs[i][j][k][2]/LZ; 
            Vstar[i][j][k]=0;
            Wstar[i][j][k]=0;
            Unew[i][j][k]=0;
            Vnew[i][j][k]=0;
            Wnew[i][j][k]=0;
            nCells++;
            //cout << nCells<<" < X, Y, Z > = <"<<centroidLocs[i][j][k][0]<<", "<<centroidLocs[i][j][k][1]<<", "<<centroidLocs[i][j][k][2]<<">\n";
            z+=dz;
        }
    y+=dy;
    }
x+=dx;
}

cout << "bBarProfile.size() = " << uBarProfile.size() <<endl;
return 0;
}

double calcHx(int i, int j, int k)
{
    double hx, up, ue, uw, un, us, uu, ud, vn, vs, wu, wd, dudxe, dudxw, dudyn, dudys, dudzu, dudzd;
    up = U[i][j][k];
    
    //Ueast
    if(i == NX-1) ue = 0.5*(up+U[0][j][k]);
    else ue = 0.5*(up+U[i+1][j][k]);       
        
    //Uwest
    if(i == 0) uw = 0.5*(up+U[NX-1][j][k]);
    else uw = 0.5*(up+U[i-1][j][k]);
    
    //Unorth
    if(j==NY-1) un = 0.5*(up+U[i][0][k]);
    else un = 0.5*(up+U[i][j+1][k]);
        
    //Usouth
    if(j == 0) us = 0.5*(up+U[i][NY-1][k]);
    else us = 0.5*(up+U[i][j-1][k]);
    
    //Uup
    if(k == NZ-1) uu = U_TOPWALL;
    else uu = 0.5*(up+U[i][j][k+1]);
        
    //Udown
    if(k == 0) ud = U_BOTTOMWALL;
    else ud = 0.5*(up+U[i][j][k-1]);
      
    //Vnorth
    if(j == NY-1) vn = 0.5*(V[i][j][k]+V[i+1][j][k]);
    else vn = 0.5*(V[i][j][k]+V[i+1][j][k]);
    
    //Vsouth
    if(j == 0) vs = 0.5*(V[i][NY-1][k]+V[i+1][NY-1][k]);
    else vs = 0.5*(V[i][j-1][k]+V[i+1][j-1][k]);

    //Wup
    if(k == NZ-1) wu = 0;
    else wu = 0.5*(W[i][j][k]+W[i+1][j][k]);

    //Wdown
    if (k == 0) wd = 0;
    else wd = 0.5*(W[i][j][k-1]+W[i+1][j][k-1]);
    
    //Viscous Stresses
    dudxe = 2*(ue-up)/dx;
    dudxw = 2*(up-uw)/dx;
    dudyn = 2*(un-up)/dy;
    dudys = 2*(up-us)/dy;
    dudwu = 2*(uu-up)/dz;
    dudwd = 2*(up-ud)/dz;
    
    hx = ()/dx+()/dy+()/dz-MU/(RHO*dx
    
    
    
}

            
            


    





