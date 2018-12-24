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

static vector<vector<vector<vector<double>>>> nodeLocs(NX+1,vector<vector<vector<double>>>(NY+1,vector<vector<double>>(NZ+1,vector<double>(3))));
static vector<vector<vector<vector<double>>>> centroidLocs(NX,vector<vector<vector<double>>>(NY,vector<vector<double>>(NZ,vector<double>(3))));
static vector<vector<vector<double>>> P(NX,vector<vector<double>>(NY,vector<double>(NZ)));
static vector<vector<vector<double>>> Pprime(NX,vector<vector<double>>(NY,vector<double>(NZ)));
static vector<vector<vector<double>>> massbal(NX,vector<vector<double>>(NY,vector<double>(NZ)));
static vector<vector<vector<double>>> T(NX,vector<vector<double>>(NY,vector<double>(NZ)));
static vector<vector<vector<double>>> Tnew(NX,vector<vector<double>>(NY,vector<double>(NZ)));
static vector<vector<vector<double>>> U(NX,vector<vector<double>>(NY,vector<double>(NZ)));
static vector<vector<vector<double>>> V(NX,vector<vector<double>>(NY,vector<double>(NZ)));
static vector<vector<vector<double>>> W(NX,vector<vector<double>>(NY,vector<double>(NZ)));
static vector<vector<vector<double>>> Ustar(NX,vector<vector<double>>(NY,vector<double>(NZ)));
static vector<vector<vector<double>>> Vstar(NX,vector<vector<double>>(NY,vector<double>(NZ)));
static vector<vector<vector<double>>> Wstar(NX,vector<vector<double>>(NY,vector<double>(NZ)));
static vector<vector<vector<double>>> Unew(NX,vector<vector<double>>(NY,vector<double>(NZ)));
static vector<vector<vector<double>>> Vnew(NX,vector<vector<double>>(NY,vector<double>(NZ)));
static vector<vector<vector<double>>> Wnew(NX,vector<vector<double>>(NY,vector<double>(NZ)));
static vector<vector<vector<double>>> Hnx(NX,vector<vector<double>>(NY,vector<double>(NZ)));
static vector<vector<vector<double>>> Hny(NX,vector<vector<double>>(NY,vector<double>(NZ)));
static vector<vector<vector<double>>> Hnz(NX,vector<vector<double>>(NY,vector<double>(NZ)));
static vector<vector<vector<double>>> Holdx(NX,vector<vector<double>>(NY,vector<double>(NZ)));
static vector<vector<vector<double>>> Holdy(NX,vector<vector<double>>(NY,vector<double>(NZ)));
static vector<vector<vector<double>>> Holdz(NX,vector<vector<double>>(NY,vector<double>(NZ)));
static vector<double> uBarProfile(NZ);
static vector<double> uRMSProfile(NZ);
static vector<double> hxBarProfile(NZ);
static vector<double> hyBarProfile(NZ);
static vector<double> hzBarProfile(NZ);
static vector<double> zLocsProfile(NZ);
static double dx = LX/NX;
static double dy = LY/NY;
static double dz = LZ/NZ;
//static double Ax = dy*dz;
//static double Ay = dx*dz;
//static double Az = dx*dy;


void calcHx();
double calcHy();
double calcHz();
void calcAvgProfiles();
void initializeNodeLocations();
void initializeScalarLocations();


int main()
{

initializeNodeLocations();
initializeScalarLocations();


calcHx();
calcAvgProfiles();


return 0;
}

void calcHx()
{
    double up, ue, uw, un, us, uu, ud, vn, vs, wu, wd, dudxe, dudxw, dudyn, dudys, dudzu, dudzd;
    int iplusone, jplusone, iminusone, jminusone;

    for (int i=0;i<NX;i++) //Update Hxyz
    {
        for(int j=0;j<NY;j++)
        {
            for(int k=0;k<NZ;k++)
            {
                //cout <<"check1 "<<" "<<i<<" "<<j<<" "<<k<<endl;
                Holdx[i][j][k]=Hnx[i][j][k];// Push Current Hnx into past
                up = U[i][j][k];

                if (i==NX-1) iplusone = 0;
                    else iplusone = i+1;
                if (j==NY-1) jplusone = 0;
                    else jplusone = j+1;
                if (i==0) iminusone = 0;
                    else iminusone = i-1;
                if (j==0) jminusone = 0;
                    else jminusone = j-1;

                //Ueast
                ue = 0.5*(up+U[iplusone][j][k]);

                //Uwest
                uw = 0.5*(up+U[iminusone][j][k]);

                //Unorth
                un = 0.5*(up+U[i][jplusone][k]);

                //Usouth
                us = 0.5*(up+U[i][jminusone][k]);

                //Uup
                if(k == NZ-1) uu = U_TOPWALL;
                else uu = 0.5*(up+U[i][j][k+1]);

                //Udown
                if(k == 0) ud = U_BOTTOMWALL;
                else ud = 0.5*(up+U[i][j][k-1]);

                //Vnorth
                vn = 0.5*(V[i][j][k]+V[iplusone][j][k]);

                //Vsouth
                vs = 0.5*(V[i][jminusone][k]+V[iplusone][jminusone][k]);

                //Wup
                wu = 0.5*(W[i][j][k]+W[iplusone][j][k]);

                //Wdown
                wd = 0.5*(W[i][j][k-1]+W[iplusone][j][k-1]);

                //Viscous Stresses
                dudxe = 2*(ue-up)/dx;
                dudxw = 2*(up-uw)/dx;
                dudyn = 2*(un-up)/dy;
                dudys = 2*(up-us)/dy;
                dudzu = 2*(uu-up)/dz;
                dudzd = 2*(up-ud)/dz;

                //Set Hnx
                Hnx[i][j][k] = (ue*ue-uw*uw)/dx+(un*vn-us*vs)/dy+(uu*wu-ud*wd)/dz-MU*(dudxe-dudxw)/(RHO*dx)-MU*(dudyn-dudys)/(RHO*dy)-MU*(dudzu-dudzd)/(RHO*dz);
            }
        }
    }
    return;
}

void calcAvgProfiles()
{
    for (int k=0;k<NZ;k++)// initialize Node Locations
    {
        uBarProfile[k]=0;
        //uRMSProfile[k]=0;
        hxBarProfile[k]=0;
        hyBarProfile[k]=0;
        hzBarProfile[k]=0;
        zLocsProfile[k]=centroidLocs[0][0][k][2];

        for (int j=0;j<NY;j++)
        {
            for (int i=0; i<NX;i++)
            {
                uBarProfile[k]+=U[i][j][k]/(NX*NY);
                hxBarProfile[k]+=Hnx[i][j][k]/(NX*NY);
            }
        }
        cout<<zLocsProfile[k]<<" "<<uBarProfile[k]<<" "<<hxBarProfile[k]<<endl;
    }
    return;

}

void initializeNodeLocations()
{
    double x=0, y=0,z=0;
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

    return;
}

void initializeScalarLocations()
{
    int nCells=0;
    double x=0.5*dx,y,z;
    for (int i=0;i<NX;i++) //initialize Scalar Locations
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
                P[i][j][k]=0+DPDX*centroidLocs[i][j][k][0]+RHO*GRAV*centroidLocs[i][j][k][2];
                Pprime[i][j][k]=0;
                massbal[i][j][k]=0;
                T[i][j][k]=T_BOTTOMWALL+(T_TOPWALL-T_BOTTOMWALL)*centroidLocs[i][j][k][2]/LZ;
                U[i][j][k]=0;//U_BOTTOMWALL+(U_TOPWALL-U_BOTTOMWALL)*centroidLocs[i][j][k][2]/LZ;
                V[i][j][k]=0;
                W[i][j][k]=0;
                Ustar[i][j][k]=0;
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
    cout << "Total Elements = " << NX*NY*NZ << " Total Volume = " << LX*LY*LZ << "  < DX, DY, DZ > = <"<<dx<<", "<<dy<<", "<<dz<<">\n";
    cout << "Maximum Aspect Ratio = " << max(dx/dy,max(dx/dz,max(dy/dx,max(dz/dx,max(dy/dz,dz/dy))))) << endl;

    return;
}









