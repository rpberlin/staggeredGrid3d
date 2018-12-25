#include <iostream>
#include <cmath>
#include <vector>

#define LX  0.05
#define LY  0.01
#define LZ  0.01

#define DT 0.00025
#define NX  600
#define NY  200
#define NZ  200

#define NMAXSTEPS 1000000

#define U_TOPWALL 1
#define U_BOTTOMWALL -1
#define BETA    0
#define TREF    300
#define RHO     1
#define CP  1000
#define MU  1.8e-7
#define KCOND .02
#define T_TOPWALL    400
#define T_BOTTOMWALL 300
#define DPDX    0

using namespace std;

static vector<double> gravity{0,0,0};
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


//InitializeFunctions
void calcHx();
void calcHy();
void calcHz();
void calcUstar();
void calcVstar();
void calcWstar();
double calcPprime(double urf);
double updatePressureField();
double runPressureCorrector();
double calcUnew();
double calcVnew();
double calcWnew();
double checkMassBal();
void calcAvgProfiles();
void initializeNodeLocations();
void initializeScalarLocations();


int main()
{

initializeNodeLocations();
initializeScalarLocations();

cout<<"Mass Bal: "<< checkMassBal()<<endl;
double pResid, pResidNorm=-1e6, xResid, xResidNorm=1e-16, yResid, yResidNorm=1e-16, zResid, zResidNorm=1e-16;

for (int nOuter =1;nOuter<NMAXSTEPS;nOuter++)
{

    //Update Convective Stresses
    calcHx();
    calcHy();
    calcHz();
    calcUstar();
    calcVstar();
    calcWstar();



    pResid=runPressureCorrector();
    updatePressureField();
    if(nOuter<5) pResidNorm = max(pResid,pResidNorm);
    xResid = calcUnew();
    if(nOuter<5) xResidNorm = max(xResid,xResidNorm);
    yResid = calcVnew();
    if(nOuter<5) yResidNorm = max(yResid,yResidNorm);
    zResid = calcWnew();
    if(nOuter<5) zResidNorm = max(xResid,zResidNorm);

    cout<<nOuter<<" "<<pResid/pResidNorm<<" "<< xResid/xResidNorm<<" "<<yResid/yResidNorm<<" "<<zResid/zResidNorm<<" "<<endl;
//    for(int i=0;i<NX;i++)
//    {
//        for(int k=NZ-1;k>=0;k--)
//        {
//            //cout<<massbal[i][0][k]<<" ";
//            //cout<<Pprime[i][0][k]<<" ";
//        }
//        cout<<endl;
//    }

}

calcAvgProfiles();

return 0;
}

double runPressureCorrector()
{
    double urf=0.8;
    int nMaxSweeps=10000;
    double l2norm,normFactor=1e-16;


    for(int nInner=0;nInner<nMaxSweeps;nInner++)
    {
        l2norm = calcPprime(urf);
        if(nInner <5) normFactor = max(normFactor,l2norm);
//        cout<<ninner<<"  "<<l2norm<<endl;
        if (l2norm/normFactor < 1e-8) break; //converged
    }
    return normFactor;
}

double calcPprime(double urf)
{

    double l2NormPprime=0, deltaPprime,invdxsq, invdysq, invdzsq, invDenom;
    double PprimeNew, Pp, Pe, Pw, Pn, Ps, Pu, Pd, ue, uw, vn, vs, wu, wd;
    invdxsq=1.0/(dx*dx);
    invdysq=1.0/(dy*dy);
    invdzsq=1.0/(dz*dz);
    invDenom=1.0/((2.0*invdxsq)+(2.0*invdysq)+(2.0*invdzsq));

    for (int i=0;i<NX;i++) //loop through all scalar cells
    {
        for(int j=0;j<NY;j++)
        {
            for(int k=0;k<NZ;k++)
            {
                Pp=Pprime[i][j][k];

                //West Face
                if(i==0)
                {
                    uw=Ustar[NX-1][j][k];
                    Pw=Pprime[NX-1][j][k];
                }
                else
                {
                    uw=Ustar[i-1][j][k];
                    Pw=Pprime[i-1][j][k];
                }

                //East Face
                ue=Ustar[i][j][k];
                if(i==NX-1) Pe=Pprime[0][j][k];
                else        Pe=Pprime[i+1][j][k];

                //South Face
                if(j==0)
                {
                    vs=Vstar[i][NY-1][k];
                    Ps=Pprime[i][NY-1][k];
                }
                else
                {
                    vs=Vstar[i][j-1][k];
                    Ps=Vstar[i][j-1][k];
                }

                //North face
                vn=Vstar[i][j][k];
                if(j==NY-1)Pn=Pprime[i][0][k];
                else       Pn=Pprime[i][j+1][k];

                //Down Face
                if(k==0) //Wall => w=0 and dP/dz=0
                {
                    Pd=Pp; //Wall
                    wd=0; //Wall
                }
                else
                {
                    Pd=Pprime[i][j][k-1];
                    wd=Wstar[i][j][k-1];
                }

                //Up Face
                if(k==NZ-1)
                {
                    Pu=Pp;
                    wu=0;
                }
                else
                {
                    Pu=Pprime[i][j][k+1];
                    wu=Wstar[i][j][k];
                }
                PprimeNew=invDenom*(invdxsq*(Pe+Pw) + invdysq*(Pn+Ps) + invdzsq*(Pu+Pd) - (RHO/DT)*(((ue-uw)/dx)+((vn-vs)/dx)+((wu-wd)/dx)));
                deltaPprime=PprimeNew-Pp;
                l2NormPprime+=deltaPprime * deltaPprime;
                Pprime[i][j][k]=Pp+deltaPprime*urf;

            }
        }
    }

return sqrt(l2NormPprime);

}

double updatePressureField()
{
    double Pref = P[0][0][0]+Pprime[0][0][0];
    for (int i=0;i<NX;i++)
    {
        for(int j=0;j<NY;j++)
        {
            for(int k=0;k<NZ;k++)
            {
                P[i][j][k]+=Pprime[i][j][k]-Pref;
            }
        }
    }

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
                // Push Current Hnx into past
                Holdx[i][j][k]=Hnx[i][j][k];
                up = U[i][j][k];

                //do periodicity in i,j directions (x,y)
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
                if(k == NZ-1) uu = U_TOPWALL; //check for wall
                else uu = 0.5*(up+U[i][j][k+1]);

                //Udown
                if(k == 0) ud = U_BOTTOMWALL; //check for wall
                else ud = 0.5*(up+U[i][j][k-1]);

                //Vnorth
                vn = 0.5*(V[i][j][k]+V[iplusone][j][k]);

                //Vsouth
                vs = 0.5*(V[i][jminusone][k]+V[iplusone][jminusone][k]);

                //Wup
                if(k==NZ-1) wu=0; //check for wall
                else wu = 0.5*(W[i][j][k]+W[iplusone][j][k]);

                //Wdown
                if(k==0) wd =0; //check for wall
                else wd = 0.5*(W[i][j][k-1]+W[iplusone][j][k-1]);

                //Viscous Stresses
                dudxe = 2*(ue-up)/dx;
                dudxw = 2*(up-uw)/dx;
                dudyn = 2*(un-up)/dy;
                dudys = 2*(up-us)/dy;
                dudzu = 2*(uu-up)/dz;
                dudzd = 2*(up-ud)/dz;

                //Set Hnx
                Hnx[i][j][k] = (MU*(dudxe-dudxw)/(RHO*dx)+MU*(dudyn-dudys)/(RHO*dy)+MU*(dudzu-dudzd)/(RHO*dz))-((ue*ue-uw*uw)/dx+(un*vn-us*vs)/dy+(uu*wu-ud*wd)/dz);
            }
        }
    }
    return;
}

void calcHy()
{
    double vp, ve, vw, vn, vs, vu, vd, ue, uw, wu, wd, dvdxe, dvdxw, dvdyn, dvdys, dvdzu, dvdzd;
    int iplusone, jplusone, iminusone, jminusone;

    for (int i=0;i<NX;i++)
    {
        for(int j=0;j<NY;j++)
        {
            for(int k=0;k<NZ;k++)
            {

                // Push Current Hnx into past
                Holdy[i][j][k]=Hny[i][j][k];
                vp = V[i][j][k];

                //do periodicity in i,j directions (x,y)
                if (i==NX-1) iplusone = 0;
                    else iplusone = i+1;
                if (j==NY-1) jplusone = 0;
                    else jplusone = j+1;
                if (i==0) iminusone = 0;
                    else iminusone = i-1;
                if (j==0) jminusone = 0;
                    else jminusone = j-1;

                //Veast
                ve = 0.5*(vp+V[iplusone][j][k]);

                //Vwest
                vw = 0.5*(vp+V[iminusone][j][k]);

                //Vnorth
                vn = 0.5*(vp+V[i][jplusone][k]);

                //Vsouth
                vs = 0.5*(vp+V[i][jminusone][k]);

                //Vup
                if(k == NZ-1) vu = 0; //check for wall
                else vu = 0.5*(vp+V[i][j][k+1]);

                //Vdown
                if(k == 0) vd = 0; //check for wall
                else vd = 0.5*(vp+V[i][j][k-1]);

                //Ueast
                ue = 0.5*(U[i][j][k]+U[i][jplusone][k]);

                //Uwest
                uw = 0.5*(U[iminusone][j][k]+U[iminusone][jplusone][k]);

                //Wup
                if(k==NZ-1) wu=0; //check for wall
                else wu = 0.5*(W[i][j][k]+W[iplusone][j][k]);

                //Wdown
                if(k==0) wd =0; //check for wall
                else wd = 0.5*(W[i][j][k-1]+W[iplusone][j][k-1]);

                //Viscous Stresses
                dvdxe = 2*(ve-vp)/dx;
                dvdxw = 2*(vp-vw)/dx;
                dvdyn = 2*(vn-vp)/dy;
                dvdys = 2*(vp-vs)/dy;
                dvdzu = 2*(vu-vp)/dz;
                dvdzd = 2*(vp-vd)/dz;

                //Set Hny
                Hny[i][j][k] = (MU*(dvdxe-dvdxw)/(RHO*dx)+MU*(dvdyn-dvdys)/(RHO*dy)+MU*(dvdzu-dvdzd)/(RHO*dz))-((ue*ve-uw*vw)/dx+(vn*vn-vs*vs)/dy+(vu*wu-vd*wd)/dz);
            }
        }
    }
    return;
}

void calcHz()
{
    double wp, we, ww, wn, ws, wu, wd, ue, uw, vn, vs, dwdxe, dwdxw, dwdyn, dwdys, dwdzu, dwdzd;
    int iplusone, jplusone, iminusone, jminusone;

    for (int i=0;i<NX;i++)
    {
        for(int j=0;j<NY;j++)
        {
            for(int k=0;k<NZ-1;k++) //Only iterate to the NZ-2 Wcell
            {

                // Push Current Hwx into past
                Holdz[i][j][k]=Hnz[i][j][k];
                wp = W[i][j][k];

                //do periodicity in i,j directions (x,y)
                if (i==NX-1) iplusone = 0;
                    else iplusone = i+1;
                if (j==NY-1) jplusone = 0;
                    else jplusone = j+1;
                if (i==0) iminusone = 0;
                    else iminusone = i-1;
                if (j==0) jminusone = 0;
                    else jminusone = j-1;

                //Weast
                we = 0.5*(wp+W[iplusone][j][k]);

                //Wwest
                ww = 0.5*(wp+W[iminusone][j][k]);

                //Wnorth
                wn = 0.5*(wp+W[i][jplusone][k]);

                //Wsouth
                ws = 0.5*(wp+W[i][jminusone][k]);

                //Wup
                if(k == NZ-2) wu = 0.5*wp; //check for wall
                else wu = 0.5*(wp+W[i][j][k+1]);

                //Wdown
                if(k == 0) wd = 0.5*wp; //check for wall
                else wd = 0.5*(wp+W[i][j][k-1]);

                //Ueast
                ue = 0.5*(U[i][j][k]+U[i][j][k+1]);

                //Uwest
                uw = 0.5*(U[iminusone][j][k]+U[iminusone][j][k+1]);

                //Vnorth
                vn = 0.5*(V[i][j][k]+V[i][j][k+1]);

                //Vsouth
                vs = 0.5*(V[i][jminusone][k]+V[i][jminusone][k+1]);


                //Viscous Stresses
                dwdxe = 2*(we-wp)/dx;
                dwdxw = 2*(wp-ww)/dx;
                dwdyn = 2*(wn-wp)/dy;
                dwdys = 2*(wp-ws)/dy;
                dwdzu = 2*(wu-wp)/dz;
                dwdzd = 2*(wp-wd)/dz;

                //Set Hnz
                Hnz[i][j][k] = (MU*(dwdxe-dwdxw)/(RHO*dx)+MU*(dwdyn-dwdys)/(RHO*dy)+MU*(dwdzu-dwdzd)/(RHO*dz))-((ue*we-uw*ww)/dx+(vn*wn-vs*ws)/dy+(wu*wu-wd*wd)/dz);
            }
        }
    }
    return;
}

void calcUstar()
{
    double Tp, gravForceX,Pperiod;
    int iplusone;
    for (int i=0;i<NX;i++)
    {
        for(int j=0;j<NY;j++)
        {
            for(int k=0;k<NZ;k++)
            {
                if(i==NX-1)
                {
                    iplusone=0;
                    Pperiod=DPDX*LX;
                }
                else
                {
                    iplusone = i+1;
                    Pperiod=0;
                }
                Tp= 0.5*(T[iplusone][j][k]+T[i][j][k]);
                gravForceX = gravity[0]*RHO*(1.0-BETA)*(Tp-TREF);
                Ustar[i][j][k]=U[i][j][k]+DT*gravForceX+ (((-1.0*DT)/(RHO*dx))*(P[iplusone][j][k]-P[i][j][k])+Pperiod) + DT*(1.5*Hnx[i][j][k]-0.5*Holdx[i][j][k]);
            }
        }
    }
}

void calcVstar()
{
    double Tp, gravForceY;
    int jplusone;
    for (int i=0;i<NX;i++)
    {
        for(int j=0;j<NY;j++)
        {
            for(int k=0;k<NZ;k++)
            {
            if(j==NY-1) jplusone=0;
            else jplusone = j+1;
            Tp= 0.5*(T[i][jplusone][k]+T[i][j][k]);
            gravForceY = gravity[1]*RHO*(1.0-BETA)*(Tp-TREF);
             Vstar[i][j][k]=V[i][j][k]+DT*gravForceY+(((-1.0*DT)/(RHO*dy))*(P[i][jplusone][k]-P[i][j][k]))+DT*(1.5*Hny[i][j][k]-0.5*Holdy[i][j][k]);
            }
        }
    }
}

void calcWstar()
{
    double Tp, gravForceZ;
    for (int i=0;i<NX;i++)
    {
        for(int j=0;j<NY;j++)
        {
            for(int k=0;k<NZ-1;k++)//Only iterate to the NZ-2 Wcell
            {
             Tp= 0.5*(T[i][j][k+1]+T[i][j][k]);
             gravForceZ = gravity[2]*RHO*(1.0-BETA)*(Tp-TREF);
             Wstar[i][j][k]=W[i][j][k]+DT*gravForceZ+(((-1.0*DT)/(RHO*dz))*(P[i][j][k+1]-P[i][j][k]))+DT*(1.5*Hnz[i][j][k]-0.5*Holdz[i][j][k]);
            }
        }
    }
}

double calcUnew()
{
    double Tp, gravForceX,Pperiod,deltaU,l2Norm=0;
    int iplusone;
    for (int i=0;i<NX;i++)
    {
        for(int j=0;j<NY;j++)
        {
            for(int k=0;k<NZ;k++)
            {
                if(i==NX-1)
                {
                    iplusone=0;
                    Pperiod=DPDX*LX;
                }
                else
                {
                    iplusone = i+1;
                    Pperiod=0;
                }
                Tp= 0.5*(T[iplusone][j][k]+T[i][j][k]);
                gravForceX = gravity[0]*RHO*(1.0-BETA)*(Tp-TREF);
                deltaU=DT*gravForceX+(((-1.0*DT)/(RHO*dx))*(P[iplusone][j][k]-P[i][j][k])+Pperiod)+DT*(1.5*Hnx[i][j][k]-0.5*Holdx[i][j][k]);
                U[i][j][k]=U[i][j][k]+deltaU;
                l2Norm+=deltaU*deltaU;
            }
        }
    }
return sqrt(l2Norm);
}

double calcVnew()
{
    double Tp, gravForceY, deltaV,l2Norm=0;
    int jplusone;
    for (int i=0;i<NX;i++)
    {
        for(int j=0;j<NY;j++)
        {
            for(int k=0;k<NZ;k++)
            {
                if(j==NY-1) jplusone=0;
                else jplusone = j+1;
                Tp= 0.5*(T[i][jplusone][k]+T[i][j][k]);
                gravForceY = gravity[1]*RHO*(1.0-BETA)*(Tp-TREF);
                deltaV = DT*gravForceY+(((-1.0*DT)/(RHO*dy))*(P[i][jplusone][k]-P[i][j][k]))+DT*(1.5*Hny[i][j][k]-0.5*Holdy[i][j][k]);
                V[i][j][k]=V[i][j][k]+deltaV;
                l2Norm=deltaV*deltaV;
            }
        }
    }
    return sqrt(l2Norm);
}

double calcWnew()
{
    double Tp, gravForceZ,deltaW,l2Norm=0;
    for (int i=0;i<NX;i++)
    {
        for(int j=0;j<NY;j++)
        {
            for(int k=0;k<NZ-1;k++)//Only iterate to the NZ-2 Wcell
            {
             Tp= 0.5*(T[i][j][k+1]+T[i][j][k]);
             gravForceZ = gravity[2]*RHO*(1.0-BETA)*(Tp-TREF);
             deltaW=DT*gravForceZ+(((-1.0*DT)/(RHO*dz))*(P[i][j][k+1]-P[i][j][k]))+DT*(1.5*Hnz[i][j][k]-0.5*Holdz[i][j][k]);
             W[i][j][k]=W[i][j][k]+deltaW;
             l2Norm=deltaW*deltaW;
            }
        }
    }
    return sqrt(l2Norm);
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
                P[i][j][k]=0+DPDX*centroidLocs[i][j][k][0];
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

double checkMassBal()
{
    double localMassImbalance,vn,vs,ue,uw, wu,wd;
    double cumulativeMassImbalance=0;
    for (int i=0;i<NX;i++) //initialize Scalar Locations
    {
        for(int j=0;j<NY;j++)
        {
            for(int k=0;k<NZ;k++)
            {

                //West Face
                if(i==0) uw=U[NX-1][j][k];
                else     uw=U[i-1][j][k];

                //East Face
                ue=U[i][j][k];


                //South Face
                if(j==0) vs=V[i][NY-1][k];
                else     vs=V[i][j-1][k];

                //North face
                vn=Vstar[i][j][k];


                //Down Face
                if(k==0)  wd=0; //Wall => w=0 and dP/dz=0
                else wd=W[i][j][k-1];


                //Up Face
                if(k==NZ-1) wu=0;
                else wu=Wstar[i][j][k];

                localMassImbalance = dy*dz*(ue-uw)+dx*dz*(vn-vs)+dx*dy*(wu-wd);
                massbal[i][j][k]=localMassImbalance;
                cumulativeMassImbalance+=fabs(localMassImbalance);
            }
        }
    }
return cumulativeMassImbalance;
}








