#ifndef ULEMGRID

#ifndef std
#include <iostream>
#endif
#ifndef vector_h
#include <vector>
#endif
#include <string>
#include <fstream>
#include<cmath>
#include <stdio.h>
#include<sstream>

#include<stdexcept>
#include <algorithm>
//#include <omp.h>
void mprint(double);
void mprint(int);

#include <fstream>


#define ULEMGRID

class ulemgrid
{
public:
    ulemgrid(int,int);//constructor takes size of grid mxn as input params
    void set(std::string,int);
    void set(std::string,double);
    void set(std::string,std::vector<double>);
    void set(std::string,double*,int);
        std::vector<double> sinkareas;

    std::vector<double> get(std::string);
    
    void erosion_fluvial();//Fastscape method, erosion
    void lakefill();//Fill lakes
    void diffuse();//"Hillslope" diffusion
    void deposit(); //Deposition
    void erosion_fluvial2();//Temporary function for Testing explicit vs implicit methods
    
    
private:
    std::vector<double> U;
    std::vector<double> test;
    std::vector<double> kval;
    bool usefd=false;
    std::vector<double> ero;
    
    double L;
    double T=0;
    double m;
    int nstep;
    int nn;
    double dt;
    std::vector<int> stack;
    std::vector<int> undercapacity;
    double dy;
    std::vector<double> adds;
    std::vector<double> Z;
    int nx=1;
    int ny=1;
    double ks=.0001;
    double dx;
    int firstcall;
    std::vector<int> slpis ;
    double kd;
    double tt;
    double n_;
    double maxareasinkfill;
    std::vector<double> BC;
    std::vector<std::vector<int>> stackij;
    std::vector<int> idx;
    
    std::vector<int> ndons;
    int numiter;
    
    int nnn;
    double sk;
    
    
    std::vector<int> olength;
    int recursivestack(int,int,int);
    
    
    std::vector<int> BCX;
    
    std::vector<double> kw;
    int tstep;
    
    
    void erode();
    void basinfill(int,int);
    void fillls();
    void fillls2();
    
    void filll();
    void getdonors();
    
    void findsteepest();
    void createstack();
    void createstack2();
    std::vector<int> olist;
    
    void reset();
    std::vector<std::vector<int> > donors;
    
    void getacc();
    std::vector<double> accgrid;
    
    std::vector<double> dists;
    
    
    void landsed();
    
    
    class cmpr
    {
        std::vector<double> *dem;
    public:
        
        bool operator() (const int &i,const int &j) const
        {
            if ((*dem)[i]>(*dem)[j])
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    };
    void lakefill2();
    std::vector<int> sinkfill;
    std::vector<double> sed;
    bool uselandsed=false;
    std::vector<bool> usedr;
    void recursivesed(int);
    double sedextra;
    double area;
    std::vector<int>sedlist;
    std::vector<int> next;
    
    int nrecur;
    int nseds;
    std::vector<double> landsurf;
    
    std::vector<double> d,dsi;
    int mres=1000;
    
    void rangesearch(int,int);
    int rcount=1;
    double sangle=.000001;//shelf angle
    int rsi = 1;
    int minrs=1;
    int maxrs=1;
    int maxr=0;
    
    double critangle=10;//cone angle
    std::vector<int> catchments;
    int itriangle(int,int,int,int);
    std::vector<std::vector<int>> anglesx;
    std::vector<std::vector<int>> anglesy;
    std::vector<std::vector<int>> anglesx2;
    std::vector<std::vector<int>> anglesy2;
    std::vector<std::vector<double>> Dist;
    std::vector<std::vector<int>> I;
    std::vector<int> ic;
    std::vector<int> ic2;
    
    std::vector<std::vector<double>> Dist2;
    
    std::vector<double> anglespolar1;
    std::vector<double> anglespolar2;
    std::vector<int> ys;
    std::vector<int> xs;
    std::vector<int> ys2;
    std::vector<int> xs2;
    
    
    std::vector<int> anglesi2;
    std::vector<int> zni;
    std::vector<double>anglesz;
    std::vector<double>angleszi;
    std::vector<double> minz,anglesz2;
    std::vector<bool> nidx;
    int landorsea=0;//0 for land diffuse, 1 for sea
    
};


ulemgrid::ulemgrid(int i, int j)
{
    ny=i;
    nx=j;
    nn=i*j;
    
    ero.resize(nn+1);
    undercapacity.resize(nn+1);
    adds.resize(nn+1);
    Z.resize(nn+1);
        usedr.resize(nn+1);

    std::fill(adds.begin(),adds.end(),0.0);
    slpis.resize(nn+1);
    std::fill(slpis.begin(),slpis.end(),0);
    kval.resize(nn+1);
    stack.resize(nn+1);
    std::fill(stack.begin(),stack.end(),0);
    BCX.resize(nn+1);
    adds.resize(nn+1);
    donors.resize(nn+1,std::vector<int>(9));
    accgrid.resize(nn+1);
    std::fill(adds.begin(),adds.end(),0.0);
    dists.resize(nn+1);
    ndons.resize(nn+1);
    

    idx.resize(8);
    idx[0]=1;idx[1]=-1;idx[2]=ny;idx[3]=-ny;
    idx[4]=ny+1;
    idx[5]=ny-1;idx[6]=-ny+1;idx[7]=-ny-1;
    
    sinkfill.resize(nn+1);
    std::fill(sinkfill.begin(),sinkfill.end(),1);
    landsurf=Z;
    U.resize(nn+1);
    std::fill(U.begin(),U.end(),0);
    
    
    anglespolar1.resize(nn);
    anglespolar2.resize(nn);
    anglesz.resize(nn);
    
    anglesz2.resize(nn);
    
    
    
    angleszi.resize(nn);
    zni.resize(nn);
    anglesx.resize(20);
    anglesy.resize(20);
    ys.resize(nn+1);
    xs.resize(nn+1);
    ys2.resize(nn+1);
    xs2.resize(nn+1);
    anglesx2.resize(20);
    anglesy2.resize(20);
    sed.resize(nn+1);
    
    I.resize(20);
    Dist.resize(20);
    Dist2.resize(20);
    ic.resize(20);
    ic2.resize(20);
    anglesi2.resize(nn/1);
    minz.resize(nn+1);
    anglesz2.resize(nn);
    catchments.resize(nn+1);
    nidx.resize(nn+1);
    
    maxareasinkfill=0;
}

void ulemgrid::set(std::string nm,int val)
{
    
}
void ulemgrid::set(std::string nm,double val)
{
    if (nm.compare("firstcall")==0)
    {
        firstcall=val;
        std::cout<<val;
    }
    else if (nm.compare("kd")==0)
    {
        kd=val;
    }
    else if (nm.compare("m")==0)
    {
        m=val;
    }
    else if (nm.compare("n")==0)
    {
        n_=val;
    }
    
    else if (nm.compare("ks")==0)
    {
        ks=val;
    }
    else if (nm.compare("l")==0)
    {
        L=val;
    }
    else if (nm.compare("dt")==0)
    {
        dt=val;
    }
    else if (nm.compare("dx")==0)
    {
        dx=val;
    }
    else if (nm.compare("landsed") == 0)
    {
        uselandsed=(bool)val;
    }
    else if (nm.compare("dy")==0)
    {
        dy=val;
    }
    
    else if (nm.compare("seadiffusion")==0)
    {
        landorsea=val;
    }
    else if (nm.compare("maxareasinkfill")==0)
    {
        maxareasinkfill=val;
    }
    else
    {
        std::cout<<nm<<": parameter not found";
    }
}
void ulemgrid::set(std::string nm,double* val,int len)
{
    if (nm.compare("bc")==0)
    {
        BC.resize(len+1);
        for (int i=1;i<=len;i++)
        {
            BC[i]=val[i-1];
        }
        for (int i=0;i<len;i++)
        {
            BCX[val[i]]=1;
        }
        
    }
    else if (len==nn)
    {
        
        if (nm.compare("z")==0)
        {
            for (int i=1;i<nn+1;i++)
            {
                Z[i]=val[i-1];
            }
            
        }
        else if (nm.compare("k")==0)
        {
            for (int i=1;i<nn+1;i++)
            {
                kval[i]=val[i-1];
            }
            
        }
        else if (nm.compare("undercapacity")==0)
        {
            for (int i=1;i<nn+1;i++)
            {
                undercapacity[i]=(int)val[i-1];
            }
            
        }
        
        else if (nm.compare("u")==0)
        {
            for (int i=1;i<nn+1;i++)
            {
                U[i]=val[i-1];
            }
            
        }
        else
        {
            std::cout<<nm<<": parameter not found";
        }
    }
    else
    {
        throw std::invalid_argument( " Error in ulemgrid, input grid is wrong size" );
    }
    
    
}
void ulemgrid::set(std::string nm,std::vector<double> val)
{
    
    int len=val.size();
    
    if (nm.compare("bc")==0)
    {
        std::fill(BCX.begin(),BCX.end(),0);
        
        for (int i=0;i<len;i++)
        {
            BCX[val[i]]=1;
        }
        
    }
    else if (len==nn)
    {
        if (nm.compare("z")==0)
        {
            for (int i=1;i<nn+1;i++)
            {
                Z[i]=val[i-1];
            }
            
        }
        if (nm.compare("k")==0)
        {
            for (int i=1;i<nn+1;i++)
            {
                kval[i]=kval[i-1];
            }
            
        }
        if (nm.compare("undercapacity")==0)
        {
            
            for (int i=1;i<nn+1;i++)
            {
                undercapacity[i]=undercapacity[i-1];
            }
            
        }
        
        
    }
    else
    {
        throw std::invalid_argument( "Input grid is wrong size" );
    }
}



void ulemgrid::findsteepest()
{
    int ij;
    double slpxp;
    double slpxm;
    double slpyp;
    double slpym;
    double slpxpyp;
    double slpxmym;
    double slpxpym;
    double slpxmyp;
    double slpi = 1;
    double diag = std::pow(std::pow((double)dx,2.0)+std::pow((double)dy,2.0),.5);
    double maxslp;
    int redo=1;
    int nredo=0;
    int nsink=0;
    test.resize(nn+1);
#pragma omp simd
    for (int i=1;i<=nn;i++)
    {
        slpis[i]=i;
        dists[i]=100000000000000;
    }
#pragma omp parallel for private(slpyp,slpym,slpxm,slpxpyp,slpxmym,maxslp,ij) firstprivate(slpi,diag,ny,nx,dy,dx)
    for(int j=2;j<nx;j++)
        
    {
        for (int i=2;i <ny;i++)
        {
            
            ij = i+(j-1)*ny;
            if (BCX[ij]==1)
            {
                
                continue;
            }
            slpyp = (Z[ij]-Z[ij+1])/dy;
            slpym = (Z[ij]-Z[ij-1])/dy;
            slpxp = (Z[ij]-Z[ij+ny])/dx;
            slpxm = (Z[ij] -Z[ij-ny])/dx;
            slpxpyp = (Z[ij]-Z[ij+1+ny])/diag;
            slpxpym = (Z[ij]-Z[ij-1+ny])/diag;
            slpxmyp = (Z[ij]-Z[ij+1-ny])/diag;
            slpxmym = (Z[ij]-Z[ij-1-ny])/diag;
            maxslp = 0;
            slpi = 0;
            
            if (slpyp>maxslp)
            {
                maxslp = slpyp;
                slpi = 1;
                dists[ij] = dy;
            };
            
            if (slpym>maxslp)
            {
                maxslp = slpym;
                slpi = -1;
                dists[ij]=dy;
            };
            if(slpxp>maxslp)
                
            {maxslp = slpxp;
             slpi = ny;
             dists[ij]=dx;
            };
            if (slpxm>maxslp)
            {
                maxslp = slpxm;
                slpi = -ny;
                dists[ij]=dx;
            };
            if (slpxpyp>maxslp)
            {
                maxslp = slpxpyp;
                slpi = ny+1;
                dists[ij]=diag;
            };
            if (slpxpym>maxslp)
            {
                maxslp = slpxpym;
                slpi = ny-1;
                dists[ij]=diag;
            };
            if(slpxmyp>maxslp)
            {
                maxslp = slpxmyp;
                slpi = -ny+1;
                dists[ij]=diag;
            };
            if(slpxmym>maxslp)
            {
                maxslp = slpxmym;
                slpi = -ny-1;
                dists[ij]=diag;
                
            };
            if (maxslp==0)
            {
                nsink++;
            }
            
            slpis[ij] = slpi+ij;
            
        };
    };
    std::cout<<nsink<<std::endl;
    //mexPrintf("%d\n",nsink);
};

void ulemgrid::getdonors()
{
#pragma omp parallel for
    for (int ij=1;ij<=nn;ij++)
    {
        if (slpis[ij]!=ij)
        {
            ndons.at(slpis[ij])++;
            donors[slpis[ij]][ndons[slpis[ij]]]=ij;
        }
    }
};

void ulemgrid::createstack()
{
    olist.resize(1);
    
#pragma omp simd
    for (int ij=1;ij<=nn;ij++)
    {
        if (slpis[ij]==ij)
        {
            olist.push_back(ij);
        }
    }
    stackij.resize(olist.size());
#pragma omp parallel for schedule(dynamic)
    for (int i=1;i<olist.size();i++)
    {
        int d=1;
        int p=1;
        int ndon=0;
        int ij=olist[i];
        
        //mexPrintf("%d\n",olist[i]);
        
        stackij[i].push_back(ij);
        p++;
        ndon++;
        
        for (int ik=1;ik<=ndons[ij];ik++)
        {
            numiter=0;
            d = donors[ij][ik];
            p = recursivestack(d,p,i);
        };
    };
    olength.resize(olist.size());
#pragma omp parallel for
    for (int i=1;i<olist.size();i++)
    {
        olength[i]=stackij[i].size();
    }
    
    nnn=olist.size();
    int cc=1;
#pragma omp simd
    
    for (int j=1;j<nnn;j++)
    {
        
        for (int i=0; i<olength[j];i++)
        {
            stack[cc]=stackij[j][i];
            cc++;
        }
    }
    // mexPrintf("%d\n",stack[0]);
};
int ulemgrid::recursivestack(int r,int p,int i)
{
    
    if (p<=nn&&numiter<1e5)
    {
        int d;
        stackij[i].push_back(r);
        p++;
#pragma omp simd
        for (int ij=1;ij<=ndons[r];ij++)
        {
            d = donors[r][ij];
            p = recursivestack(d,p,i);
        }
        return p;
    }
    else
    {
        
        return 1;
    }
    
};
void ulemgrid::createstack2()
{
    olist.resize(1);
    
#pragma omp simd
    for (int ij=1;ij<=nn;ij++)
    {
        if (slpis[ij]==ij)
        {
            olist.push_back(ij);
        }
    }
    
    stackij.resize(olist.size());
    
#pragma omp parallel for schedule(dynamic)
    for (int i=1;i<olist.size();i++)
    {
        int d=1;
        int p=1;
        int ndon=0;
        int ij=olist[i];
        
        //mexPrintf("%d\n",olist[i]);
        
        stackij[i].push_back(ij);
        
        for (int u=0;u<=ndon;u++)
        {
            
            ij=stackij[i][u];
            
            for (int ik=1;ik<=ndons.at(ij);ik++)
            {
                ndon++;
                d = donors[ij][ik];
                stackij[i].push_back(d);
                
            };
        };
    };
    olength.resize(olist.size());
#pragma omp parallel for
    for (int i=1;i<olist.size();i++)
    {
        olength[i]=stackij[i].size();
    }
    nnn=olist.size();
    int cc=1;
#pragma omp simd
    for (int j=1;j<nnn;j++)
    {
        
        for (int i=0; i<olength[j];i++)
        {
            stack[cc]=stackij[j][i];
            cc++;
        }
    }
    // mexPrintf("%d\n",stack[0]);
};

void ulemgrid::getacc()
{
    
#pragma omp simd
    for (int i=1;i<=nn;i++)
    {
        accgrid[i]=1.0;
    }
    int cc;
#pragma omp parallel for
    for (int j=1;j<nnn;j++)
    {
#pragma omp simd
        for (int i=olength[j]-1;i>0;i--)
        {
            
            accgrid[slpis[stackij[j][i]]]=accgrid[slpis[stackij[j][i]]]+accgrid[stackij[j][i]];
        }
    }
    
}

void ulemgrid::erode()
{
    double f;
    double fs=0;
    
    double x;
    int ni=1,ni2=1;
    if (n_>1)
    {
        ni=1;
        ni2=5;
    }
    if (n_<1)
    {
        ni=1;
        ni2=5;
    }
    double f2=0;
    std::vector<double> Zi;
    std::vector<double>Zii;
    double Zi1;
    int first=1;
    Zi=Z;
    double sp1 = std::pow(dx*dy,m)*dt;
    //std::cout<<std::endl<<"dt= "<<dt<<std::endl;
    if (n_>=1)
    {
        int ncs=0;
        for (int l=1;l<=ni;l++)
        {
            
#pragma omp parallel for private(f) firstprivate(x,ncs,first) schedule(dynamic)
            for (int j=1;j<nnn;j++)
            {
#pragma omp simd
                for (int i=1;i<olength[j];i++)
                {
                    if (undercapacity[stackij[j][i]]==0)
                    {
                        if (Z[stackij[j][i]]>=0)
                        {
                            
                            f = kval[stackij[j][i]]/std::pow(dists[stackij[j][i]],n_) *
                                    std::pow(accgrid[stackij[j][i]],m)*sp1*
                                    std::pow(Zi[stackij[j][i]]-std::max(0.0,Z[slpis[stackij[j][i]]]),n_-1.0);
                            f2=0;//dt/(Zi[stackij[j][i]]-std::max(0.0,Z[slpis[stackij[j][i]]])); Gives weird result, probably wrong
                        }
                        x=1;
                        
                        for (int ll=1;ll<=5;ll++)
                        {
                            
                            
                            
                            x=(x-(x-1+f*std::pow(x,n_))/(1+(n_)*f*std::pow(x,n_-1.0)));
                        }
                        
                        Zi1=Z[stackij[j][i]];
                        if (Z[stackij[j][i]]>=0)
                        {///Z[stackij[j][i]]=(f*Z[slpis[stackij[j][i]]]+Z[stackij[j][i]])/(1+f);
                            Z[stackij[j][i]] = std::max(0.0,Z[slpis[stackij[j][i]]]) +x*
                                    std::max(0.0,(Zi[stackij[j][i]]-std::max(0.0,Z[slpis[stackij[j][i]]])));
                            
                        }
                        // Z[stackij[j][i]]=(Zi[stackij[j][i]]+f*Z[slpis[stackij[j][i]]]+dt*U[stackij[j][i]])/(1+f);}
                        if (kval[stackij[j][i]]>0)
                        {
                            ero[stackij[j][i]]=Zi1-Z[stackij[j][i]];
                        }
                    }
                }
                
                
            }
        }
    }
}

void ulemgrid::reset()

{

    if (usefd==false)
    {
        stackij.clear();
    }
    std::fill(accgrid.begin(),accgrid.end(),0);
    std::fill(sed.begin(),sed.end(),0);
    
    std::fill(ndons.begin(),ndons.end(),0);
    std::fill(BCX.begin(),BCX.end(),0);

    std::fill(ero.begin(),ero.end(),0.0);
#pragma omp parallel for
    for (int i=1;i<BC.size();i++)
    {
        
        BCX[BC[i]]=1;
    }

    for (int i=1;i<=BC.size();i++)
    {
        // std::cout<<i<<" ";
        BCX[BC[i]]=1;
    }
    
    
    
    
}


void ulemgrid::filll()
{
    if (stackij.size()>0)
    {
        
#pragma omp parallel for
        for (int j=1;j<nnn;j++)
        {
#pragma omp simd
            for (int i=0;i<olength[j];i++)
            {
                if (Z[stackij[j][i]]<=Z[slpis[stackij[j][i]]])
                {
                    Z[stackij[j][i]]=Z[slpis[stackij[j][i]]]+1e-6;
                }
            }
        }
    }
}

void ulemgrid::fillls()
{
    std::vector<double> cero=ero;
    std::vector<double> Zn=Z;
    double Lt;
    double ssum=0;int nj=0;double sedsum;
    double dti=dt;
    double A2;
    double di;
    double f;
    double d2;
    std::vector<double> Zi=Z;
    
    if (stackij.size()>0)
    {
        
#pragma omp parallel for schedule(dynamic) firstprivate(dti,ssum,sedsum,nj,dt,dx,dy,nn,ks,L,T)
        
        for (int j=1;j<nnn;j++)
        {
            
            double erol;
            
            double lsum;
            lsum=1;
            int y;
            int eroli=1;
            int dts=1;
            dt=dti;
            int y2=1;
            
            int nt1=1;
            
            while (nt1<=dts)
            {
                y=1;
                nt1++;
                
                
                sedsum=1.00000001e-3;
                
                while (sedsum>1e-3)
                {
                    if (y==50)
                    {
                        dts*=2;
                        dt/=2;
                        for (int i=1;i<olength[j];i++)
                        {
                            ero[stackij[j][i]]/=2;
                        }
                        y=1;
                    }
                    y2++;
                    if (y2>5000)
                    {
                        break;
                        std::cout<<"Not Converging";
                    }
                    
                    
                    lsum=sedsum;
                    
                    sedsum=0;
                    
                    y++;
                    
                    for (int i=olength[j]-1;i>=0;i--)
                    {
                        Zn[stackij[j][i]]=Z[stackij[j][i]];
                        
                        cero[stackij[j][i]]=ero[stackij[j][i]];
                        
                    }
                    
                    Zn[slpis[stackij[j][0]]]=Z[slpis[stackij[j][0]]];
                    cero[slpis[stackij[j][0]]]=ero[slpis[stackij[j][0]]];
                    for (int i=olength[j]-1;i>0;i--)
                    {
                        if (stackij[j][i]!=slpis[stackij[j][i]])
                        {
                            
                            cero[slpis[stackij[j][i]]]+=cero[stackij[j][i]];
                            
                        }
                    }
                    
                    
                    double tl;
                    if (m==1)
                    {
                        
                        f= dt*ks;
                    }
                    else
                    {
                        f=dt*kval[stackij[j][1]]*std::pow(dx*dy,m);
                    }
                    
                    double A;
                    
                    for (int i=1;i<=olength[j]-1;i++)
                    {
                        
                        if (undercapacity[stackij[j][i]]>=1&&Zn[stackij[j][i]]>=0)
                        {
                            
                            A=std::pow(accgrid[stackij[j][i]],m);
                            
                            if (m==1)//If uncercapacity
                            {
                                Lt=L;
                                A2=1.0;
                                di=dists[stackij[j][i]];
                                d2=1.0;
                            }
                            else//If yuan et al
                            {
                                Lt=1.0;
                                A2=accgrid[stackij[j][i]];
                                di=1.0;
                                d2=std::pow(dists[stackij[j][i]],n_);
                            }
                            
           
                            tl=((cero[stackij[j][i]]-ero[stackij[j][i]])/L/A2*di+Zn[stackij[j][i]]+
                                    std::max(0.0,Zn[slpis[stackij[j][i]]])*
                                    f/Lt/d2*A)/(f/Lt/d2*A+1.0);
                            
          
                            
                            erol=ero[stackij[j][i]];
                            ero[stackij[j][i]]=(Zn[stackij[j][i]]-tl);
                            
                            Zn[stackij[j][i]]=tl;
                            
                            if (sedsum<std::abs(erol-ero[stackij[j][i]]))
                            {sedsum=std::abs(erol-ero[stackij[j][i]]);eroli=i;}
                            
                        }
                        
                        
                    }
                }
                
                for (int i=0;i<=olength[j]-1;i++)
                {
                    if (undercapacity[stackij[j][i]]>=1)
                    {Z[stackij[j][i]]=Zn[stackij[j][i]];
                     
                    }
                }
            }
            
            
            
            {ssum+=cero.at(stackij[j][0]);
            }
            
        }


    }
    dt=dti;
    
    
    
}

void ulemgrid::fillls2()
{
    std::vector<double> cero=ero;
    std::vector<double> Zn=Z;
    
    double ssum=0;int nj=0;double sedsum;
    double dti=dt;
    
    if (stackij.size()>0)
    {
        int dts=1;
        
        dts=200;
        dt/=200;
        int no=1;
        std::vector<double> cerol;
        cerol.resize(nn+1);
#pragma omp parallel for schedule(dynamic) firstprivate(dti,ssum,sedsum,nj,dt,dx,dy,nn,ks,L,T)
        
        for (int j=1;j<nnn;j++)
        {
            double erol;
            
            double lsum;
            lsum=1;
            int y;
            int eroli=1;
            int y2=1;
            
            int nt1=1;
            
            
            while (nt1<dts)
            {
                y=1;
                nt1++;
                
                sedsum=1.0001e-3;
                
                
                {
                    
                    
                    lsum=sedsum;
                    
                    sedsum=0;
                    
                    y++;
                    
                    
                    Zn[slpis[stackij[j][0]]]=Z[slpis[stackij[j][0]]];
                    
                    
                    double tl;
                    double f= dt*ks;
                    double A;
                    for (int i=olength[j]-1;i>0;i--)
                    {
                        cerol[stackij[j][i]]=cero[stackij[j][i]];
                        cerol[slpis[stackij[j][i]]]=cero[slpis[stackij[j][i]]];
                        
                        cero[stackij[j][i]]=0;
                        cero[slpis[stackij[j][i]]]=0;
                        
                    }
                    cero[stackij[j][olength[j]-1]]=0;
                    cero[stackij[j][0]]=0;
                    for (int i=olength[j]-1;i>0;i--)
                    {
                        if (nt1==2)
                            
                        {                        no++;
                        }
                        if (undercapacity[stackij[j][i]]>=1)
                        {
                            
                            A=std::pow(accgrid[stackij[j][i]],1.0);
                            
                            tl=1/L*dists[stackij[j][i]]*(cero[stackij[j][i]]-
                                    ks*dt*accgrid[stackij[j][i]]*
                                    (Z[stackij[j][i]]-std::max(0.0,Z[slpis[stackij[j][i]]]))/dists[stackij[j][i]]);
                            
                            erol=ero[stackij[j][i]];
                            
                            Zn[stackij[j][i]]+=tl;
                            cero[stackij[j][i]]-=tl;
                            cero[slpis[stackij[j][i]]]+=cero[stackij[j][i]];
                            
                            
                            if (sedsum<std::abs(erol-ero[stackij[j][i]]))
                            {sedsum=std::abs(erol-ero[stackij[j][i]]);eroli=i;}
                        }
                        
                        
                    }
                }
                for (int i=1;i<=olength[j]-1;i++)
                {
                    if (undercapacity[stackij[j][i]]>=1)
                    {Z[stackij[j][i]]=Zn[stackij[j][i]];
                     
                    }
                }
            }
            
            
            
            {ssum+=cero.at(stackij[j][0]);
            }
        }
        
        
    }
    dt=dti;
    
    
}

void ulemgrid::basinfill(int ij,int ij1)
{
    double min=Z.at(ij);
    for (int i=0;i<8;i++)
    {
        if (Z[ij+idx[i]]<min)
        {
            min=Z[ij+idx[i]];
        }
        
    }
    if (Z[ij]<=min)
    {
        Z[ij]=min+.01;
        adds[ij1]+=min+.01-Z[ij];
        
        basinfill(ij,ij1);
        for (int i=0;i<8;i++)
        {
            if (BCX[ij+idx[i]]==0)
            {basinfill(ij+idx[i],ij1);}
        }
    }
}
void ulemgrid::erosion_fluvial()
{    

    reset();
        

    if (usefd==false)
    {
        findsteepest();
            

        getdonors();
            
        createstack2();
        
    }
    getacc();
    
    erode();
    std::cout<<ks<<std::endl;
    if (ks>1e-100)
    {
        fillls();
    }
    for (int i=1;i<=ny;i++)
    {
        
        if (BCX[i]==0)
        {
            Z[i]=Z[i+ny];
            ero[i]=0;
        }
    }
#pragma omp parallel for
    for (int i=nn-ny+1;i<=nn;i++)
    {
        if (BCX[i]==0)
        {
            Z[i]=Z[i-ny];
            ero[i]=0;
        }
        
    }
#pragma omp parallel for
    for (int i=1;i<=nn;i+=ny)
    {
        if (BCX[i]==0)
        {
            Z[i]=Z[i+1];
            ero[i]=0;
        }
        
    }
    for (int i=ny;i<=nn;i+=ny)
    {
        if (BCX[i]==0)
        {
            Z[i]=Z[i-1];
            ero[i]=0;
        }
    }
    
}

void ulemgrid::erosion_fluvial2()
{
    reset();
    if (usefd==false)
    {
        findsteepest();

        getdonors();

        createstack2();

    }
    getacc();
    
    erode();
    if (ks>1e-100)
    {
        fillls2();
    }
        std::cout<<'here5';

}

std::vector<double> ulemgrid::get(std::string nm)
{
    std::vector<double> val;
    val.resize(nn);
    if (nm.compare("bc")==0)
    {
        for (int i=1;i<nn;i++)
        {
            val[i-1]=BCX[i];
        }
    }
    else if (nm.compare("acc")==0)
    {
        for (int i=1;i<nn;i++)
        {
            val[i-1]=accgrid[i];
        }
    }
    else if (nm.compare("z")==0)
    {
        for (int i=1;i<nn+1;i++)
        {
            val[i-1]=Z[i];
        }
        
    }
    else if (nm.compare("sinkareas")==0)
    {
        for (int i=1;i<nn+1;i++)
        {
            val[i-1]=sinkareas[i];
        }
    }

    else if (nm.compare("k")==0)
    {
        for (int i=1;i<nn+1;i++)
        {
            val[i-1]=kval[i];
        }
        
    }
    else if (nm.compare("test")==0)
    {
        for (int i=1;i<nn+1;i++)
        {
            val[i-1]=test[i];
        }
        
    }
    else if (nm.compare("ero")==0)
    {
        for (int i=1;i<nn+1;i++)
        {
            val[i-1]=ero[i];
        }
        
    }
    else if (nm.compare("stack")==0)
    {
         for (int i=1;i<nn+1;i++)
        {
            val[i-1]=stack[i];
        }
    }
    else if (nm.compare("rec")==0)
    {
         for (int i=1;i<nn+1;i++)
        {
            val[i-1]=slpis[i];
        }
    }
    else
    {std::cout<< "Parameter not found";}
    
    return val;
}
#ifndef PQ
#define PQ
class priorityq
{
public:
    int queue;
    priorityq(std::vector<double>&);
    std::vector<int> left;
    std::vector<double> z;
    void push(int);
    int pop();
    int top();
    int numel=0;
    
private:
    int nn;
    int tmp;
    int u;
    int uu;
    int ul;
    
};

#endif
priorityq::priorityq(std::vector<double>& grid)
{
    nn=grid.size();
    left.resize(nn+1);
    z=grid;
}

int priorityq::top()
{
    return left[1];
    
}
int priorityq::pop()
{
    
    
    uu =left[1];
    
    left[1]=left[numel];
    left[numel]=0;
    u=2;
    ul=u/2;
    while (u<numel-2)
    {
        if (z[left[u]]<z[left[u+1]])
        {
            
            if (z[left[ul]]>z[left[u]])
            {
                
                std::swap(left[ul],left[u]);
                
                ul=u;
                u=2*u;
                
                
            }
            else
            {
                break;
            }
        }
        else if(z[left[ul]]>z[left[u+1]])
        {
            std::swap(left[ul],left[u+1]);
            
            u=2*(u+1);
            ul=u/2;
        }
        else
        {
            break;
        }
        
        
    }
    numel--;
    
    return uu;
    
}



void priorityq::push(int i)

{
    numel++;
    
    
    u=numel;
    ul=u/2;
    
    left[u]=i;
    while (ul>0)
    {
        
        if (z[left[ul]]>z[left[u]])
        {
            
            std::swap(left[ul],left[u]);
            
            
            
        }
        else
        {
            break;
        }
        u=u/2;
        ul=u/2;
        
        
        
        
    }
    
}

void ulemgrid::lakefill2()
{
    double pittop;
    int c=0;
    int p = 0;
    std::vector<bool> closed;
    std::vector<int> pit;
    closed.resize(nn+1);
    //std::vector<int> pq;
    //pq.resize(nn+1);
    //std:fill(pq.begin(),pq.end(),-9999);
    pit.resize(nn+1);
    // typedef std::priority_queue<int,std::vector<int>,cmpr> mpr;
    // mpr open ((cmpr(&grid)));
    priorityq open(Z);
    if (!uselandsed)
    {
        
        for (int i=0;i<BC.size();i++)
        {
            open.push(BC[i]);
            closed[BC[i]]=true;
            c++;
        }
    }
    
    for (int i=1;i<=ny;i++)
    {
        if (closed[i]==false)
        {
            closed[i]=true;
            
            open.push(i);
            
            c++;
        }
        
    }
    for (int i=nn-ny+1;i<=nn;i++)
    {
        if (closed[i]==false)
        {
            closed[i]=true;
            
            open.push(i);
            
            c++;
        }
    }
    for (int i=ny;i<=nn;i+=ny)
    {
        if (closed[i]==false)
        {
            closed[i]=true;
            
            open.push(i);
            
            c++;
        }
    }
    for (int i=ny+1;i<=nn;i+=ny)
    {
        if (closed[i]==false)
        {
            closed[i]=true;
            
            open.push(i);
            
            c++;
        }
    }
    int s;
    int si=1;
    int ij;
    int ii;
    int jj;
    int ci=1;
    pittop=-9999;
    while (c>0||p>0)
    {
        if (p>0&&c>0&&pit[p-1]==-9999)//pq[c-1]==pit[p-1])
        {
            //s=open.top();
            s=open.pop();
            c--;
            pittop=-9999;
        }
        else if (p>0)
        {
            
            s=pit[p-1];
            pit[p-1]=-9999;
            p--;
            if (pittop==-9999)
            {
                pittop=Z[s];
            }
        }
        else
        {
            //s=open.top();
            s=open.pop();
            c--;
            pittop=-9999;
        }
        
        //  #pragma omp parallel for
        
        for (int i=0;i<=7;i++)
        {
            
            ij = idx[i]+s;
            
            ii= (ij-1)%ny+1;
            jj = (int)((ij-1)/ny)+1;
            if ((ii>=1)&&(jj>=1)&&(ii<=ny)&&(jj<=nx)&&closed[ij]==false)//)
            {
                
                closed[ij]=true;
                
                if (Z[ij]<=Z[s])
                {
                    
                    Z[ij]=Z[s]+1e-8;
                    
                    pit[p]=ij;
                    p++;
                    
                }
                
                else
                {
                    open.push(ij);
                    c++;
                }
                
            }
        }
    }
    
    
    
}


void ulemgrid::landsed()
{
    if (maxareasinkfill>0)
    {
        sinkareas.resize(nn+1);
        std::fill(sinkareas.begin(),sinkareas.end(),0.0);
    }

    sed = ero;
    if (sed.size()<nn){sed.resize(nn+1);}
    std::vector<double> sedsub;
    sedsub.resize(nn+1);
    next.resize(nn+1);
    std::fill(sedsub.begin(),sedsub.end(),0);
    sedlist.resize(nn+1);
    std::fill(usedr.begin(),usedr.end(),false);

    double tsed=0;double tsed2=0;
    for (int i=1;i<=ny;i++)
    {
        usedr[i]=true;
    }
    
    for (int i=nn-ny+1;i<=nn;i++)
    {
        usedr[i]=true;
    }
    
    for (int i=ny;i<=nn;i+=ny)
    {
        usedr[i]=true;
    }
    for (int i=ny+1;i<=nn;i+=ny)
    {
        usedr[i]=true;
    }
    //build cumulative sediment from ero and stream network
    //zero sediment everywhere but outlet to not count it twice
    //#pragma omp parallel for

    for (int i =nn;i>=1;i--)
    {
        tsed+=ero[stack[i]];
        
        sed[slpis[stack[i]]] +=sed[stack[i]];
    }
    int nsink=0;
    //#pragma omp parallel for
    for (int i =1;i<=nn;i++)
    {
        if (slpis[stack[i]]!=stack[i])
        {
            
            sed[stack[i]]=0;
            
        }
        else
        {
            nsink++;
        }
    }

    for (int i =1;i<=nn;i++)
    {
        if( (landsurf[i]<Z[i])&&(usedr[i])==false)
        {
            
            nseds=0;
            sedextra=0;
            area=0;
            nrecur=0;
            //std::fill(sedlist.begin(),sedlist.end(),0);
            recursivesed(i);
            double fact;
            if (sedextra>=area)
            {
                fact=1;
            }
            else if (area>0)
            {
                fact = sedextra/area;
            }
            
            else
            {
                fact=0;
            }
            
            double addsed=0;
            for (int j=1;j<=nseds;j++)
            {
                if (sedlist[j]>nn||sedlist[j]<0)
                {std::cout<<"warning: bad access in recursive sedfill"<<std::endl;
                 break;
                }
                
                addsed = (Z[sedlist[j]]-landsurf[sedlist[j]])*fact;
                landsurf[sedlist[j]]+=addsed;
                if (maxareasinkfill>0)
                {
                    sinkareas[sedlist[j]]=area;
                }
                tsed2+=addsed;
                
            }
            
        }
    }
    if (maxareasinkfill>0)
    {
        for (int i=1;i<nn+1;i++)
        {
            if (sinkareas[i]<maxareasinkfill)
            {
                landsurf[i] = Z[i];
            }
        }
    }
    Z=landsurf;
}
void ulemgrid::recursivesed(int ij)
{
    usedr[ij]=true;
    int c=0;
    bool go = true;
    while (go)
    {
        
        nrecur++;
        nseds++;
        
        sedlist[nseds]=ij;
        for (int i=0;i<8;i++)
        {
            
            if( (landsurf[idx[i]+ij]<Z[idx[i]+ij] )&&(usedr[idx[i]+ij]==false))
            {
                next[c]=idx[i]+ij;
                usedr[ij+idx[i]]=true;
                
                c++;
            }
            
        }
        
        area+=(Z[ij]-landsurf[ij]);
        sedextra+=sed[ij];
        c--;
        
        
        
        if (c>=0)
        {
            ij=next[c];
            
        }
        else
        {
            go=false;
        }
        
    }
}
void ulemgrid::lakefill()
{
    if (uselandsed==1)
    {
        uselandsed=true;
        landsurf=Z;
        
        lakefill2();
        
        landsed();
    }
    else
    {
        uselandsed=false;
        landsurf = Z;
        lakefill2();
        std::cout<<"ere";
        landsurf = Z;
    }
}
void ulemgrid::deposit()
{
    double sumsedij=0;
    std::vector<double> d,dsi;
    int ij;
    double sumer=0;
    //compute sediment using stack and erosion
    sed=ero;
    for (int i = nn;i>=1;i--)
    {
        if (slpis[stack[i]]!=stack[i])
        {
            
            
            sumer+=ero[stack[i]];
            sed[slpis[stack[i]]]+=sed[stack[i]];
        }
        
    }
    dsi.resize(nn+1);
    int cat=2;
    std::vector<bool> ends;
    ends.resize(nn+1);
    std::vector<double> circsurf;
    std::vector<double> seds;
    
    
    seds.resize(nn+1);
    catchments.resize(nn+1);
    int tic;
    int ri=1;
    //compute catchments
    std::fill(catchments.begin(),catchments.end(),1);
    
    for (int i = 1;i<=nn;i++)
    {
        //int i=1;
        
        
        catchments[stack[i]]=catchments[slpis[stack[i]]];
        if (catchments[stack[i]]==1)
        {
            
            catchments[stack[i]]=cat;
            ends[stack[i]]=true;
            dsi[slpis[stack[i]]]=sed[slpis[stack[i]]];
            
            cat++;
        }
    }
    //std::cout<<sumsedij<<std::endl;
    int rs;
    double lastsed=1e10;
    std::cout<<"nx= "<< nx << std::endl;
    std::cout<<"ny= "<< ny << std::endl;
    for (int i=1;i<=ny;i++)
    {

        for (int j=1;j<=nx;j++)
        {
            ij = i+(j-1)*ny;
            if (dsi.at(ij)>0&&(((i>1&&i<ny&&j>1&&j<nx))&&Z[ij]<0))
            {
                rs=minrs;
                test[ij]=1;
                double extrased=1;
                int rsii=rsi;
                double bslvl=0;
                ri=0;
                while (extrased>0)
                {
                    
                    tic=itriangle(rs,ri,j,i);
                    ri++;
                    
                    for (int l=1;l<tic;l++)
                    {
                        
                        if (!nidx[l])
                        {
                            
                            if (anglesz[l]>0)
                            {
                                anglesz[l]=0;
                            }
                        }
                    }
                    
                    double sumseds=0;
                    double circsurf;
                    
                    
                    
                    for (int l=1;l<tic;l++)
                    {
                        if (!nidx[l])
                        {
                            circsurf=anglesz[l];
                            seds[l]=circsurf-Z[zni[l]];
                            seds[l]=(seds[l]+std::fabs(seds[l]))/2.0;
                            sumseds+=seds[l];
                            
                        }
                    }
                    
                    extrased=sed[ij]-sumseds;
                    extrased=(extrased+std::fabs(extrased))/2.0;
                    
                    
                    
                    if (extrased==0)
                    {
                        
                        double sedx;
                        for (int l=1;l<tic;l++)
                        {
                            if (!nidx[l])
                            {
                                sedx=sed.at(ij)/sumseds*seds.at(l);
                                Z.at(zni.at(l))=Z[zni[l]]+sedx;
                                sumsedij+=sed[ij]/sumseds*sedx;;
                                
                            }
                        }
                    }
                    
                    rs=rs+rsii;
                    rsii=rsii+(int)(rsii/2+1);
                    
                    
                    if ((3.2*rs*rs>nn&&extrased>0))
                    {
                        double sedx;
                        
                        for (int l=1;l<tic;l++)
                        {
                            
                            if (!nidx[l])
                            {
                                
                                
                                sedx=seds[l];
                                
                                if (sumseds<=sed[ij])
                                {
                                    Z[zni[l]]=Z[zni[l]]+sedx;
                                }
                                else
                                {
                                    Z[zni[l]]=Z[zni[l]]+sed[ij]/sumseds*sedx;
                                }
                                
                            }
                        }
                        
                        break;
                    }
                    lastsed=extrased;
                    
                }
                
                
                
                
            }
            
        }
    }
    
}




int ulemgrid::itriangle(int r,int ri, int xi, int yi)
{
    double I1;
    
    int id = (int)(r/rsi);
    //if geometry of given radius has not previously been solved for...
    
    if (ri>maxr)
    {
        
        //x,y locations of area within circle
        anglesx[ri].resize(4*r*r+4);
        anglesy[ri].resize(4*r*r+4);
        
        //x,y locations on edge of circle
        anglesx2[ri].resize(4*r*r+4);
        anglesy2[ri].resize(4*r*r+4);
        
        
        I[ri].resize(4*r*r+4);
        Dist[ri].resize(4*r*r+4);
        Dist2[ri].resize(4*r*r+4);
        
        //Find cells within circle and on edge
        rangesearch(r,ri);
        
        
        //Find the theta value of cells within circle
        for (int i=1;i<ic[ri];i++)
        {
            
            anglespolar1[i] = std::atan2((double)anglesy[ri][i]-(double)mres/2,(double)anglesx[ri][i]-(double)mres/2);
            
        }
        
        //Find the theta value of cells on edge
        for (int i=1;i<ic2[ri];i++)
        {
            anglespolar2[i] = std::atan2((double)anglesy2[ri][i]-(double)mres/2,(double)anglesx2[ri][i]-(double)mres/2);
            
        }
        
        //Find the closest point on edge of circle, to that within the circle
        //based on theta value
        for (int i=1;i<ic[ri];i++)
        {
            double mI1=9999999;
            for (int j=1;j<ic2[ri];j++)
            {
                
                I1=std::fabs((double) anglespolar1[i]-(double)anglespolar2[j]);
                
                if (I1<mI1)
                {
                    mI1=I1;
                    I[ri][i]=j;
                }
                
            }
        }
        
        //calculate the distance from each interior point to closest edge point
        for (int i=1;i<ic[ri];i++)
        {
            Dist[ri][i]=std::sqrt(std::pow((double)anglesy2[ri][I[ri][i]]-
                    (double)anglesy[ri][i],2.0)+
                    std::pow((double)anglesx2[ri][I[ri][i]]-
                    (double)anglesx[ri][i],2.0));
        }
        
        //make dimension-independent
        for (int i=1;i<ic2[ri];i++)
        {
            anglesy2[ri][i]-=(int)mres/2;
            anglesx2[ri][i]-=(int)mres/2;
        }
        //Not sure .
        for (int i=1;i<ic[ri];i++)
        {
            Dist2[ri][i]=std::sqrt(std::pow(anglesy2[ri][I[ri][i]],2.0) +std::pow(anglesx2[ri][I[ri][i]],2.0));
            
        }
        
        //Make dimension independent
        for (int i=1;i<ic[ri];i++)
        {
            anglesy[ri][i]-=(int)mres/2;
            anglesx[ri][i]-=(int)mres/2;
        }
        
        
        maxr=ri;
        rcount++;
        std::cout<<"maxr="<<r<<std::endl;
        
    }
    
    //Number of points within circle
    int tic= ic[ri];
    
    //Number of points on edge
    int tic2=ic2[ri];
    
    
    bool xsi=false;
    bool xsi2=false;
    bool ysi=false;
    bool ysi2=false;
    
    
    
    //Get circle locations
    for (int l=1;l<ic[ri];l++)
    {
        ys[l]=anglesy[ri][l];
        xs[l]=anglesx[ri][l];
    }
    
    //Get circle edge locations
    for (int l=1;l<ic2[ri];l++)
    {
        ys2[l]=anglesy2[ri][l];
        xs2[l]=anglesx2[ri][l];
    }
    
    
    bool nidxy=false;
    bool nidxx=false;
    
    //Edge correct - don't use values which go out of bounds of the grid
    for (int i=1;i<tic2;i++)
    {
        if (ys2[i]+yi<1)
        {
            ys2[i]=-yi+1;
        }
        if (xs2[i]+xi<1)
        {
            xs2[i]=-xi+1;
        }
        if (ys2[i]+yi>ny)
        {
            ys2[i]=ny-yi;
        }
        if (xs2[i]+xi>nx)
        {
            xs2[i]=nx-xi;
        }
    }
    
    int c=1;
    for (int i=1;i<tic;i++)
    {
        ysi=false;
        xsi=false;
        ysi2=false;
        xsi2=false;
        nidxx=false;
        nidxy=false;
        nidx[i]=false;
        if (ys[i]+yi<1)
        {
            ysi=true;
        }
        if (xs[i]+xi<1)
        {
            xsi=true;
        }
        if(ys[i]+yi>ny)
        {
            ysi2=true;
        }
        if (xs[i]+xi>nx)
        {
            xsi2=true;
        }
        
        if (xsi||xsi2)
        {
            nidxx=true;
        }
        if (ysi||ysi2)
        {
            nidxy=true;
        }
        if (nidxy||nidxx)
        {
            nidx[i]=true;
        }
    }
    
    
    //Get linear indices on edge
    for (int i=1;i<tic2;i++)
    {
        anglesi2[i]=ys2[i]+yi+(xs2[i]+xi-1)*ny;
        
    }
    
    double mz=99999999;
    
    //Find lowest point on edge
    for (int i=1;i<tic2;i++)
    {
        if (Z[anglesi2[i]]<mz)
        {
            mz=Z[anglesi2[i]];
        }
        
    }
    double anglez2;
    
    //Calculate cone shape based on shelf and cone angle
    for (int i=1;i<ic[ri];i++)
    {
        
        anglesz[i]=Dist[ri][i]*std::tan(critangle*0.0174533)*dx+mz;
        anglez2=0-(r-Dist[ri][i])*std::tan(sangle*0.0174533)*dx;
        
        if (anglez2<anglesz[i])
        {
            anglesz[i]=anglez2;
            
        }
        
    }
    c=1;
    
    //Calculate linear indices of cone
    for(int i=1;i<tic;i++)
    {
        angleszi[i]=ys[i]+r+1+(2*r+1)*(xs[i]+r);
    }
    //Calculate linear indices of cone on grid.
    for(int i=1;i<tic;i++)
    {
        zni[i]=yi+ys[i]+ny*(xi+xs[i]-1);
    }
    
    return tic;
    
}
void ulemgrid::rangesearch(int r,int ri)
{
    double I1,I2;
    
    int idx1=1;
    int idx2=1;
    int ij;
    double x=(double)(mres/2.0);
    double y=(double)(mres/2.0);
    double r2=(double)r;
    
    for (int i=(int)(mres/2)-r-1;i<=(int)(mres/2)+r+1;i++)
    {
        for (int j=(int)(mres/2)-r-1;j<=(int)(mres/2)+r+1;j++)
        {
            double i2=(double)i;
            double j2=(double)j;
            if (std::pow(i2-y,2.0)+std::pow(j2-x,2.0)<=std::pow(r2,2.0))
            {
                anglesy[ri][idx1]=i;
                anglesx[ri][idx1]=j;
                idx1++;
                
                if (std::pow(i2-y,2.0)+std::pow(j2-x,2.0)>std::pow(r2-1.41,2.0))
                {
                    anglesx2[ri][idx2]=j;
                    anglesy2[ri][idx2]=i;
                    idx2++;
                }
            }
        }
    }
    ic[ri]=idx1;
    ic2[ri]=idx2;
    
}

void ulemgrid::diffuse()
{
    std::vector<double> diffuarray;
    diffuarray.resize(nn+1);
    diffuarray=Z;
    int ij;
    double diffy;
    double diffx;
    if (landorsea==0)
    {
        for (int i=1;i<ny*nx+1;i++)
        {
            if (diffuarray[i]<0)
            {
                diffuarray[i]=0;
            }
        }
    }
    
    
    
    else
    {
        for (int i=1;i<ny*nx+1;i++)
        {
            if (diffuarray[i]>0)
            {
                diffuarray[i]=0;
                
            }
        }
    }
    int mult =1;
    double dti=dt;
    while (dt>.2*dx*dx/kd)
    {
        dt = dt/2;
        mult = 2*mult;
        
    }
    
    std::vector<double> arrayii;
    for (int l=1;l<=mult;l++)
    {
        arrayii=diffuarray;
        for (int i=2;i<ny;i++)
        {
            
            for (int j=2;j<nx;j++)
            {
                
                ij=(j-1)*ny+i;
                //   if (sidx[ij])
                {
                    diffx=(arrayii[ij-ny]+arrayii[ij+ny]-2*arrayii[ij])/(dx*dx);
                    diffy=(arrayii[ij-1]+arrayii[ij+1]-2*arrayii[ij])/(dy*dy);
                    Z[ij]+=(diffx+diffy)*kd*dt;
                    diffuarray[ij]+=(diffx+diffy)*kd*dt;
                    
                }
                
            }
        }
        
        
    }
    dt=dti;
}






#endif


