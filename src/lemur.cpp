
#include <string>
#include <fstream>
#include<cmath>
#include<sstream>
#include "lemur.h"

//#include <omp.h>



lemur::lemur(int i, int j)
{
    ny=i;
    nx=j;
    nn=i*j;
    
    ero.resize(nn+1);
    undercapacity.resize(nn+1);
    adds.resize(nn+1);
    Z.resize(nn+1);
        usedr.resize(nn+1);
    runoff.resize(nn+1);

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
    idx.at(0)=1;idx.at(1)=-1;idx.at(2)=ny;idx.at(3)=-ny;
    idx.at(4)=ny+1;
    idx.at(5)=ny-1;idx.at(6)=-ny+1;idx.at(7)=-ny-1;
    
    sinkfill.resize(nn+1);
    std::fill(sinkfill.begin(),sinkfill.end(),1);
    U.resize(nn+1);
    std::fill(U.begin(),U.end(),0);
    
    
    anglespolar1.resize(nn);
    anglespolar2.resize(nn);
    anglesz.resize(nn);
    
    anglesz2.resize(nn);
    
    
    watertot.resize(nn+1);

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


void lemur::set(std::string nm,double val)
{

    if (nm.compare("firstcall")==0)
    {
        firstcall=val;
    }
    else if (nm.compare("kd")==0)
    {
        kd=val;
    }
    else if (nm.compare("precip")==0)
    {
        precip = (double) val;
        
    }
    else if (nm.compare("evaprate")==0)
    {
        evaprate = (double) val;
        
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
        std::cout<<" 2 " << std::endl;
    }
    else if (nm.compare("uselandsed") == 0)
    {
        uselandsed=val;
    }
    else if (nm.compare("dy")==0)
    {
        dy=val;
        
                std::cout<<" 3 " << std::endl;

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
void lemur::set(std::string nm,double* val,int len)
{
    if (nm.compare("bc")==0)
    {
        BC.resize(len+1);
        for (int i=1;i<=len;i++)
        {
            BC.at(i)=val[i-1];
        }
        for (int i=0;i<len;i++)
        {
            BCX.at(val[i])=1;
        }
        
    }
    else if (len==nn)
    {
        
        if (nm.compare("z")==0)
        {
            for (int i=1;i<nn+1;i++)
            {
                Z.at(i)=val[i-1];
            }
            
        }
        else if (nm.compare("k")==0)
        {
            for (int i=1;i<nn+1;i++)
            {
                kval.at(i)=val[i-1];
            }
            
        }
        else if (nm.compare("undercapacity")==0)
        {
            for (int i=1;i<nn+1;i++)
            {
                undercapacity.at(i)=(int)val[i-1];
            }
            
        }
        
        else if (nm.compare("u")==0)
        {
            for (int i=1;i<nn+1;i++)
            {
                U.at(i)=val[i-1];
            }
            
        }
        else
        {
            std::cout<<nm<<": parameter not found";
        }
    }
    else
    {
        throw std::invalid_argument( " Error in lemur, input grid is wrong size" );
    }
    
    
}
void lemur::set(std::string nm,std::vector<double> val)
{
    
    int len=val.size();
    
    if (nm.compare("bc")==0)
    {
        std::fill(BCX.begin(),BCX.end(),0);
        
        for (int i=0;i<len;i++)
        {
            BCX.at(val[i])=1;
        }
        
    }
    else if (len==nn)
    {
        if (nm.compare("z")==0)
        {
            for (int i=1;i<nn+1;i++)
            {
                Z.at(i)=val[i-1];
            }
            
        }
        if (nm.compare("k")==0)
        {
            for (int i=1;i<nn+1;i++)
            {
                kval.at(i)=val[i-1];
            }
            //std::cout.at(kval.at(2))
            
        }
        if (nm.compare("undercapacity")==0)
        {
            
            for (int i=1;i<nn+1;i++)
            {
                undercapacity.at(i)=val[i-1];
            }
            
        }
        
        
    }
    else
    {
        throw std::invalid_argument( "Input grid is wrong size" );
    }
}



void lemur::findsteepest()
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
    double slpi;
    double diag = std::pow(std::pow((double)dx,2.0)+std::pow((double)dy,2.0),.5);
    double maxslp;
    int nsink=0;
    test.resize(nn+1);
#pragma omp simd
    for (int i=1;i<=nn;i++)
    {
        slpis.at(i)=i;
        dists.at(i)=100000000000000;
    }
#pragma omp parallel for private(slpyp,slpym,slpxm,slpxpyp,slpxmym,maxslp,ij) firstprivate(slpi,diag,ny,nx,dy,dx)
    for(int j=2;j<nx;j++)
        
    {
        for (int i=2;i <ny;i++)
        {
            
            ij = i+(j-1)*ny;
            if (BCX.at(ij)==1)
            {
                
                continue;
            }
            slpyp = (Z.at(ij)-Z.at(ij+1))/dy;
            slpym = (Z.at(ij)-Z.at(ij-1))/dy;
            slpxp = (Z.at(ij)-Z.at(ij+ny))/dx;
            slpxm = (Z.at(ij) -Z.at(ij-ny))/dx;
            slpxpyp = (Z.at(ij)-Z.at(ij+1+ny))/diag;
            slpxpym = (Z.at(ij)-Z.at(ij-1+ny))/diag;
            slpxmyp = (Z.at(ij)-Z.at(ij+1-ny))/diag;
            slpxmym = (Z.at(ij)-Z.at(ij-1-ny))/diag;
            maxslp = 0;
            slpi = 0;
            
            if (slpyp>maxslp)
            {
                maxslp = slpyp;
                slpi = 1;
                dists.at(ij) = dy;
            }
            
            if (slpym>maxslp)
            {
                maxslp = slpym;
                slpi = -1;
                dists.at(ij)=dy;
            }
            if(slpxp>maxslp)
                
            {maxslp = slpxp;
             slpi = ny;
             dists.at(ij)=dx;
            }
            if (slpxm>maxslp)
            {
                maxslp = slpxm;
                slpi = -ny;
                dists.at(ij)=dx;
            }
            if (slpxpyp>maxslp)
            {
                maxslp = slpxpyp;
                slpi = ny+1;
                dists.at(ij)=diag;
            }
            if (slpxpym>maxslp)
            {
                maxslp = slpxpym;
                slpi = ny-1;
                dists.at(ij)=diag;
            }
            if(slpxmyp>maxslp)
            {
                maxslp = slpxmyp;
                slpi = -ny+1;
                dists.at(ij)=diag;
            }
            if(slpxmym>maxslp)
            {
                maxslp = slpxmym;
                slpi = -ny-1;
                dists.at(ij)=diag;
                
            }
            if (maxslp==0)
            {
                nsink++;
            }
            
            slpis.at(ij) = slpi+ij;
            
        }
    }
    std::cout<<nsink<<std::endl;
    //mexPrintf("%d\n",nsink);
}

void lemur::getdonors()
{
#pragma omp parallel for
    for (int ij=1;ij<=nn;ij++)
    {
        if (slpis.at(ij)!=ij)
        {
            ndons.at(slpis.at(ij))++;
            donors.at(slpis.at(ij)).at(ndons.at(slpis.at(ij)))=ij;
        }
    }
}

void lemur::createstack()
{
    olist.resize(1);
    
#pragma omp simd
    for (int ij=1;ij<=nn;ij++)
    {
        if (slpis.at(ij)==ij)
        {
            olist.push_back(ij);
        }
    }
    stackij.resize(olist.size());
#pragma omp parallel for schedule(dynamic)
    for (int i=1;i<olist.size();i++)
    {
        int d;
        int p=1;
        int ij=olist.at(i);
        

        stackij.at(i).push_back(ij);
        p++;

        for (int ik=1;ik<=ndons.at(ij);ik++)
        {
            numiter=0;
            d = donors.at(ij).at(ik);
            p = recursivestack(d,p,i);
        }
    }
    olength.resize(olist.size());
#pragma omp parallel for
    for (int i=1;i<olist.size();i++)
    {
        olength.at(i)=stackij.at(i).size();
    }
    
    nnn=olist.size();
    int cc=1;
#pragma omp simd
    
    for (int j=1;j<nnn;j++)
    {
        
        for (int i=0; i<olength.at(j);i++)
        {
            stack.at(cc)=stackij.at(j).at(i);
            cc++;
        }
    }
    // mexPrintf("%d\n",stack.at(0));
}
int lemur::recursivestack(int r,int p,int i)
{
    
    if (p<=nn&&numiter<1e5)
    {
        int d;
        stackij.at(i).push_back(r);
        p++;
#pragma omp simd
        for (int ij=1;ij<=ndons.at(r);ij++)
        {
            d = donors.at(r).at(ij);
            p = recursivestack(d,p,i);
        }
        return p;
    }
    else
    {
        
        return 1;
    }
    
}
void lemur::createstack2()
{
    olist.resize(1);
    
#pragma omp simd
    for (int ij=1;ij<=nn;ij++)
    {
        if (slpis.at(ij)==ij)
        {
            olist.push_back(ij);
        }
    }
    
    stackij.resize(olist.size());
    
#pragma omp parallel for schedule(dynamic)
    for (int i=1;i<olist.size();i++)
    {
        int d;
        int ndon=0;
        int ij=olist.at(i);
        

        stackij.at(i).push_back(ij);
        
        for (int u=0;u<=ndon;u++)
        {
            
            ij=stackij.at(i).at(u);
            
            for (int ik=1;ik<=ndons.at(ij);ik++)
            {
                ndon++;
                d = donors.at(ij).at(ik);
                stackij.at(i).push_back(d);
                
            }
        }
    }
    olength.resize(olist.size());
#pragma omp parallel for
    for (int i=1;i<olist.size();i++)
    {
        olength.at(i)=stackij.at(i).size();
    }
    nnn=olist.size();
    int cc=1;
#pragma omp simd
    for (int j=1;j<nnn;j++)
    {
        
        for (int i=0; i<olength.at(j);i++)
        {
            stack.at(cc)=stackij.at(j).at(i);
            cc++;
        }
    }
    // mexPrintf("%d\n",stack.at(0));
}

void lemur::getacc()
{
    
#pragma omp simd
    for (int i=1;i<=nn;i++)
    {
        accgrid.at(i)=1.0;
    }
#pragma omp parallel for
    for (int j=1;j<nnn;j++)
    {
#pragma omp simd
        for (int i=olength.at(j)-1;i>0;i--)
        {
            
            accgrid.at(slpis.at(stackij.at(j).at(i)))=accgrid.at(slpis.at(stackij.at(j).at(i)))+accgrid.at(stackij.at(j).at(i));
        }
    }
    
}

void lemur::erode()
{
    double f;

    double x;
    int ni=1;
    if (n_>1)
    {
        ni=1;
    }
    if (n_<1)
    {
        ni=1;
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
                for (int i=1;i<olength.at(j);i++)
                {
                    if (undercapacity.at(stackij.at(j).at(i))==0)
                    {
                        if (Z.at(stackij.at(j).at(i))>=0)
                        {
                            
                            f = kval.at(stackij.at(j).at(i))/std::pow(dists.at(stackij.at(j).at(i)),n_) *
                                    std::pow(accgrid.at(stackij.at(j).at(i)),m)*sp1*
                                    std::pow(Zi.at(stackij.at(j).at(i))-std::max(0.0,Z.at(slpis.at(stackij.at(j).at(i)))),n_-1.0);
                            f2=0;//dt/(Zi.at(stackij.at(j).at(i))-std::max(0.0,Z.at(slpis.at(stackij.at(j).at(i))))); Gives weird result, probably wrong
                        }
                        x=1;
                        
                        for (int ll=1;ll<=5;ll++)
                        {
                            
                            
                            
                            x=(x-(x-1+f*std::pow(x,n_))/(1+(n_)*f*std::pow(x,n_-1.0)));
                        }
                        
                        Zi1=Z.at(stackij.at(j).at(i));
                        if (Z.at(stackij.at(j).at(i))>=0)
                        {///Z.at(stackij.at(j).at(i))=(f*Z.at(slpis.at(stackij.at(j).at(i)))+Z.at(stackij.at(j).at(i)))/(1+f);
                            Z.at(stackij.at(j).at(i)) = std::max(0.0,Z.at(slpis.at(stackij.at(j).at(i)))) +x*
                                    std::max(0.0,(Zi.at(stackij.at(j).at(i))-std::max(0.0,Z.at(slpis.at(stackij.at(j).at(i))))));
                            
                        }
                        // Z.at(stackij.at(j).at(i))=(Zi.at(stackij.at(j).at(i))+f*Z.at(slpis.at(stackij.at(j).at(i)))+dt*U.at(stackij.at(j).at(i)))/(1+f);}
                        if (kval.at(stackij.at(j).at(i))>0)
                        {
                            ero.at(stackij.at(j).at(i))=Zi1-Z.at(stackij.at(j).at(i));
                        }
                    }
                }
                
                
            }
        }
    }
}

void lemur::reset()

{

    if ( !usefd )
    {
        stackij.clear();
    }
    std::fill(accgrid.begin(),accgrid.end(),0);
   
    std::fill(sed.begin(),sed.end(),0);
    
    std::fill(ndons.begin(),ndons.end(),0);
    std::fill(BCX.begin(),BCX.end(),0);

    std::fill(ero.begin(),ero.end(),0.0);
#pragma omp parallel for
    for (int i=0;i<BC.size();i++)
    {
        
        BCX.at(BC.at(i))=1;
    }

  
    
    
    
    
}


void lemur::filll()
{
    if (stackij.size()>0)
    {
        
#pragma omp parallel for
        for (int j=1;j<nnn;j++)
        {
#pragma omp simd
            for (int i=0;i<olength.at(j);i++)
            {
                if (Z.at(stackij.at(j).at(i))<=Z.at(slpis.at(stackij.at(j).at(i))))
                {
                    Z.at(stackij.at(j).at(i))=Z.at(slpis.at(stackij.at(j).at(i)))+1e-6+rand()/RAND_MAX*1e-6;
                }
            }
        }
    }
}

void lemur::fillls()
{

    std::vector<double> cero=ero;
    std::vector<double> Zn=Z;
    double Lt;
    double ssum=0;double sedsum;
    double dti=dt;
    double A2;
    double di;
    double f;
    double d2;

    if (stackij.size()>0)
    {
        
#pragma omp parallel for schedule(dynamic) firstprivate(dti,ssum,sedsum,nj,dt,dx,dy,nn,ks,L,T)
        
        for (int j=1;j<nnn;j++)
        {
            
            double erol;
            

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
                        for (int i=1;i<olength.at(j);i++)
                        {
                            ero.at(stackij.at(j).at(i))/=2;
                        }
                        y=1;
                    }
                    y2++;

                    if (y2>5000)
                    {
                        break;
                    }
                    sedsum=0;
                    y++;
                    for (int i=olength.at(j)-1;i>=0;i--)
                    {
                        Zn.at(stackij.at(j).at(i))=Z.at(stackij.at(j).at(i));
                        
                        cero.at(stackij.at(j).at(i))=ero.at(stackij.at(j).at(i));
                        
                    }

                    Zn.at(slpis.at(stackij.at(j).at(0)))=Z.at(slpis.at(stackij.at(j).at(0)));
                    cero.at(slpis.at(stackij.at(j).at(0)))=ero.at(slpis.at(stackij.at(j).at(0)));
                    for (int i=olength.at(j)-1;i>0;i--)
                    {
                        if (stackij.at(j).at(i)!=slpis.at(stackij.at(j).at(i)))
                        {
                            cero.at(slpis.at(stackij.at(j).at(i)))+=cero.at(stackij.at(j).at(i));
                            
                        }
                    }

                    double tl;
                    if (m==1)
                    {
                        
                        f= dt * ks;
                    }
                    else
                    {
                        //std::cout<<kval.at(stackij.at(j).at(1));
                        if (stackij.at(j).size()>1){

                        f=dt*kval.at(stackij.at(j).at(1))*std::pow(dx*dy,m);
                        
                        }
                    }

                    double A;
                    
                    for (int i=1;i<=olength.at(j)-1;i++)
                    {                             std::cout<<i<<std::endl;

                        
                        if (undercapacity.at(stackij.at(j).at(i)) >= 1 && Zn.at(stackij.at(j).at(i)) >= 0 )
                        {
                            
                            A=std::pow(accgrid.at(stackij.at(j).at(i)),m);
                            
                            if (m==1)//If uncercapacity
                            {
                                Lt=L;
                                A2=1.0;
                                di=dists.at(stackij.at(j).at(i));
                                d2=1.0;
                            }
                            else//If yuan et al
                            {
                                Lt=1.0;
                                A2=accgrid.at(stackij.at(j).at(i));
                                di=1.0;
                                d2=std::pow(dists.at(stackij.at(j).at(i)),n_);

                            }
                            
           
                            tl=((cero.at(stackij.at(j).at(i))-ero.at(stackij.at(j).at(i)))/L/A2*di+Zn.at(stackij.at(j).at(i))+
                                    std::max(0.0,Zn.at(slpis.at(stackij.at(j).at(i))))*
                                    f/Lt/d2*A)/(f/Lt/d2*A+1.0); 
                            
          
                            
                            erol=ero.at(stackij.at(j).at(i));
                            ero.at(stackij.at(j).at(i))=(Zn.at(stackij.at(j).at(i))-tl);
                            
                            Zn.at(stackij.at(j).at(i))=tl;
                            
                            if (sedsum<std::abs(erol-ero.at(stackij.at(j).at(i))))
                            {sedsum=std::abs(erol-ero.at(stackij.at(j).at(i)));}
                            
                        }
                        
                        
                    }
                }
                
                for (int i=0;i<=olength.at(j)-1;i++)
                {
                    if (undercapacity.at(stackij.at(j).at(i))>=1)
                    {Z.at(stackij.at(j).at(i))=Zn.at(stackij.at(j).at(i));
                     
                    }
                }
            }
            
            
            
            {ssum+=cero.at(stackij.at(j).at(0));
            }
            
        }


    }
    dt=dti;
    
    
    
}

void lemur::fillls2()
{
    /**
     *
     */
    std::vector<double> cero=ero;
    std::vector<double> Zn=Z;
    
    double ssum=0;double sedsum;
    double dti=dt;
    
    if (stackij.size()>0)
    {
        int dts=200;
        dt/=200;
        int no=1;
        std::vector<double> cerol;
        cerol.resize(nn+1);
#pragma omp parallel for schedule(dynamic) firstprivate(dti,ssum,sedsum,nj,dt,dx,dy,nn,ks,L,T)
        
        for (int j=1;j<nnn;j++)
        {
            double erol;
            
            int y;
            
            int nt1=1;
            
            
            while (nt1<dts)
            {
                y=1;
                nt1++;
                
                sedsum=1.0001e-3;
                
                
                {
                    
                    

                    sedsum=0;
                    
                    y++;
                    
                    
                    Zn.at(slpis.at(stackij.at(j).at(0)))=Z.at(slpis.at(stackij.at(j).at(0)));
                    
                    
                    double tl;
                    double f= dt*ks;
                    double A;
                    for (int i=olength.at(j)-1;i>0;i--)
                    {
                        cerol.at(stackij.at(j).at(i))=cero.at(stackij.at(j).at(i));
                        cerol.at(slpis.at(stackij.at(j).at(i)))=cero.at(slpis.at(stackij.at(j).at(i)));
                        
                        cero.at(stackij.at(j).at(i))=0;
                        cero.at(slpis.at(stackij.at(j).at(i)))=0;
                        
                    }
                    cero.at(stackij.at(j).at(olength.at(j)-1))=0;
                    cero.at(stackij.at(j).at(0))=0;
                    for (int i=olength.at(j)-1;i>0;i--)
                    {
                        if (nt1==2)
                            
                        {                        no++;
                        }
                        if (undercapacity.at(stackij.at(j).at(i))>=1)
                        {
                            
                            A=std::pow(accgrid.at(stackij.at(j).at(i)),1.0);
                            
                            tl=1/L*dists.at(stackij.at(j).at(i))*(cero.at(stackij.at(j).at(i))-
                                    ks*dt*accgrid.at(stackij.at(j).at(i))*
                                    (Z.at(stackij.at(j).at(i))-std::max(0.0,Z.at(slpis.at(stackij.at(j).at(i)))))/dists.at(stackij.at(j).at(i)));
                            
                            erol=ero.at(stackij.at(j).at(i));
                            
                            Zn.at(stackij.at(j).at(i))+=tl;
                            cero.at(stackij.at(j).at(i))-=tl;
                            cero.at(slpis.at(stackij.at(j).at(i)))+=cero.at(stackij.at(j).at(i));
                            
                            
                            if (sedsum<std::abs(erol-ero.at(stackij.at(j).at(i))))
                            {sedsum=std::abs(erol-ero.at(stackij.at(j).at(i)));}
                        }
                        
                        
                    }
                }
                for (int i=1;i<=olength.at(j)-1;i++)
                {
                    if (undercapacity.at(stackij.at(j).at(i))>=1)
                    {Z.at(stackij.at(j).at(i))=Zn.at(stackij.at(j).at(i));
                     
                    }
                }
            }
            
            
            
            {ssum+=cero.at(stackij.at(j).at(0));
            }
        }
        
        
    }
    dt=dti;
    
    
}

void lemur::basinfill(int ij,int ij1)
{
    double min=Z.at(ij);
    for (int i=0;i<8;i++)
    {
        if (Z.at(ij+idx.at(i))<min)
        {
            min=Z.at(ij+idx.at(i));
        }
        
    }
    if (Z.at(ij)<=min)
    {
        Z.at(ij)=min+.01;
        adds.at(ij1)+=min+.01-Z.at(ij);
        
        basinfill(ij,ij1);
        for (int i=0;i<8;i++)
        {
            if (BCX.at(ij+idx.at(i))==0)
            {basinfill(ij+idx.at(i),ij1);}
        }
    }
}
void lemur::erosion_fluvial()
{    
if (firstcall == 0)
{reset();
}
        
    for (int i=1;i<=nn;i++)
    {
        Z.at(i) += watertot.at(i);
    }

    if (!usefd)
    {
        findsteepest();

        getdonors();

        createstack2();

    }

    getacc();
        runoff=accgrid;

    for (int i=1;i<=nn;i++)
    {
        Z.at(i)-=watertot.at(i);
    }

    erode();

    if (ks>1e-100)
    {
        fillls();
    }
    for (int i=1;i<=ny;i++)
    {
        
        if (BCX.at(i)==0)
        {
            Z.at(i)=Z.at(i+ny);
            ero.at(i)=0;
        }
    }
#pragma omp parallel for
    for (int i=nn-ny+1;i<=nn;i++)
    {
        if (BCX.at(i)==0)
        {
            Z.at(i)=Z.at(i-ny);
            ero.at(i)=0;
        }
        
    }
#pragma omp parallel for
    for (int i=1;i<=nn;i+=ny)
    {
        if (BCX.at(i)==0)
        {
            Z.at(i)=Z.at(i+1);
            ero.at(i)=0;
        }
        
    }
    for (int i=ny;i<=nn;i+=ny)
    {
        if (BCX.at(i)==0)
        {
            Z.at(i)=Z.at(i-1);
            ero.at(i)=0;
        }
    }
    firstcall = 0;

}

void lemur::erosion_fluvial2()
{
    reset();
    if (!usefd)
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

}

std::vector<double> lemur::get(std::string nm)
{
    std::vector<double> val;
    val.resize(nn);
    if (nm.compare("bcx")==0)
    {
        for (int i=1;i<nn;i++)
        {
            val[i-1]=BCX.at(i);
        }
    }
    else if (nm.compare("acc")==0)
    {
        for (int i=1;i<nn;i++)
        {
            val[i-1]=accgrid.at(i);
        }
    }
    else if (nm.compare("bc")==0)
    {
        for (int i=0;i<BC.size();i++)
        {
            val[i-1]=accgrid.at(i);
        }
    }
    else if (nm.compare("landsurf")==0)
    {
        for (int i=1;i<nn;i++)
        {
            val[i-1]=landsurf.at(i);
        }
    }
    else if (nm.compare("z")==0)
    {
        std::cout<<"z";
 
        for (int i=1;i<nn+1;i++)
        {
            val[i-1]=Z.at(i);
        }
        
    }
    else if (nm.compare("sinkareas")==0)
    {
        for (int i=1;i<nn+1;i++)
        {
            val[i-1]=sinkareas.at(i);
        }
        
    }

    else if (nm.compare("k")==0)
    {
//std::cout<<kval.at(2);
        for (int i=1;i<nn+1;i++)
        {
            val[i-1]=kval.at(i);
        }
        
    }
    else if (nm.compare("watertot")==0)
    {
        for (int i=1;i<nn+1;i++)
        {
            val[i-1]=watertot.at(i);
        }
    }
        else if (nm.compare("k")==0)
    {
        for (int i=1;i<nn+1;i++)
        {
            val[i-1]=kval.at(i);
        }
        
    }
    else if (nm.compare("undercapacity")==0)
    {
        for (int i=1;i<nn+1;i++)
        {
            val[i-1]=undercapacity.at(i);
        }
        
    }

    else if (nm.compare("stack")==0)
    {
         for (int i=1;i<nn+1;i++)
        {
            val[i-1]=stack.at(i);
        }
    }
    else if (nm.compare("rec")==0)
    {
         for (int i=1;i<nn+1;i++)
        {
            val[i-1]=slpis.at(i);
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
    priorityq(std::vector<double>&);

    void push(int);
    int pop();
    int top();
    
private:
    int nn;
    int numel=0;

    int queue;
    int tmp;
    int u;
    int uu;
    int ul;
    std::vector<int> left;
    std::vector<double> z;
    
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
    return left.at(1);
    
}
int priorityq::pop()
{
    
    
    uu =left.at(1);
    
    left.at(1)=left.at(numel);
    left.at(numel)=0;
    u=2;
    ul=u/2;
    while (u<numel-2)
    {
        if (z.at(left.at(u))<z.at(left.at(u+1)))
        {
            
            if (z.at(left.at(ul))>z.at(left.at(u)))
            {
                
                std::swap(left.at(ul),left.at(u));
                
                ul=u;
                u=2*u;
                
                
            }
            else
            {
                break;
            }
        }
        else if(z.at(left.at(ul))>z.at(left.at(u+1)))
        {
            std::swap(left.at(ul),left.at(u+1));
            
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
    
    left.at(u)=i;
    while (ul>0)
    {
        
        if (z.at(left.at(ul))>z.at(left.at(u)))
        {
            
            std::swap(left.at(ul),left.at(u));
            
            
            
        }
        else
        {
            break;
        }
        u=u/2;
        ul=u/2;
        
        
        
        
    }
    
}

void lemur::lakefill2()
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
    if (uselandsed==0)
    {
        
        for (int i=0;i<BC.size();i++)
        {
            open.push(BC.at(i));
            closed.at(BC.at(i))=true;
            c++;
        }
    }
    
    for (int i=1;i<=ny;i++)
    {
        if (!closed.at(i))
        {
            closed.at(i)=true;
            
            open.push(i);
            
            c++;
        }
        
    }
    for (int i=nn-ny+1;i<=nn;i++)
    {
        if (!closed.at(i))
        {
            closed.at(i)=true;
            
            open.push(i);
            
            c++;
        }
    }
    for (int i=ny;i<=nn;i+=ny)
    {
        if (!closed.at(i))
        {
            closed.at(i)=true;
            
            open.push(i);
            
            c++;
        }
    }
    for (int i=ny+1;i<=nn;i+=ny)
    {
        if (!closed.at(i))
        {
            closed.at(i)=true;
            
            open.push(i);
            
            c++;
        }
    }
    int s;
    int ij;
    int ii;
    int jj;
    pittop=-9999;
    while (c>0||p>0)
    {
        if (p>0&&c>0&&pit.at(p-1)==-9999)//pq.at(c-1)==pit.at(p-1))
        {
            //s=open.top();
            s=open.pop();
            c--;
            pittop=-9999;
        }
        else if (p>0)
        {
            
            s=pit.at(p-1);
            pit.at(p-1)=-9999;
            p--;
            if (pittop==-9999)
            {
                pittop=Z.at(s);
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
            
            ij = idx.at(i)+s;
            
            ii= (ij-1)%ny+1;
            jj = (int)((ij-1)/ny)+1;
            if ((ii>=1)&&(jj>=1)&&(ii<=ny)&&(jj<=nx) && !closed.at(ij))
            {
                
                closed.at(ij)=true;
                
                if (Z.at(ij)<=Z.at(s))
                {
                    
                    Z.at(ij)=Z.at(s)+1e-6 + rand()/RAND_MAX*1e-6;
                    
                    pit.at(p)=ij;
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


void lemur::landsed()
{

    sinkareas.resize(nn+1);
    std::fill(sinkareas.begin(),sinkareas.end(),0.0);
    std::fill(watertot.begin(),watertot.end(),0.0);

    sed = ero;
    double addprecip;
    double fact2;
    double fact;

    if (sed.size()<nn){sed.resize(nn+1);}
    std::vector<double> sedsub;
    sedsub.resize(nn+1);
    next.resize(nn+1);
    std::fill(sedsub.begin(),sedsub.end(),0);
    sedlist.resize(nn+1);
    std::fill(usedr.begin(),usedr.end(),false);

    double tsed=0;
    for (int i=1;i<=ny;i++)
    {
        usedr.at(i)=true;
    }
    
    for (int i=nn-ny+1;i<=nn;i++)
    {
        usedr.at(i)=true;
    }
    
    for (int i=ny;i<=nn;i+=ny)
    {
        usedr.at(i)=true;
    }
    for (int i=ny+1;i<=nn;i+=ny)
    {
        usedr.at(i)=true;
    }
    std::fill(runoff.begin(),runoff.end(),precip);
    //build cumulative sediment from ero and stream network
    //zero sediment everywhere but outlet to not count it twice
    //#pragma omp parallel for

    for (int i =nn;i>=1;i--)
    {
        
        tsed+=ero.at(stack.at(i));
        if (slpis.at(stack.at(i))!=stack.at(i))
        {
            sed.at(slpis.at(stack.at(i))) +=sed.at(stack.at(i));
            runoff.at(slpis.at(stack.at(i))) += runoff.at(stack.at(i));
        }
    }
    int nsink=0;
    //#pragma omp parallel for
    for (int i =1;i<=nn;i++)
    {
        if (slpis.at(stack.at(i))!=stack.at(i))
        {
            
            sed.at(stack.at(i))=0;
            runoff.at(stack.at(i))=0;
            
        }
        else
        {
            nsink++;
        }
    }
    std::cout<<"here1x"<<std::endl;

    for (int i = 1;i <= nn;i++)
    {
        if( (landsurf.at(i) < Z.at(i)) && !(usedr.at(i)) )
        {

            nseds=0;
            massextra=0;
            massextra_precip=0;
            area=0;
            volume=0;
            nrecur=0;
            //std::fill(sedlist.begin(),sedlist.end(),0);
            recursivesed(i);
            fact=0;
            fact2=0;
            if (massextra>=volume)
            {
                fact=1;
            }
            else if (volume>0)
            {
                fact = massextra/volume;
            }
            
            if (massextra_precip*dt+massextra-evaprate*area*dt>=volume)
            {
                fact2=1;
            }
            else if (volume-massextra>0)
            {
                fact2 = std::max(0.0,(massextra_precip*dt-evaprate*area*dt)/(volume-massextra));
            }
            

            double addsed;
            for (int j=1;j<=nseds;j++)
            {
                if (sedlist.at(j)>nn||sedlist.at(j)<0)
                {std::cout<<"warning: bad access in recursive sedfill"<<std::endl;
                 break;
                }
                
                addsed = (Z.at(sedlist.at(j))-landsurf.at(sedlist.at(j)))*fact;
                landsurf.at(sedlist.at(j)) += addsed;
                addprecip = (Z.at(sedlist.at(j))-(landsurf.at(sedlist.at(j)))) * fact2;
                
                //landsurf.at(sedlist.at(j)) += addprecip;
                watertot.at(sedlist.at(j)) += addprecip;
                sinkareas.at(sedlist.at(j)) = area;
                
            }
            
        }

    }
    std::cout<<"here2x"<< std::endl;

    if (maxareasinkfill>=0)
    {
        for (int i=1;i<nn+1;i++)
        {
            if (sinkareas.at(i)<maxareasinkfill)
            {
                landsurf.at(i) = Z.at(i);
            }
        }
    }
    
    Z=landsurf;
}
void lemur::recursivesed(int ij)
{
    usedr.at(ij)=true;
    int c=0;
    bool go = true;
    while (go)
    {
        
        nrecur++;
        nseds++;
        
        sedlist.at(nseds)=ij;
        for (int i=0;i<8;i++)
        {
            
            if( ( landsurf.at(idx.at(i)+ij)<Z.at(idx.at(i)+ij) ) && (!usedr.at(idx.at(i)+ij)) )
            {
                next.at(c)=idx.at(i)+ij;
                usedr.at(ij+idx.at(i))=true;
                c++;
            }
            
        }
        
        volume+=(Z.at(ij)-landsurf.at(ij));
        area+=1;
        if (uselandsed==1)
        {
            massextra+=sed.at(ij);
        }
        else if (uselandsed ==2)
        {
            massextra_precip+=runoff.at(ij)*precip;
            
        }
        else if (uselandsed == 3)
        {
            massextra += sed.at(ij);
            massextra_precip += runoff.at(ij)*precip;
        }
        c--;
        if (c>=0)
        {
            ij=next.at(c);
            
        }
        else
        {
            go=false;
        }
        
    }
}
void lemur::checkparam(std::string name, std::vector<double> &var)
{
    if (var.size() < nn + 1 )
    {
        std::cout<< "\n" << name << " is not the correct sizef "<< "\n";

        return;
    }
    else
    {
        for (int i = 0 ;i <= nn; i ++)
        {
            if ((var.at(i) <-1e15) || (var.at(i) > 1e15))
            {
                std::cout<< "\n" << name << " not set properly "<< "\n";
                return;
            }
        }
    }
     
}

void lemur::checkparam(std::string name, std::vector<int> &var)
{
    if (var.size() < nn + 1 && !name.compare("bc"))
    {
        std::cout<< "\n" << name << " is not the correct size "<< "\n";
        return;
    }
    else if (var.size() < nn + 1 && name.compare("bc"))
    {
        for (int i = 0 ;i < BC.size(); i ++)
        {
            if (var.at(i) < 0) 
            {
                std::cout<< "\n" << name << " not set properly -"<< "\n";
            }
            if ((var.at(i) > nn+1))
            {
                std::cout<< "\n" << name << " not set properly +"<< "\n";

            }
        }
    }
    else
    {
        for (int i = 0 ;i <= nn; i ++)
        {
            if ((var.at(i) < 0) || (var.at(i) > 10))
            {
                std::cout<< "\n" << name << " not set properly "<< "\n";

                return;

            }
        }
    }
     
}

void lemur::checkparams()
{
    std::string p = "Z";
    checkparam( p, Z );
    p = "BC";
    checkparam( p, BC );

}



void lemur::lakefill()
{
    std::cout<<"111"<<std::endl;
    checkparams();
    if (uselandsed>0)
    {
       // if (firstcall==1)
      //  {landsurf=Z;}
        landsurf=Z;
        if (firstcall ==0)
        {
         reset();
        }
         findsteepest();
            

        getdonors();
            
        createstack2();
        
        lakefill2();
        
        landsed();
    }
    else
    {
        uselandsed=0;
        landsurf = Z;
        lakefill2();
        landsurf = Z;
    }
}
void lemur::deposit()
{
    double sumsedij=0;
    std::vector<double> dsi;
    int ij;
    double sumer=0;
    //compute sediment using stack and erosion
    sed=ero;
    for (int i = nn;i>=1;i--)
    {
        if (slpis.at(stack.at(i))!=stack.at(i))
        {
            
            
            sumer+=ero.at(stack.at(i));
            sed.at(slpis.at(stack.at(i)))+=sed.at(stack.at(i));
        }
        
    }
    dsi.resize(nn+1);
    int cat=2;
    std::vector<bool> ends;
    ends.resize(nn+1);
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
        
        
        catchments.at(stack.at(i))=catchments.at(slpis.at(stack.at(i)));
        if (catchments.at(stack.at(i))==1)
        {
            
            catchments.at(stack.at(i))=cat;
            ends.at(stack.at(i))=true;
            dsi.at(slpis.at(stack.at(i)))=sed.at(slpis.at(stack.at(i)));
            
            cat++;
        }
    }
    int rs;
    std::cout<<"nx= "<< nx << std::endl;
    std::cout<<"ny= "<< ny << std::endl;
    for (int i=1;i<=ny;i++)
    {

        for (int j=1;j<=nx;j++)
        {
            ij = i+(j-1)*ny;
            if (dsi.at(ij)>0&&(((i>1&&i<ny&&j>1&&j<nx))&&Z.at(ij)<0))
            {
                rs=minrs;
                test.at(ij)=1;
                double extrased=1;
                int rsii=rsi;
                ri=0;
                while (extrased>0)
                {
                    
                    tic=itriangle(rs,ri,j,i);
                    ri++;
                    
                    for (int l=1;l<tic;l++)
                    {
                        
                        if (!nidx.at(l))
                        {
                            
                            if (anglesz.at(l)>0)
                            {
                                anglesz.at(l)=0;
                            }
                        }
                    }
                    
                    double sumseds=0;
                    double circsurf;
                    
                    
                    
                    for (int l=1;l<tic;l++)
                    {
                        if (!nidx.at(l))
                        {
                            circsurf=anglesz.at(l);
                            seds.at(l)=circsurf-Z.at(zni.at(l));
                            seds.at(l)=(seds.at(l)+std::fabs(seds.at(l)))/2.0;
                            sumseds+=seds.at(l);
                            
                        }
                    }
                    
                    extrased=sed.at(ij)-sumseds;
                    extrased=(extrased+std::fabs(extrased))/2.0;
                    
                    
                    
                    if (extrased==0)
                    {
                        
                        double sedx;
                        for (int l=1;l<tic;l++)
                        {
                            if (!nidx.at(l))
                            {
                                sedx=sed.at(ij)/sumseds*seds.at(l);
                                Z.at(zni.at(l))=Z.at(zni.at(l))+sedx;
                                sumsedij+=sed.at(ij)/sumseds*sedx;
                                
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
                            
                            if (!nidx.at(l))
                            {
                                
                                
                                sedx=seds.at(l);
                                
                                if (sumseds<=sed.at(ij))
                                {
                                    Z.at(zni.at(l))=Z.at(zni.at(l))+sedx;
                                }
                                else
                                {
                                    Z.at(zni.at(l))=Z.at(zni.at(l))+sed.at(ij)/sumseds*sedx;
                                }
                                
                            }
                        }
                        
                        break;
                    }

                }
                
                
                
                
            }
            
        }
    }
    
}




int lemur::itriangle(int r,int ri, int xi, int yi)
{
    double I1;
    
    //if geometry of given radius has not previously been solved for...
    
    if (ri>maxr)
    {
        
        //x,y locations of area within circle
        anglesx.at(ri).resize(4*r*r+4);
        anglesy.at(ri).resize(4*r*r+4);
        
        //x,y locations on edge of circle
        anglesx2.at(ri).resize(4*r*r+4);
        anglesy2.at(ri).resize(4*r*r+4);
        
        
        I.at(ri).resize(4*r*r+4);
        Dist.at(ri).resize(4*r*r+4);
        Dist2.at(ri).resize(4*r*r+4);
        
        //Find cells within circle and on edge
        rangesearch(r,ri);
        
        
        //Find the theta value of cells within circle
        for (int i=1;i<ic.at(ri);i++)
        {
            
            anglespolar1.at(i) = std::atan2((double)anglesy.at(ri).at(i)-(double)mres/2,(double)anglesx.at(ri).at(i)-(double)mres/2);
            
        }
        
        //Find the theta value of cells on edge
        for (int i=1;i<ic2.at(ri);i++)
        {
            anglespolar2.at(i) = std::atan2((double)anglesy2.at(ri).at(i)-(double)mres/2,(double)anglesx2.at(ri).at(i)-(double)mres/2);
            
        }
        
        //Find the closest point on edge of circle, to that within the circle
        //based on theta value
        for (int i=1;i<ic.at(ri);i++)
        {
            double mI1=9999999;
            for (int j=1;j<ic2.at(ri);j++)
            {
                
                I1=std::fabs((double) anglespolar1.at(i)-(double)anglespolar2.at(j));
                
                if (I1<mI1)
                {
                    mI1=I1;
                    I.at(ri).at(i)=j;
                }
                
            }
        }
        
        //calculate the distance from each interior point to closest edge point
        for (int i=1;i<ic.at(ri);i++)
        {
            Dist.at(ri).at(i)=std::sqrt(std::pow((double)anglesy2.at(ri).at(I.at(ri).at(i))-
                    (double)anglesy.at(ri).at(i),2.0)+
                    std::pow((double)anglesx2.at(ri).at(I.at(ri).at(i))-
                    (double)anglesx.at(ri).at(i),2.0));
        }
        
        //make dimension-independent
        for (int i=1;i<ic2.at(ri);i++)
        {
            anglesy2.at(ri).at(i)-=(int)mres/2;
            anglesx2.at(ri).at(i)-=(int)mres/2;
        }
        //Not sure .
        for (int i=1;i<ic.at(ri);i++)
        {
            Dist2.at(ri).at(i)=std::sqrt(std::pow(anglesy2.at(ri).at(I.at(ri).at(i)),2.0) +std::pow(anglesx2.at(ri).at(I.at(ri).at(i)),2.0));
            
        }
        
        //Make dimension independent
        for (int i=1;i<ic.at(ri);i++)
        {
            anglesy.at(ri).at(i)-=(int)mres/2;
            anglesx.at(ri).at(i)-=(int)mres/2;
        }
        
        
        maxr=ri;
        rcount++;
        std::cout<<"maxr="<<r<<std::endl;
        
    }


    //Number of points within circle
    int tic= ic.at(ri);
    
    //Number of points on edge
    int tic2=ic2.at(ri);
    
    
    bool xsi;
    bool xsi2;
    bool ysi;
    bool ysi2;
    
    
    
    //Get circle locations
    for (int l=1;l<ic.at(ri);l++)
    {
        ys.at(l)=anglesy.at(ri).at(l);
        xs.at(l)=anglesx.at(ri).at(l);
    }
    
    //Get circle edge locations
    for (int l=1;l<ic2.at(ri);l++)
    {
        ys2.at(l)=anglesy2.at(ri).at(l);
        xs2.at(l)=anglesx2.at(ri).at(l);
    }
    
    
    bool nidxy;
    bool nidxx;
    
    //Edge correct - don't use values which go out of bounds of the grid
    for (int i=1;i<tic2;i++)
    {
        if (ys2.at(i)+yi<1)
        {
            ys2.at(i)=-yi+1;
        }
        if (xs2.at(i)+xi<1)
        {
            xs2.at(i)=-xi+1;
        }
        if (ys2.at(i)+yi>ny)
        {
            ys2.at(i)=ny-yi;
        }
        if (xs2.at(i)+xi>nx)
        {
            xs2.at(i)=nx-xi;
        }
    }
    

    for (int i=1;i<tic;i++)
    {
        ysi=false;
        xsi=false;
        ysi2=false;
        xsi2=false;
        nidxx=false;
        nidxy=false;
        nidx.at(i)=false;
        if (ys.at(i)+yi<1)
        {
            ysi=true;
        }
        if (xs.at(i)+xi<1)
        {
            xsi=true;
        }
        if(ys.at(i)+yi>ny)
        {
            ysi2=true;
        }
        if (xs.at(i)+xi>nx)
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
            nidx.at(i)=true;
        }
    }
    
    
    //Get linear indices on edge
    for (int i=1;i<tic2;i++)
    {
        anglesi2.at(i)=ys2.at(i)+yi+(xs2.at(i)+xi-1)*ny;
        
    }
    
    double mz=99999999;
    
    //Find lowest point on edge
    for (int i=1;i<tic2;i++)
    {
        if (Z.at(anglesi2.at(i))<mz)
        {
            mz=Z.at(anglesi2.at(i));
        }
        
    }
    double anglez2;
    
    //Calculate cone shape based on shelf and cone angle
    for (int i=1;i<ic.at(ri);i++)
    {
        
        anglesz.at(i)=Dist.at(ri).at(i)*std::tan(critangle*0.0174533)*dx+mz;
        anglez2=0-(r-Dist.at(ri).at(i))*std::tan(sangle*0.0174533)*dx;
        
        if (anglez2<anglesz.at(i))
        {
            anglesz.at(i)=anglez2;
            
        }
        
    }
    
    //Calculate linear indices of cone
    for(int i=1;i<tic;i++)
    {
        angleszi.at(i)=ys.at(i)+r+1+(2*r+1)*(xs.at(i)+r);
    }
    //Calculate linear indices of cone on grid.
    for(int i=1;i<tic;i++)
    {
        zni.at(i)=yi+ys.at(i)+ny*(xi+xs.at(i)-1);
    }
    
    return tic;
    
}
void lemur::rangesearch(int r,int ri)
{

    int idx1=1;
    int idx2=1;
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
                anglesy.at(ri).at(idx1)=i;
                anglesx.at(ri).at(idx1)=j;
                idx1++;
                
                if (std::pow(i2-y,2.0)+std::pow(j2-x,2.0)>std::pow(r2-1.41,2.0))
                {
                    anglesx2.at(ri).at(idx2)=j;
                    anglesy2.at(ri).at(idx2)=i;
                    idx2++;
                }
            }
        }
    }
    ic.at(ri)=idx1;
    ic2.at(ri)=idx2;
    
}

void lemur::diffuse()
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
            if (diffuarray.at(i)<0)
            {
                diffuarray.at(i)=0;
            }
        }
    }
    
    
    
    else
    {
        for (int i=1;i<ny*nx+1;i++)
        {
            if (diffuarray.at(i)>0)
            {
                diffuarray.at(i)=0;
                
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
                //   if (sidx.at(ij))
                {
                    diffx=(arrayii.at(ij-ny)+arrayii.at(ij+ny)-2*arrayii.at(ij))/(dx*dx);
                    diffy=(arrayii.at(ij-1)+arrayii.at(ij+1)-2*arrayii.at(ij))/(dy*dy);
                    Z.at(ij)+=(diffx+diffy)*kd*dt;
                    diffuarray.at(ij)+=(diffx+diffy)*kd*dt;
                    
                }
                
            }
        }
        
        
    }
    dt=dti;
}



