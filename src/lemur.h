
#ifndef LEMUR_H
#define LEMUR_H

#include <iostream>
#include <vector>
//#include "lemur.cpp"
#include <fstream>


class lemur
{
public:
    lemur(int,int);//constructor takes size of grid mxn as input params
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
    int firstcall = 1;
    std::vector<int> slpis ;
    double kd;
    double tt;
    double n_;
    double maxareasinkfill;
    std::vector<double> BC;
    std::vector<std::vector<int> > stackij;
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
    double volume;
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
    int uselandsed;
    std::vector<bool> usedr;
    void recursivesed(int);
    double massextra;
    double area;
    std::vector<int>sedlist;
    std::vector<int> next;
    double precip = 1;
    int nrecur;
    int nseds;
    double evaprate =0;
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
    std::vector<std::vector<int> > anglesx;
    std::vector<std::vector<int> > anglesy;
    std::vector<std::vector<int> > anglesx2;
    std::vector<std::vector<int> > anglesy2;
    std::vector<std::vector<double> > Dist;
    std::vector<std::vector<int> > I;
    std::vector<int> ic;
    std::vector<int> ic2;
    double massextra_precip = 0;

    std::vector<double> watertot;
    std::vector<std::vector<double> > Dist2;
    
    std::vector<double> anglespolar1;
    std::vector<double> anglespolar2;
    std::vector<int> ys;
    std::vector<int> xs;
    std::vector<int> ys2;
    std::vector<int> xs2;
    std::vector<double> runoff;
    
    std::vector<int> anglesi2;
    std::vector<int> zni;
    std::vector<double>anglesz;
    std::vector<double>angleszi;
    std::vector<double> minz,anglesz2;
    std::vector<bool> nidx;
    int landorsea=0;//0 for land diffuse, 1 for sea
    
};
#endif



// Change number numbering system of rivers
// Look over esurf
// Quantitative metrics when describing correlations/durations - tea leaves, discussion
// Cumulative effect on eroded sediment
// Other rivers - not s1 - s12, unaffected by landslides (control rivers)
// Same 12 rivers, show it is a localized effect - spatially (comparison to rivers hit hard by ls vs 
// those that weren't 
// flood distribution 
// comparison of different methods
// how long do effects last
// a and b before morakot - stable mean before / after (?) 
