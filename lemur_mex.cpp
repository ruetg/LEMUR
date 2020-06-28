

#include <vector>

#include<mex.h>
#include<matrix.h>
//#include <omp.h>
#include<ctime>
#include "lemur.cpp"
void mexFunction(
        int          nlhs,
        mxArray      *plhs[],
        int          nrhs,
        const mxArray *prhs[]
        )
{
    static bool first=true;
    int ny=1;
    int nx=1;
    if (first||nrhs==2)
    {
        if (nrhs==2&&first)
        {
            ny=(int)mxGetScalar(prhs[0]);
            nx=(int)mxGetScalar(prhs[1]);
            
            
        }
        else if (nrhs==0)
        {
            
            mexErrMsgTxt("LEM grid object not itinialized");
            
        }
        
        
        
    }
    
    static ulemgrid grid(ny,nx);
    
    if (first==false&&nrhs>=2)
    {
        std::string s;
        s=mxArrayToString(prhs[0]);
        
        if (s.compare("set")==0)
        {
            for (int i=1;i<nrhs;i+=2)
            {
                std::string nm=mxArrayToString(prhs[i]);
                int m,n;
                m=mxGetM(prhs[i+1]);
                n=mxGetN(prhs[i+1]);
                int siz=n*m;
                
                if (mxIsDouble(prhs[i+1])&&siz==1)
                {
                    double val=mxGetScalar(prhs[i+1]);
                    grid.set(nm,val);
                }
                else if (mxIsDouble(prhs[i+1])&&siz>1)
                {
                    double *val = mxGetPr(prhs[i+1]);
                    grid.set(nm,val,m*n);
                    
                }
                
                else if (~mxIsDouble(prhs[i+1])&&siz==1)
                {
                    
                }
            }
        }
        if (s.compare("get")==0)
        {
            for (int i=1;i<nrhs;i+=2)
            {
                std::string nm=mxArrayToString(prhs[i]);
                
                
                
                std::vector<double> val=grid.get(nm);
                
                plhs[0] = mxCreateNumericMatrix(1, val.size(),mxDOUBLE_CLASS, mxREAL);
                double* data = (double *) mxGetData(plhs[0]);
                for (int i=0;i<val.size();i++)
                {
                    data[i]=val[i];
                }
                
                
            }
        }
        else if(s.compare("run")==0)
        {
            std::string nm=mxArrayToString(prhs[1]);
            if (nm.compare("erode")==0)
            {
                grid.erosion_fluvial();
                
            }
            if (nm.compare("erode_fluvial")==0)
            {
                grid.erosion_fluvial();
            }
            if (nm.compare("erode_fluvial2")==0)
            {
                grid.erosion_fluvial2();
            }
            if (nm.compare("erode_diffusion")==0)
            {
                grid.diffuse();
            }
            if (nm.compare("lakefill")==0)
            {
                grid.lakefill();
            }
            if (nm.compare("deposition")==0)
            {
                grid.deposit();
            }
        }
    }
    
    first=false;
    
}