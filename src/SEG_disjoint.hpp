#ifndef SEG_DISJOINT_HPP_INCLUDED
#define SEG_DISJOINT_HPP_INCLUDED

#include "_help.hpp"

#define SQR(x) (x)*(x)
#define X3(ix,iy,iz) ((iz)*iNy+(iy))*iNx+(ix)

namespace image
{
namespace segmentation
{

    
// The code in this function is inspired by Matlab code by Tom Goldstein & Xavier Bresson [1]
// [1] Tom Goldstein et al. "Geometric applications of the split Bregman method: segmentation and surface reconstruction." J. of Sc. Comp. 45.1-3 (2010): 272-293.
template<typename pixel_type,unsigned int dimension>
void SEG_disjoint(std::vector<image::basic_image<pixel_type,dimension> >& totalIn,
            std::vector<image::basic_image<pixel_type,dimension> >& totalOut,
            image::basic_image<pixel_type,dimension>& alpha,
         
            image::basic_image<pixel_type,dimension>& brainMask,
            std::vector<image::basic_image<pixel_type,dimension> >& u,
        
            pixel_type pMin,
        
            int iMaxIter,
            pixel_type pEps,
            
            pixel_type pLambda,
            pixel_type pMu,
            
            bool bUpdateAlpha,
            
            pixel_type pCoupling,
        
            std::string basename,
            bool verbose,
            
            int iInclusion_toWhich,
            image::basic_image<pixel_type,dimension>& inclusion){

    
    // DECLARE INT
    int iNrOfChannels = u.size();
    
    geometry<dimension> geo = totalIn[0].geometry();
    
    int iNx = geo[0];
    int iNy = geo[1];
    int iNz = geo[2];
    
    std::vector<std::vector<pixel_type> > lesionIndicator;
    int iNrOfLesionStates = image::fillLesionIndicator(lesionIndicator, iNrOfChannels);
    
    int ix;
    int iy;
    int iz;
    int iX;
    std::vector<int> iGS(iNrOfChannels), iMeanGS(iNrOfChannels), iCptGS(iNrOfChannels);
    
    
    // DECLARE PIXEL_TYPE
    std::vector<pixel_type> pSumU(iNrOfChannels), pSumUold(iNrOfChannels), pDeltaSumU(iNrOfChannels);
    
    pixel_type fct1, fct2, fctST, fG, fDxu, fDyu, fDzu, fs, fTemp;
    pixel_type fct1b, fct2b, fct1e, fct2e, fct1c, fct2c, fInvMu, fInvMu2, f1, f2, f3, fNyxz;
    pixel_type fDiffNew2, fSumU, fError, fSumDiff, fSumUold, fDiffFirst2;
    
    
    // DECLARE BASIC_IMAGE
    std::vector<basic_image<pixel_type,dimension> > uOld(iNrOfChannels,basic_image<pixel_type,dimension>(geo));
    basic_image<pixel_type,dimension> uOld2(geo);
    
    image::basic_image<pixel_type,dimension> pfGb(geo);
    image::basic_image<pixel_type,dimension> pfGb2(geo);
    std::vector<image::basic_image<pixel_type,dimension> > pfHr(iNrOfChannels,basic_image<pixel_type,dimension>(geo));
    
    
    // DECLARE STD::VECTOR<BASIC_IMAGE>
    std::vector<image::basic_image<pixel_type,dimension> > pfdx(iNrOfChannels,basic_image<pixel_type,dimension>(geo));
    std::vector<image::basic_image<pixel_type,dimension> > pfdy(iNrOfChannels,basic_image<pixel_type,dimension>(geo));
    std::vector<image::basic_image<pixel_type,dimension> > pfdz(iNrOfChannels,basic_image<pixel_type,dimension>(geo));
    std::vector<image::basic_image<pixel_type,dimension> > pfbx(iNrOfChannels,basic_image<pixel_type,dimension>(geo));
    std::vector<image::basic_image<pixel_type,dimension> > pfby(iNrOfChannels,basic_image<pixel_type,dimension>(geo));
    std::vector<image::basic_image<pixel_type,dimension> > pfbz(iNrOfChannels,basic_image<pixel_type,dimension>(geo));
    

    // INITIALIZE INT
    std::fill(iMeanGS.begin(),iMeanGS.end(),0);
    std::fill(iCptGS.begin(),iCptGS.end(),0);
    
    
    // INITIALIZE PIXEL_TYPE
    fInvMu = 1./ pMu;
    fInvMu2 = SQR(fInvMu);
    fct1 = 1./6.;
    fct2 = pLambda/(6.0*pMu);
    fct1b = 1./5.;
    fct2b = pLambda/(5.0*pMu);
    fct1e = 1./4.;
    fct2e = pLambda/(4.0*pMu);
    fct1c = 1./3.;
    fct2c = pLambda/(3.0*pMu);
    fNyxz = (pixel_type)(iNy*iNx*iNz);
    
    
    // INITIALIZE BASIC_IMAGE
    for(int l=0; l<iNrOfChannels; l++){
        std::fill(pfdx[l].begin(),pfdx[l].end(),0.0);
        std::fill(pfdy[l].begin(),pfdy[l].end(),0.0);
        std::fill(pfdz[l].begin(),pfdz[l].end(),0.0);
        
        std::fill(pfbx[l].begin(),pfbx[l].end(),0.0);
        std::fill(pfby[l].begin(),pfby[l].end(),0.0);
        std::fill(pfbz[l].begin(),pfbz[l].end(),0.0);
    }
    
    std::fill(pfGb.begin(),pfGb.end(),1.0);
    std::fill(pfGb2.begin(),pfGb2.end(),1.0);
    
    
    // INITIALIZE FORCE TERM
    for(int l=0; l<iNrOfChannels; l++){
        
        std::fill(pfHr[l].begin(),pfHr[l].end(),0.0);
       
        for (int iz=0; iz< iNz; iz++)
            for (int iy=0; iy< iNy; iy++)
                for (int ix=0; ix< iNx; ix++)
                {
                        if(alpha[X3(ix,iy,iz)]<pMin) alpha[X3(ix,iy,iz)] = pMin;
                        pixel_type beta = 1.0-alpha[X3(ix,iy,iz)];
                        if(beta<pMin) beta = pMin;
                        
                        pixel_type total = pCoupling*std::log(beta) - pCoupling*std::log(alpha[X3(ix,iy,iz)]);
                        pfHr[l][X3(ix,iy,iz)] = total;
                    
                        for(int t=0; t<iNrOfLesionStates; t++){
                            
                            pixel_type total = totalIn[t][X3(ix,iy,iz)]*totalOut[t][X3(ix,iy,iz)];
                            
                            if((iInclusion_toWhich>=0) && (inclusion[X3(ix,iy,iz)]<0.1)){
                                if(t==0) total = 0.9;
                                else if(t==iNrOfLesionStates-1) total = 0.1;
                                else total = 0.0;
                            }
                            
                            if(total<pMin) total=pMin;
                            
                            if(lesionIndicator[t][l]==0) pfHr[l][X3(ix,iy,iz)] += std::log(total);
                            if(lesionIndicator[t][l]==1) pfHr[l][X3(ix,iy,iz)] -= std::log(total);
                        }
                }
    }
    
    
    // START SEG LOOP
    std::cout << "-------------------------SEG----------------------------" << std::endl;
    
    int k=0;
    while(k<iMaxIter){
        

        // STORE U_OLD
        for(int l=0; l<iNrOfChannels; l++){
            
            for (int iz=0; iz< iNz; iz++)
                for (int iy=0; iy< iNy; iy++)
                    for (int ix=0; ix< iNx; ix++)
                    {
                        uOld[l][X3(ix,iy,iz)] = u[l][X3(ix,iy,iz)];
                    }
        }
        
        
        // UPDATE SEGMENTATION
        for(int l=0; l<iNrOfChannels; l++){
        
            iGS[l]=0;
            fError = 1e10;
            while ( fError>1e-2 && iGS[l]<50 )
            {
                
                for (iz=1; iz< iNz-1; iz++)
                    for (iy=1; iy< iNy-1; iy++)
                        for (ix=1; ix< iNx-1; ix++)
                        {
                            uOld2[X3(ix,iy,iz)] = u[l][X3(ix,iy,iz)];
                        }
                
                /* Center */ // *1center - each 6 neighbours*
                for (iz=1; iz< iNz-1; iz++)
                    for (iy=1; iy< iNy-1; iy++)
                        for (ix=1; ix< iNx-1; ix++)
                        {
                            iX = X3(ix,iy,iz);
                            fG = u[l][X3(ix+1,iy,iz)] + u[l][X3(ix-1,iy,iz)] + u[l][X3(ix,iy+1,iz)] + u[l][X3(ix,iy-1,iz)] + u[l][X3(ix,iy,iz+1)] + u[l][X3(ix,iy,iz-1)];
                            fG += pfdx[l][X3(ix-1,iy,iz)] - pfdx[l][iX] - pfbx[l][X3(ix-1,iy,iz)] + pfbx[l][iX];
                            fG += pfdy[l][X3(ix,iy-1,iz)] - pfdy[l][iX] - pfby[l][X3(ix,iy-1,iz)] + pfby[l][iX];
                            fG += pfdz[l][X3(ix,iy,iz-1)] - pfdz[l][iX] - pfbz[l][X3(ix,iy,iz-1)] + pfbz[l][iX];
                            fG *= fct1;
                            fG -= fct2* pfHr[l][iX];
                            if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                            u[l][iX] = fG;
                        }
                
                /* Borders */ // *6 borders - each 5 neighbours*
                ix=0;
                for (iz=1; iz< iNz-1; iz++)
                    for (iy=1; iy< iNy-1; iy++)
                    {
                        iX = X3(ix,iy,iz);
                        fG = u[l][X3(ix+1,iy,iz)] + u[l][X3(ix,iy+1,iz)] + u[l][X3(ix,iy-1,iz)] + u[l][X3(ix,iy,iz+1)] + u[l][X3(ix,iy,iz-1)];
                        fG += - pfdx[l][iX] + pfbx[l][iX];
                        fG += pfdy[l][X3(ix,iy-1,iz)] - pfdy[l][iX] - pfby[l][X3(ix,iy-1,iz)] + pfby[l][iX];
                        fG += pfdz[l][X3(ix,iy,iz-1)] - pfdz[l][iX] - pfbz[l][X3(ix,iy,iz-1)] + pfbz[l][iX];
                        fG *= fct1b;
                        fG -= fct2b* pfHr[l][iX];
                        if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                        u[l][iX] = fG;
                    }
                
                ix=iNx-1;
                for (iz=1; iz< iNz-1; iz++)
                    for (iy=1; iy< iNy-1; iy++)
                    {
                        iX = X3(ix,iy,iz);
                        fG = u[l][X3(ix-1,iy,iz)] + u[l][X3(ix,iy+1,iz)] + u[l][X3(ix,iy-1,iz)] + u[l][X3(ix,iy,iz+1)] + u[l][X3(ix,iy,iz-1)];
                        fG += pfdx[l][X3(ix-1,iy,iz)] - pfdx[l][iX] - pfbx[l][X3(ix-1,iy,iz)] + pfbx[l][iX];
                        fG += pfdy[l][X3(ix,iy-1,iz)] - pfdy[l][iX] - pfby[l][X3(ix,iy-1,iz)] + pfby[l][iX];
                        fG += pfdz[l][X3(ix,iy,iz-1)] - pfdz[l][iX] - pfbz[l][X3(ix,iy,iz-1)] + pfbz[l][iX];
                        fG *= fct1b;
                        fG -= fct2b* pfHr[l][iX];
                        if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                        u[l][iX] = fG;
                    }
                
                iy=0;
                for (iz=1; iz< iNz-1; iz++)
                    for (ix=1; ix< iNx-1; ix++)
                    {
                        iX = X3(ix,iy,iz);
                        fG = u[l][X3(ix+1,iy,iz)] + u[l][X3(ix-1,iy,iz)] + u[l][X3(ix,iy+1,iz)] + u[l][X3(ix,iy,iz+1)] + u[l][X3(ix,iy,iz-1)];
                        fG += pfdx[l][X3(ix-1,iy,iz)] - pfdx[l][iX] - pfbx[l][X3(ix-1,iy,iz)] + pfbx[l][iX];
                        fG += - pfdy[l][iX] + pfby[l][iX];
                        fG += pfdz[l][X3(ix,iy,iz-1)] - pfdz[l][iX] - pfbz[l][X3(ix,iy,iz-1)] + pfbz[l][iX];
                        fG *= fct1b;
                        fG -= fct2b* pfHr[l][iX];
                        if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                        u[l][iX] = fG;
                    }
                
                iy=iNy-1;
                for (iz=1; iz< iNz-1; iz++)
                    for (ix=1; ix< iNx-1; ix++)
                    {
                        iX = X3(ix,iy,iz);
                        
                        fG = u[l][X3(ix+1,iy,iz)] + u[l][X3(ix-1,iy,iz)] + u[l][X3(ix,iy-1,iz)] + u[l][X3(ix,iy,iz+1)] + u[l][X3(ix,iy,iz-1)];
                        fG += pfdx[l][X3(ix-1,iy,iz)] - pfdx[l][iX] - pfbx[l][X3(ix-1,iy,iz)] + pfbx[l][iX];
                        fG += pfdy[l][X3(ix,iy-1,iz)] - pfdy[l][iX] - pfby[l][X3(ix,iy-1,iz)] + pfby[l][iX];
                        fG += pfdz[l][X3(ix,iy,iz-1)] - pfdz[l][iX] - pfbz[l][X3(ix,iy,iz-1)] + pfbz[l][iX];
                        fG *= fct1b;
                        fG -= fct2b* pfHr[l][iX];
                        if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                        u[l][iX] = fG;
                    }
                
                iz=0;
                for (iy=1; iy< iNy-1; iy++)
                    for (ix=1; ix< iNx-1; ix++)
                    {
                        iX = X3(ix,iy,iz);
                        fG = u[l][X3(ix+1,iy,iz)] + u[l][X3(ix-1,iy,iz)] + u[l][X3(ix,iy+1,iz)] + u[l][X3(ix,iy-1,iz)] + u[l][X3(ix,iy,iz+1)];
                        fG += pfdx[l][X3(ix-1,iy,iz)] - pfdx[l][iX] - pfbx[l][X3(ix-1,iy,iz)] + pfbx[l][iX];
                        fG += pfdy[l][X3(ix,iy-1,iz)] - pfdy[l][iX] - pfby[l][X3(ix,iy-1,iz)] + pfby[l][iX];
                        fG += - pfdz[l][iX] + pfbz[l][iX];
                        fG *= fct1b;
                        fG -= fct2b* pfHr[l][iX];
                        if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                        u[l][iX] = fG;
                    }
                
                iz=iNz-1;
                for (iy=1; iy< iNy-1; iy++)
                    for (ix=1; ix< iNx-1; ix++)
                    {
                        iX = X3(ix,iy,iz);
                        fG = u[l][X3(ix+1,iy,iz)] + u[l][X3(ix-1,iy,iz)] + u[l][X3(ix,iy-1,iz)] + u[l][X3(ix,iy+1,iz)] + u[l][X3(ix,iy,iz-1)] ;
                        fG += pfdx[l][X3(ix-1,iy,iz)] - pfdx[l][iX] - pfbx[l][X3(ix-1,iy,iz)] + pfbx[l][iX];
                        fG += pfdy[l][X3(ix,iy-1,iz)] - pfdy[l][iX] - pfby[l][X3(ix,iy-1,iz)] + pfby[l][iX];
                        fG += pfdz[l][X3(ix,iy,iz-1)] - pfdz[l][iX] - pfbz[l][X3(ix,iy,iz-1)] + pfbz[l][iX];
                        fG *= fct1b;
                        fG -= fct2b* pfHr[l][iX];
                        if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                        u[l][iX] = fG;
                    }
                
                /* edges */ // *12 edges - each 4 neighbours*
                ix=0; iy=0;
                for (iz=1; iz< iNz-1; iz++)
                {
                    iX = X3(ix,iy,iz);
                    
                    fG = u[l][X3(ix+1,iy,iz)] + u[l][X3(ix,iy+1,iz)] + u[l][X3(ix,iy,iz+1)] + u[l][X3(ix,iy,iz-1)];
                    fG += - pfdx[l][iX] + pfbx[l][iX];
                    fG += - pfdy[l][iX] + pfby[l][iX];
                    fG += pfdz[l][X3(ix,iy,iz-1)] - pfdz[l][iX] - pfbz[l][X3(ix,iy,iz-1)] + pfbz[l][iX];
                    fG *= fct1e;
                    fG -= fct2e* pfHr[l][iX];
                    if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                    u[l][iX] = fG;
                }
                
                ix=iNx-1; iy=0;
                for (iz=1; iz< iNz-1; iz++)
                {
                    iX = X3(ix,iy,iz);
                    
                    fG = u[l][X3(ix-1,iy,iz)] + u[l][X3(ix,iy+1,iz)] + u[l][X3(ix,iy,iz+1)] + u[l][X3(ix,iy,iz-1)];
                    fG += pfdx[l][X3(ix-1,iy,iz)] - pfdx[l][iX] - pfbx[l][X3(ix-1,iy,iz)] + pfbx[l][iX];
                    fG += - pfdy[l][iX] + pfby[l][iX];
                    fG += pfdz[l][X3(ix,iy,iz-1)] - pfdz[l][iX] - pfbz[l][X3(ix,iy,iz-1)] + pfbz[l][iX];
                    fG *= fct1e;
                    fG -= fct2e* pfHr[l][iX];
                    if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                    u[l][iX] = fG;
                }
                
                ix=0; iy=iNy-1;
                for (iz=1; iz< iNz-1; iz++)
                {
                    iX = X3(ix,iy,iz);
                    
                    fG = u[l][X3(ix+1,iy,iz)] + u[l][X3(ix,iy-1,iz)] + u[l][X3(ix,iy,iz+1)] + u[l][X3(ix,iy,iz-1)];
                    fG += - pfdx[l][iX] + pfbx[l][iX];
                    fG += pfdy[l][X3(ix,iy-1,iz)] - pfdy[l][iX] - pfby[l][X3(ix,iy-1,iz)] + pfby[l][iX];
                    fG += pfdz[l][X3(ix,iy,iz-1)] - pfdz[l][iX] - pfbz[l][X3(ix,iy,iz-1)] + pfbz[l][iX];
                    fG *= fct1e;
                    fG -= fct2e* pfHr[l][iX];
                    if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                    u[l][iX] = fG;
                }
                
                ix=iNx-1; iy=iNy-1;
                for (iz=1; iz< iNz-1; iz++)
                {
                    iX = X3(ix,iy,iz);
                    
                    fG = u[l][X3(ix-1,iy,iz)] + u[l][X3(ix,iy-1,iz)] + u[l][X3(ix,iy,iz+1)] + u[l][X3(ix,iy,iz-1)];
                    fG += pfdx[l][X3(ix-1,iy,iz)] - pfdx[l][iX] - pfbx[l][X3(ix-1,iy,iz)] + pfbx[l][iX];
                    fG += pfdy[l][X3(ix,iy-1,iz)] - pfdy[l][iX] - pfby[l][X3(ix,iy-1,iz)] + pfby[l][iX];
                    fG += pfdz[l][X3(ix,iy,iz-1)] - pfdz[l][iX] - pfbz[l][X3(ix,iy,iz-1)] + pfbz[l][iX];
                    fG *= fct1e;
                    fG -= fct2e* pfHr[l][iX];
                    if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                    u[l][iX] = fG;
                }
                
                iy=0; iz=0;
                for (ix=1; ix< iNx-1; ix++)
                {
                    iX = X3(ix,iy,iz);
                    
                    fG = u[l][X3(ix-1,iy,iz)] + u[l][X3(ix+1,iy,iz)] + u[l][X3(ix,iy+1,iz)] + u[l][X3(ix,iy,iz+1)];
                    fG += pfdx[l][X3(ix-1,iy,iz)] - pfdx[l][iX] - pfbx[l][X3(ix-1,iy,iz)] + pfbx[l][iX];
                    fG += - pfdy[l][iX] + pfby[l][iX];
                    fG += - pfdz[l][iX] + pfbz[l][iX];
                    fG *= fct1e;
                    fG -= fct2e* pfHr[l][iX];
                    if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                    u[l][iX] = fG;
                }
                
                iy=0; iz=iNz-1;
                for (ix=1; ix< iNx-1; ix++)
                {
                    iX = X3(ix,iy,iz);
                    
                    fG = u[l][X3(ix-1,iy,iz)] + u[l][X3(ix+1,iy,iz)] + u[l][X3(ix,iy+1,iz)] + u[l][X3(ix,iy,iz-1)];
                    fG += pfdx[l][X3(ix-1,iy,iz)] - pfdx[l][iX] - pfbx[l][X3(ix-1,iy,iz)] + pfbx[l][iX];
                    fG += - pfdy[l][iX] + pfby[l][iX];
                    fG += pfdz[l][X3(ix,iy,iz-1)] - pfdz[l][iX] - pfbz[l][X3(ix,iy,iz-1)] + pfbz[l][iX];
                    fG *= fct1e;
                    fG -= fct2e* pfHr[l][iX];
                    if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                    u[l][iX] = fG;
                }
                
                iy=iNy-1; iz=0;
                for (ix=1; ix< iNx-1; ix++)
                {
                    iX = X3(ix,iy,iz);
                    
                    fG = u[l][X3(ix-1,iy,iz)] + u[l][X3(ix+1,iy,iz)] + u[l][X3(ix,iy-1,iz)] + u[l][X3(ix,iy,iz+1)];
                    fG += pfdx[l][X3(ix-1,iy,iz)] - pfdx[l][iX] - pfbx[l][X3(ix-1,iy,iz)] + pfbx[l][iX];
                    fG += pfdy[l][X3(ix,iy-1,iz)] - pfdy[l][iX] - pfby[l][X3(ix,iy-1,iz)] + pfby[l][iX];
                    fG += - pfdz[l][iX] + pfbz[l][iX];
                    fG *= fct1e;
                    fG -= fct2e* pfHr[l][iX];
                    if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                    u[l][iX] = fG;
                }
                
                iy=iNy-1; iz=iNz-1;
                for (ix=1; ix< iNx-1; ix++)
                {
                    iX = X3(ix,iy,iz);
                    
                    fG = u[l][X3(ix-1,iy,iz)] + u[l][X3(ix+1,iy,iz)] + u[l][X3(ix,iy-1,iz)] + u[l][X3(ix,iy,iz-1)];
                    fG += pfdx[l][X3(ix-1,iy,iz)] - pfdx[l][iX] - pfbx[l][X3(ix-1,iy,iz)] + pfbx[l][iX];
                    fG += pfdy[l][X3(ix,iy-1,iz)] - pfdy[l][iX] - pfby[l][X3(ix,iy-1,iz)] + pfby[l][iX];
                    fG += pfdz[l][X3(ix,iy,iz-1)] - pfdz[l][iX] - pfbz[l][X3(ix,iy,iz-1)] + pfbz[l][iX];
                    fG *= fct1e;
                    fG -= fct2e* pfHr[l][iX];
                    if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                    u[l][iX] = fG;
                }
                
                ix=0; iz=0;
                for (iy=1; iy< iNy-1; iy++)
                {
                    iX = X3(ix,iy,iz);
                    
                    fG = u[l][X3(ix+1,iy,iz)] + u[l][X3(ix,iy-1,iz)] + u[l][X3(ix,iy+1,iz)] + u[l][X3(ix,iy,iz+1)];
                    fG += - pfdx[l][iX] + pfbx[l][iX];
                    fG += pfdy[l][X3(ix,iy-1,iz)] - pfdy[l][iX] - pfby[l][X3(ix,iy-1,iz)] + pfby[l][iX];
                    fG += - pfdz[l][iX] + pfbz[l][iX];
                    fG *= fct1e;
                    fG -= fct2e* pfHr[l][iX];
                    if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                    u[l][iX] = fG;
                }
                
                ix=iNx-1; iz=0;
                for (iy=1; iy< iNy-1; iy++)
                {
                    iX = X3(ix,iy,iz);
                    
                    fG = u[l][X3(ix-1,iy,iz)] + u[l][X3(ix,iy-1,iz)] + u[l][X3(ix,iy+1,iz)] + u[l][X3(ix,iy,iz+1)];
                    fG += pfdx[l][X3(ix-1,iy,iz)] - pfdx[l][iX] - pfbx[l][X3(ix-1,iy,iz)] + pfbx[l][iX];
                    fG += pfdy[l][X3(ix,iy-1,iz)] - pfdy[l][iX] - pfby[l][X3(ix,iy-1,iz)] + pfby[l][iX];
                    fG += - pfdz[l][iX] + pfbz[l][iX];
                    fG *= fct1e;
                    fG -= fct2e* pfHr[l][iX];
                    if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                    u[l][iX] = fG;
                }
                
                ix=0; iz=iNz-1;
                for (iy=1; iy< iNy-1; iy++)
                {
                    iX = X3(ix,iy,iz);
                    
                    fG = u[l][X3(ix+1,iy,iz)] + u[l][X3(ix,iy-1,iz)] + u[l][X3(ix,iy+1,iz)] + u[l][X3(ix,iy,iz-1)];
                    fG += - pfdx[l][iX] + pfbx[l][iX];
                    fG += pfdy[l][X3(ix,iy-1,iz)] - pfdy[l][iX] - pfby[l][X3(ix,iy-1,iz)] + pfby[l][iX];
                    fG += pfdz[l][X3(ix,iy,iz-1)] - pfdz[l][iX] - pfbz[l][X3(ix,iy,iz-1)] + pfbz[l][iX];
                    fG *= fct1e;
                    fG -= fct2e* pfHr[l][iX];
                    if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                    u[l][iX] = fG;
                }
                
                ix=iNx-1; iz=iNz-1;
                for (iy=1; iy< iNy-1; iy++)
                {
                    iX = X3(ix,iy,iz);
                    
                    fG = u[l][X3(ix-1,iy,iz)] + u[l][X3(ix,iy-1,iz)] + u[l][X3(ix,iy+1,iz)] + u[l][X3(ix,iy,iz-1)];
                    fG += pfdx[l][X3(ix-1,iy,iz)] - pfdx[l][iX] - pfbx[l][X3(ix-1,iy,iz)] + pfbx[l][iX];
                    fG += pfdy[l][X3(ix,iy-1,iz)] - pfdy[l][iX] - pfby[l][X3(ix,iy-1,iz)] + pfby[l][iX];
                    fG += pfdz[l][X3(ix,iy,iz-1)] - pfdz[l][iX] - pfbz[l][X3(ix,iy,iz-1)] + pfbz[l][iX];
                    fG *= fct1e;
                    fG -= fct2e* pfHr[l][iX];
                    if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                    u[l][iX] = fG;
                }
                
                /* Corners */ // *8 corners - each 3 neighbours*
                ix=0; iy=0; iz=0;
                iX = X3(ix,iy,iz);
                fG = u[l][X3(ix+1,iy,iz)] + u[l][X3(ix,iy+1,iz)] + u[l][X3(ix,iy,iz+1)];
                fG += - pfdx[l][iX] + pfbx[l][iX];
                fG += - pfdy[l][iX] + pfby[l][iX];
                fG += - pfdz[l][iX] + pfbz[l][iX];
                fG *= fct1c;
                fG -= fct2c* pfHr[l][iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                u[l][iX] = fG;
                
                ix=iNx-1; iy=0; iz=0;
                iX = X3(ix,iy,iz);
                fG = u[l][X3(ix-1,iy,iz)] + u[l][X3(ix,iy+1,iz)] + u[l][X3(ix,iy,iz+1)];
                fG += pfdx[l][X3(ix-1,iy,iz)] - pfdx[l][iX] - pfbx[l][X3(ix-1,iy,iz)] + pfbx[l][iX];
                fG += - pfdy[l][iX] + pfby[l][iX];
                fG += - pfdz[l][iX] + pfbz[l][iX];
                fG *= fct1c;
                fG -= fct2c* pfHr[l][iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                u[l][iX] = fG;
                
                ix=0; iy=iNy-1; iz=0;
                iX = X3(ix, iy, iz);
                fG = u[l][X3(ix+1, iy, iz)] + u[l][X3(ix, iy-1, iz)] + u[l][X3(ix, iy, iz+1)];
                fG += - pfdx[l][iX] + pfbx[l][iX];
                fG += pfdy[l][X3(ix,iy-1,iz)] - pfdy[l][iX] - pfby[l][X3(ix, iy-1,iz)] + pfby[l][iX];
                fG += - pfdz[l][iX] + pfbz[l][iX];
                fG *= fct1c;
                fG -= fct2c* pfHr[l][iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                u[l][iX] = fG;
                
                ix=iNx-1; iy=iNy-1; iz=0;
                iX = X3(ix,iy,iz);
                fG = u[l][X3(ix-1,iy,iz)] + u[l][X3(ix,iy-1,iz)]+ u[l][X3(ix, iy, iz+1)];
                fG += pfdx[l][X3(ix-1,iy,iz)] - pfdx[l][iX] - pfbx[l][X3(ix-1,iy,iz)] + pfbx[l][iX];
                fG += pfdy[l][X3(ix,iy-1,iz)] - pfdy[l][iX] - pfby[l][X3(ix,iy-1,iz)] + pfby[l][iX];
                fG += - pfdz[l][iX] + pfbz[l][iX];
                fG *= fct1c;
                fG -= fct2c* pfHr[l][iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                u[l][iX] = fG;
                
                ix=0; iy=0; iz=iNz-1;
                iX = X3(ix,iy,iz);
                fG = u[l][X3(ix+1,iy,iz)] + u[l][X3(ix,iy+1,iz)] + u[l][X3(ix,iy,iz-1)];
                fG += - pfdx[l][iX] + pfbx[l][iX];
                fG += - pfdy[l][iX] + pfby[l][iX];
                fG += pfdz[l][X3(ix,iy,iz-1)] - pfdz[l][iX] - pfbz[l][X3(ix,iy,iz-1)] + pfbz[l][iX];
                fG *= fct1c;
                fG -= fct2c* pfHr[l][iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                u[l][iX] = fG;
                
                ix=iNx-1; iy=0; iz=iNz-1;
                iX = X3(ix,iy,iz);
                fG = u[l][X3(ix-1,iy,iz)] + u[l][X3(ix,iy+1,iz)] + u[l][X3(ix,iy,iz-1)];
                fG += pfdx[l][X3(ix-1,iy,iz)] - pfdx[l][iX] - pfbx[l][X3(ix-1,iy,iz)] + pfbx[l][iX];
                fG += - pfdy[l][iX] + pfby[l][iX];
                fG += pfdz[l][X3(ix,iy,iz-1)] - pfdz[l][iX] - pfbz[l][X3(ix,iy,iz-1)] + pfbz[l][iX];
                fG *= fct1c;
                fG -= fct2c* pfHr[l][iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                u[l][iX] = fG;
                
                ix=0; iy=iNy-1; iz=iNz-1;
                iX = X3(ix, iy, iz);
                fG = u[l][X3(ix+1, iy, iz)] + u[l][X3(ix, iy-1, iz)] + u[l][X3(ix, iy, iz-1)];
                fG += - pfdx[l][iX] + pfbx[l][iX];
                fG += pfdy[l][X3(ix,iy-1,iz)] - pfdy[l][iX] - pfby[l][X3(ix, iy-1,iz)] + pfby[l][iX];
                fG += pfdz[l][X3(ix,iy,iz-1)] - pfdz[l][iX] - pfbz[l][X3(ix,iy,iz-1)] + pfbz[l][iX];
                fG *= fct1c;
                fG -= fct2c* pfHr[l][iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                u[l][iX] = fG;
                
                ix=iNx-1; iy=iNy-1; iz=iNz-1;
                iX = X3(ix,iy,iz);
                fG = u[l][X3(ix-1,iy,iz)] + u[l][X3(ix,iy-1,iz)]+ u[l][X3(ix, iy, iz-1)];
                fG += pfdx[l][X3(ix-1,iy,iz)] - pfdx[l][iX] - pfbx[l][X3(ix-1,iy,iz)] + pfbx[l][iX];
                fG += pfdy[l][X3(ix,iy-1,iz)] - pfdy[l][iX] - pfby[l][X3(ix,iy-1,iz)] + pfby[l][iX];
                fG += pfdz[l][X3(ix,iy,iz-1)] - pfdz[l][iX] - pfbz[l][X3(ix,iy,iz-1)] + pfbz[l][iX];
                fG *= fct1c;
                fG -= fct2c* pfHr[l][iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                u[l][iX] = fG;
                /* end Borders */
                
                
                /* Compute diff ( u - uold ) */
                fSumDiff = 0.0;
                fSumU = 0.0;
                fSumUold = 0.0;
                for (iz=0; iz< iNz; iz++)
                    for (iy=0; iy< iNy; iy++)
                        for (ix=0; ix< iNx; ix++)
                            fSumDiff += SQR(u[l][X3(ix,iy,iz)]-uOld2[X3(ix,iy,iz)]);
                fDiffNew2 = fSumDiff/ fNyxz;
                if ( iGS[l]==0 )
                {
                    fDiffFirst2 = fDiffNew2;
                    fDiffNew2 = 1e10;
                    fError = 1e10;
                }
                else
                    fError = 1.0 - std::abs(fDiffNew2-fDiffFirst2)/fDiffFirst2;
                iGS[l]++;
            } //End gauss seidel
            
            iMeanGS[l] += iGS[l];
            iCptGS[l]++;
            
            /* Center */
            for (iz=0; iz< iNz-1; iz++)
                for (iy=0; iy< iNy-1; iy++)
                    for (ix=0; ix< iNx-1; ix++)
                    {
                        iX = X3(ix,iy,iz);
                        /* d */
                        fDxu = u[l][X3(ix+1,iy,iz)] - u[l][iX];
                        fDyu = u[l][X3(ix,iy+1,iz)] - u[l][iX];
                        fDzu = u[l][X3(ix,iy,iz+1)] - u[l][iX];
                        f1 = fDxu+pfbx[l][iX];
                        f2 = fDyu+pfby[l][iX];
                        f3 = fDzu+pfbz[l][iX];
                        fs = SQR(f1)+SQR(f2)+SQR(f3);
                        fctST = fInvMu2* pfGb2[iX];
                        if ( fs<fctST ) { pfdx[l][iX]=0.0; pfdy[l][iX]=0.0; pfdz[l][iX]=0.0; }
                        else {
                            fs = std::sqrt(fs);
                            fctST = std::sqrt(fctST);
                            fTemp = fs-fctST; fTemp /= fs;
                            pfdx[l][iX] = fTemp* f1;
                            pfdy[l][iX] = fTemp* f2;
                            pfdz[l][iX] = fTemp* f3;}
                        /* b */
                        pfbx[l][iX] += fDxu - pfdx[l][iX];
                        pfby[l][iX] += fDyu - pfdy[l][iX];
                        pfbz[l][iX] += fDzu - pfdz[l][iX];
                    }
            
            /* Borders */
            ix=iNx-1;
            for (iz=1; iz< iNz-1; iz++)
                for (iy=1; iy< iNy-1; iy++)
                {
                    iX = X3(ix,iy,iz);
                    fDxu = 0.0;
                    fDyu = u[l][X3(ix,iy+1,iz)] - u[l][iX];
                    fDzu = u[l][X3(ix,iy,iz+1)] - u[l][iX];
                    f1 = fDxu+pfbx[l][iX]; f2 = fDyu+pfby[l][iX]; f3 = fDzu+pfbz[l][iX];
                    fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
                    if ( fs<fctST ) { pfdx[l][iX]=0.0; pfdy[l][iX]=0.0; pfdz[l][iX]=0.0; }
                    else {
                        fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                        fTemp = fs-fctST; fTemp /= fs;
                        pfdx[l][iX] = fTemp* f1;
                        pfdy[l][iX] = fTemp* f2;
                        pfdz[l][iX] = fTemp* f3;}
                    pfbx[l][iX] += fDxu - pfdx[l][iX];
                    pfby[l][iX] += fDyu - pfdy[l][iX];
                    pfbz[l][iX] += fDzu - pfdz[l][iX];
                }
            
            iy=iNy-1;
            for (iz=1; iz< iNz-1; iz++)
                for (ix=1; ix< iNx-1; ix++)
                {
                    iX = X3(ix,iy,iz);
                    fDxu = u[l][X3(ix+1,iy,iz)] - u[l][iX];
                    fDyu = 0.0;
                    fDzu = u[l][X3(ix,iy,iz+1)] - u[l][iX];
                    f1 = fDxu+pfbx[l][iX]; f2 = fDyu+pfby[l][iX]; f3 = fDzu+pfbz[l][iX];
                    fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
                    if ( fs<fctST ) { pfdx[l][iX]=0.0; pfdy[l][iX]=0.0; pfdz[l][iX]=0.0; }
                    else {
                        fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                        fTemp = fs-fctST; fTemp /= fs;
                        pfdx[l][iX] = fTemp* f1;
                        pfdy[l][iX] = fTemp* f2;
                        pfdz[l][iX] = fTemp* f3;}
                    pfbx[l][iX] += fDxu - pfdx[l][iX];
                    pfby[l][iX] += fDyu - pfdy[l][iX];
                    pfbz[l][iX] += fDzu - pfdz[l][iX];
                }
            
            iz=iNz-1;
            for (iy=1; iy< iNy-1; iy++)
                for (ix=1; ix< iNx-1; ix++)
                {
                    iX = X3(ix,iy,iz);
                    fDxu = u[l][X3(ix+1,iy,iz)] - u[l][iX];
                    fDyu = u[l][X3(ix,iy+1,iz)] - u[l][iX];
                    fDzu = 0.0;
                    f1 = fDxu+pfbx[l][iX]; f2 = fDyu+pfby[l][iX]; f3 = fDzu+pfbz[l][iX];
                    fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
                    if ( fs<fctST ) { pfdx[l][iX]=0.0; pfdy[l][iX]=0.0; pfdz[l][iX]=0.0; }
                    else {
                        fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                        fTemp = fs-fctST; fTemp /= fs;
                        pfdx[l][iX] = fTemp* f1;
                        pfdy[l][iX] = fTemp* f2;
                        pfdz[l][iX] = fTemp* f3;}
                    pfbx[l][iX] += fDxu - pfdx[l][iX];
                    pfby[l][iX] += fDyu - pfdy[l][iX];
                    pfbz[l][iX] += fDzu - pfdz[l][iX];
                }
            
            /*edges*/
            ix=iNx-1; iz=0;
            for (iy=1; iy< iNy-1; iy++)
            {
                iX = X3(ix,iy,iz);
                fDxu = 0.0;
                fDyu = u[l][X3(ix,iy+1,iz)] - u[l][iX];
                fDzu = u[l][X3(ix,iy,iz+1)] - u[l][iX];
                f1 = fDxu+pfbx[l][iX]; f2 = fDyu+pfby[l][iX]; f3 = fDzu+pfbz[l][iX];
                fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
                if ( fs<fctST ) { pfdx[l][iX]=0.0; pfdy[l][iX]=0.0; pfdz[l][iX]=0.0; }
                else {
                    fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                    fTemp = fs-fctST; fTemp /= fs;
                    pfdx[l][iX] = fTemp* f1;
                    pfdy[l][iX] = fTemp* f2;
                    pfdz[l][iX] = fTemp* f3;}
                pfbx[l][iX] += fDxu - pfdx[l][iX];
                pfby[l][iX] += fDyu - pfdy[l][iX];
                pfbz[l][iX] += fDzu - pfdz[l][iX];
            }
            
            ix=0; iz=iNz-1;
            for (iy=1; iy< iNy-1; iy++)
            {
                iX = X3(ix,iy,iz);
                fDxu = u[l][X3(ix+1,iy,iz)] - u[l][iX];
                fDyu = u[l][X3(ix,iy+1,iz)] - u[l][iX];
                fDzu = 0.0;
                f1 = fDxu+pfbx[l][iX]; f2 = fDyu+pfby[l][iX]; f3 = fDzu+pfbz[l][iX];
                fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
                if ( fs<fctST ) { pfdx[l][iX]=0.0; pfdy[l][iX]=0.0; pfdz[l][iX]=0.0; }
                else {
                    fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                    fTemp = fs-fctST; fTemp /= fs;
                    pfdx[l][iX] = fTemp* f1;
                    pfdy[l][iX] = fTemp* f2;
                    pfdz[l][iX] = fTemp* f3;}
                pfbx[l][iX] += fDxu - pfdx[l][iX];
                pfby[l][iX] += fDyu - pfdy[l][iX];
                pfbz[l][iX] += fDzu - pfdz[l][iX];
            }
            
            ix=iNx-1; iz=iNz-1;
            for (iy=1; iy< iNy-1; iy++)
            {
                iX = X3(ix,iy,iz);
                fDxu = 0.0;
                fDyu = u[l][X3(ix,iy+1,iz)] - u[l][iX];
                fDzu = 0.0;
                f1 = fDxu+pfbx[l][iX]; f2 = fDyu+pfby[l][iX]; f3 = fDzu+pfbz[l][iX];
                fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
                if ( fs<fctST ) { pfdx[l][iX]=0.0; pfdy[l][iX]=0.0; pfdz[l][iX]=0.0; }
                else {
                    fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                    fTemp = fs-fctST; fTemp /= fs;
                    pfdx[l][iX] = fTemp* f1;
                    pfdy[l][iX] = fTemp* f2;
                    pfdz[l][iX] = fTemp* f3;}
                pfbx[l][iX] += fDxu - pfdx[l][iX];
                pfby[l][iX] += fDyu - pfdy[l][iX];
                pfbz[l][iX] += fDzu - pfdz[l][iX];
            }
            
            
            iy=0; iz=iNz-1;
            for (ix=1; ix< iNx-1; ix++)
            {
                iX = X3(ix,iy,iz);
                fDxu = u[l][X3(ix+1,iy,iz)] - u[l][iX];
                fDyu = u[l][X3(ix,iy+1,iz)] - u[l][iX];
                fDzu = 0.0;
                f1 = fDxu+pfbx[l][iX]; f2 = fDyu+pfby[l][iX]; f3 = fDzu+pfbz[l][iX];
                fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
                if ( fs<fctST ) { pfdx[l][iX]=0.0; pfdy[l][iX]=0.0; pfdz[l][iX]=0.0; }
                else {
                    fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                    fTemp = fs-fctST; fTemp /= fs;
                    pfdx[l][iX] = fTemp* f1;
                    pfdy[l][iX] = fTemp* f2;
                    pfdz[l][iX] = fTemp* f3;}
                pfbx[l][iX] += fDxu - pfdx[l][iX];
                pfby[l][iX] += fDyu - pfdy[l][iX];
                pfbz[l][iX] += fDzu - pfdz[l][iX];
            }
            
            iy=iNy-1; iz=0;
            for (ix=1; ix< iNx-1; ix++)
            {
                iX = X3(ix,iy,iz);
                fDxu = u[l][X3(ix+1,iy,iz)] - u[l][iX];
                fDyu = 0.0;
                fDzu = u[l][X3(ix,iy,iz+1)] - u[l][iX];
                f1 = fDxu+pfbx[l][iX]; f2 = fDyu+pfby[l][iX]; f3 = fDzu+pfbz[l][iX];
                fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
                if ( fs<fctST ) { pfdx[l][iX]=0.0; pfdy[l][iX]=0.0; pfdz[l][iX]=0.0; }
                else {
                    fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                    fTemp = fs-fctST; fTemp /= fs;
                    pfdx[l][iX] = fTemp* f1;
                    pfdy[l][iX] = fTemp* f2;
                    pfdz[l][iX] = fTemp* f3;}
                pfbx[l][iX] += fDxu - pfdx[l][iX];
                pfby[l][iX] += fDyu - pfdy[l][iX];
                pfbz[l][iX] += fDzu - pfdz[l][iX];
            }
            
            iy=iNy-1; iz=iNz-1;
            for (ix=1; ix< iNx-1; ix++)
            {
                iX = X3(ix,iy,iz);
                fDxu = u[l][X3(ix+1,iy,iz)] - u[l][iX];
                fDyu = 0.0;
                fDzu = 0.0;
                f1 = fDxu+pfbx[l][iX]; f2 = fDyu+pfby[l][iX]; f3 = fDzu+pfbz[l][iX];
                fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
                if ( fs<fctST ) { pfdx[l][iX]=0.0; pfdy[l][iX]=0.0; pfdz[l][iX]=0.0; }
                else {
                    fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                    fTemp = fs-fctST; fTemp /= fs;
                    pfdx[l][iX] = fTemp* f1;
                    pfdy[l][iX] = fTemp* f2;
                    pfdz[l][iX] = fTemp* f3;}
                pfbx[l][iX] += fDxu - pfdx[l][iX];
                pfby[l][iX] += fDyu - pfdy[l][iX];
                pfbz[l][iX] += fDzu - pfdz[l][iX];
            }
            
            ix=iNx-1; iy=0;
            for (iz=1; iz< iNz-1; iz++)
            {
                iX = X3(ix,iy,iz);
                fDxu = 0.0;
                fDyu = u[l][X3(ix,iy+1,iz)] - u[l][iX];
                fDzu = u[l][X3(ix,iy,iz+1)] - u[l][iX];
                f1 = fDxu+pfbx[l][iX]; f2 = fDyu+pfby[l][iX]; f3 = fDzu+pfbz[l][iX];
                fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
                if ( fs<fctST ) { pfdx[l][iX]=0.0; pfdy[l][iX]=0.0; pfdz[l][iX]=0.0; }
                else {
                    fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                    fTemp = fs-fctST; fTemp /= fs;
                    pfdx[l][iX] = fTemp* f1;
                    pfdy[l][iX] = fTemp* f2;
                    pfdz[l][iX] = fTemp* f3;}
                pfbx[l][iX] += fDxu - pfdx[l][iX];
                pfby[l][iX] += fDyu - pfdy[l][iX];
                pfbz[l][iX] += fDzu - pfdz[l][iX];
            }
            
            ix=0; iy=iNy-1;
            for (iz=1; iz< iNz-1; iz++)
            {
                iX = X3(ix,iy,iz);
                fDxu = u[l][X3(ix+1,iy,iz)] - u[l][iX];
                fDyu = 0.0;
                fDzu = u[l][X3(ix,iy,iz+1)] - u[l][iX];
                f1 = fDxu+pfbx[l][iX]; f2 = fDyu+pfby[l][iX]; f3 = fDzu+pfbz[l][iX];
                fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
                if ( fs<fctST ) { pfdx[l][iX]=0.0; pfdy[l][iX]=0.0; pfdz[l][iX]=0.0; }
                else {
                    fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                    fTemp = fs-fctST; fTemp /= fs;
                    pfdx[l][iX] = fTemp* f1;
                    pfdy[l][iX] = fTemp* f2;
                    pfdz[l][iX] = fTemp* f3;}
                pfbx[l][iX] += fDxu - pfdx[l][iX];
                pfby[l][iX] += fDyu - pfdy[l][iX];
                pfbz[l][iX] += fDzu - pfdz[l][iX];
            }
            
            ix=iNx-1; iy=iNy-1;
            for (iz=1; iz< iNz-1; iz++)
            {
                iX = X3(ix,iy,iz);
                fDxu = 0.0;
                fDyu = 0.0;
                fDzu = u[l][X3(ix,iy,iz+1)] - u[l][iX];
                f1 = fDxu+pfbx[l][iX]; f2 = fDyu+pfby[l][iX]; f3 = fDzu+pfbz[l][iX];
                fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
                if ( fs<fctST ) { pfdx[l][iX]=0.0; pfdy[l][iX]=0.0; pfdz[l][iX]=0.0; }
                else {
                    fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                    fTemp = fs-fctST; fTemp /= fs;
                    pfdx[l][iX] = fTemp* f1;
                    pfdy[l][iX] = fTemp* f2;
                    pfdz[l][iX] = fTemp* f3;}
                pfbx[l][iX] += fDxu - pfdx[l][iX];
                pfby[l][iX] += fDyu - pfdy[l][iX];
                pfbz[l][iX] += fDzu - pfdz[l][iX];
            }
            
            /*corners*/
            ix=iNx-1; iy=0; iz=0;
            iX = X3(ix,iy,iz);
            fDxu = 0.0;
            fDyu = u[l][X3(ix,iy+1,iz)] - u[l][iX];
            fDzu = u[l][X3(ix,iy,iz+1)] - u[l][iX];
            f1 = fDxu+pfbx[l][iX]; f2 = fDyu+pfby[l][iX]; f3 = fDzu+pfbz[l][iX];
            fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
            if ( fs<fctST ) { pfdx[l][iX]=0.0; pfdy[l][iX]=0.0; pfdz[l][iX]=0.0; }
            else {
                fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                fTemp = fs-fctST; fTemp /= fs;
                pfdx[l][iX] = fTemp* f1;
                pfdy[l][iX] = fTemp* f2;
                pfdz[l][iX] = fTemp* f3;}
            pfbx[l][iX] += fDxu - pfdx[l][iX];
            pfby[l][iX] += fDyu - pfdy[l][iX];
            pfbz[l][iX] += fDzu - pfdz[l][iX];
            
            ix=0; iy=iNy-1; iz=0;
            iX = X3(ix,iy,iz);
            fDxu = u[l][X3(ix+1,iy,iz)] - u[l][iX];
            fDyu = 0.0;
            fDzu = u[l][X3(ix,iy,iz+1)] - u[l][iX];
            f1 = fDxu+pfbx[l][iX]; f2 = fDyu+pfby[l][iX]; f3 = fDzu+pfbz[l][iX];
            fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
            if ( fs<fctST ) { pfdx[l][iX]=0.0; pfdy[l][iX]=0.0; pfdz[l][iX]=0.0; }
            else {
                fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                fTemp = fs-fctST; fTemp /= fs;
                pfdx[l][iX] = fTemp* f1;
                pfdy[l][iX] = fTemp* f2;
                pfdz[l][iX] = fTemp* f3;}
            pfbx[l][iX] += fDxu - pfdx[l][iX];
            pfby[l][iX] += fDyu - pfdy[l][iX];
            pfbz[l][iX] += fDzu - pfdz[l][iX];
            
            ix=iNx-1; iy=iNy-1; iz=0;
            iX = X3(ix,iy,iz);
            fDxu = 0.0;
            fDyu = 0.0;
            fDzu = u[l][X3(ix,iy,iz+1)] - u[l][iX];
            f1 = fDxu+pfbx[l][iX]; f2 = fDyu+pfby[l][iX]; f3 = fDzu+pfbz[l][iX];
            fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
            if ( fs<fctST ) { pfdx[l][iX]=0.0; pfdy[l][iX]=0.0; pfdz[l][iX]=0.0;}
            else {
                fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                fTemp = fs-fctST; fTemp /= fs;
                pfdx[l][iX] = fTemp* f1;
                pfdy[l][iX] = fTemp* f2;
                pfdz[l][iX] = fTemp* f3;}
            pfbx[l][iX] += fDxu - pfdx[l][iX];
            pfby[l][iX] += fDyu - pfdy[l][iX];
            pfbz[l][iX] += fDzu - pfdz[l][iX];
            
            ix=iNx-1; iy=0; iz=iNz-1;
            iX = X3(ix,iy,iz);
            fDxu = 0.0;
            fDyu = u[l][X3(ix,iy+1,iz)] - u[l][iX];
            fDzu = 0.0;
            f1 = fDxu+pfbx[l][iX]; f2 = fDyu+pfby[l][iX]; f3 = fDzu+pfbz[l][iX];
            fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
            if ( fs<fctST ) { pfdx[l][iX]=0.0; pfdy[l][iX]=0.0; pfdz[l][iX]=0.0; }
            else {
                fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                fTemp = fs-fctST; fTemp /= fs;
                pfdx[l][iX] = fTemp* f1;
                pfdy[l][iX] = fTemp* f2;
                pfdz[l][iX] = fTemp* f3;}
            pfbx[l][iX] += fDxu - pfdx[l][iX];
            pfby[l][iX] += fDyu - pfdy[l][iX];
            pfbz[l][iX] += fDzu - pfdz[l][iX];
            
            ix=0; iy=iNy-1; iz=iNz-1;
            iX = X3(ix,iy,iz);
            fDxu = u[l][X3(ix+1,iy,iz)] - u[l][iX];
            fDyu = 0.0;
            fDzu = 0.0;
            f1 = fDxu+pfbx[l][iX]; f2 = fDyu+pfby[l][iX]; f3 = fDzu+pfbz[l][iX];
            fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
            if ( fs<fctST ) { pfdx[l][iX]=0.0; pfdy[l][iX]=0.0; pfdz[l][iX]=0.0; }
            else {
                fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                fTemp = fs-fctST; fTemp /= fs;
                pfdx[l][iX] = fTemp* f1;
                pfdy[l][iX] = fTemp* f2;
                pfdz[l][iX] = fTemp* f3;}
            pfbx[l][iX] += fDxu - pfdx[l][iX];
            pfby[l][iX] += fDyu - pfdy[l][iX];
            pfbz[l][iX] += fDzu - pfdz[l][iX];
            
            ix=iNx-1; iy=iNy-1; iz=iNz-1;
            iX = X3(ix,iy,iz);
            fDxu = 0.0;
            fDyu = 0.0;
            fDzu = 0.0;
            f1 = fDxu+pfbx[l][iX]; f2 = fDyu+pfby[l][iX]; f3 = fDzu+pfbz[l][iX];
            fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
            if ( fs<fctST ) { pfdx[l][iX]=0.0; pfdy[l][iX]=0.0; pfdz[l][iX]=0.0;}
            else {
                fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                fTemp = fs-fctST; fTemp /= fs;
                pfdx[l][iX] = fTemp* f1;
                pfdy[l][iX] = fTemp* f2;
                pfdz[l][iX] = fTemp* f3;}
            pfbx[l][iX] += fDxu - pfdx[l][iX];
            pfby[l][iX] += fDyu - pfdy[l][iX];
            pfbz[l][iX] += fDzu - pfdz[l][iX];
        
            ix=0; iy=0; iz=iNz-1;
            iX = X3(ix,iy,iz);
            fDxu = u[l][X3(ix+1,iy,iz)] - u[l][iX];
            fDyu = u[l][X3(ix,iy+1,iz)] - u[l][iX];
            fDzu = 0.0;
            f1 = fDxu+pfbx[l][iX]; f2 = fDyu+pfby[l][iX]; f3 = fDzu+pfbz[l][iX];
            fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
            if ( fs<fctST ) { pfdx[l][iX]=0.0; pfdy[l][iX]=0.0; pfdz[l][iX]=0.0;}
            else {
                fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                fTemp = fs-fctST; fTemp /= fs;
                pfdx[l][iX] = fTemp* f1;
                pfdy[l][iX] = fTemp* f2;
                pfdz[l][iX] = fTemp* f3;}
            pfbx[l][iX] += fDxu - pfdx[l][iX];
            pfby[l][iX] += fDyu - pfdy[l][iX];
            pfbz[l][iX] += fDzu - pfdz[l][iX];
            
            /* Borders */
        
            
            // CHECK CONVERGENCE
            pDeltaSumU[l] = 0.0;
            pSumUold[l] = 0.0;
            pSumU[l] = 0.0;
            
            for (int iz=0; iz< iNz; iz++)
                for (int iy=0; iy< iNy; iy++)
                    for (int ix=0; ix< iNx; ix++)
                    {
                        pDeltaSumU[l] += SQR(u[l][X3(ix,iy,iz)]-uOld[l][X3(ix,iy,iz)]);
                        pSumU[l] += SQR(u[l][X3(ix,iy,iz)]);
                        pSumUold[l] += SQR(uOld[l][X3(ix,iy,iz)]);
                    }
            
            pDeltaSumU[l] = pDeltaSumU[l]/ (pSumU[l]*pSumUold[l]);
        
        }
        
        
        // SCREEN FEEDBACK
        std::stringstream ss; ss << k;
        std::cout << " End SEG iteration: " << ss.str() << " | Relative change u: " << pDeltaSumU << std::endl;
        
        bool converged = true;
        for(int l=0; l<iNrOfChannels; l++) if(std::abs(pDeltaSumU[l])>pEps) converged = false;
        if(converged) break;
        k++;
        
    }
    
    // UPDATE ALPHA
    if(bUpdateAlpha){
        
        pixel_type sum1 = 0.0;
        pixel_type sum2 = 0.0;
        for(int l=0; l<iNrOfChannels; l++){
            for (int iz=0; iz< iNz; iz++)
            for (int iy=0; iy< iNy; iy++)
            for (int ix=0; ix< iNx; ix++){
                if(brainMask[X3(ix,iy,iz)]){
                    sum1 += u[l][X3(ix,iy,iz)]/(pixel_type) iNrOfChannels;
                    sum2 += 1.0;
                }
            }
        }
        std::fill(alpha.begin(),alpha.end(),sum1/sum2);
    }

    return;

}

}
}
#endif //SEG_DISJOINT_HPP_INCLUDED
