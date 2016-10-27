#ifndef SEG_JOINT_HPP_INCLUDED
#define SEG_JOINT_HPP_INCLUDED

#define SQR(x) (x)*(x)
#define X3(ix,iy,iz) ((iz)*iNy+(iy))*iNx+(ix)

namespace image
{
namespace segmentation
{

    
// The code in this function is inspired by Matlab code by Tom Goldstein & Xavier Bresson [1]
// [1] Tom Goldstein et al. "Geometric applications of the split Bregman method: segmentation and surface reconstruction." J. of Sc. Comp. 45.1-3 (2010): 272-293.
template<typename pixel_type,unsigned int dimension>
void SEG_joint(image::basic_image<pixel_type,dimension>& totalIn,
            image::basic_image<pixel_type,dimension>& totalOut,
            image::basic_image<pixel_type,dimension>& alpha,
         
            image::basic_image<pixel_type,dimension>& brainMask,
            image::basic_image<pixel_type,dimension>& u,
        
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
    geometry<dimension> geo = totalIn.geometry();
    
    int iNx = geo[0];
    int iNy = geo[1];
    int iNz = geo[2];
    
    int ix;
    int iy;
    int iz;
    int iX;
    int iGS, iMeanGS, iCptGS;
    
    
    // DECLARE PIXEL_TYPE
    pixel_type pSumU,pSumUold,pDeltaSumU;
    
    pixel_type fct1, fct2, fctST, fG, fDxu, fDyu, fDzu, fs, fTemp;
    pixel_type fct1b, fct2b, fct1e, fct2e, fct1c, fct2c, fInvMu, fInvMu2, f1, f2, f3, fNyxz;
    pixel_type fDiffNew2, fSumU, fError, fSumDiff, fSumUold, fDiffFirst2;

    
    // DECLARE BASIC_IMAGE
    basic_image<pixel_type,dimension> uOld(geo);
    basic_image<pixel_type,dimension> uOld2(geo);
    
    image::basic_image<pixel_type,dimension> pfGb(geo);
    image::basic_image<pixel_type,dimension> pfGb2(geo);
    image::basic_image<pixel_type,dimension> pfHr(geo);
    
    image::basic_image<pixel_type,dimension> pfdx(geo);
    image::basic_image<pixel_type,dimension> pfdy(geo);
    image::basic_image<pixel_type,dimension> pfdz(geo);
    image::basic_image<pixel_type,dimension> pfbx(geo);
    image::basic_image<pixel_type,dimension> pfby(geo);
    image::basic_image<pixel_type,dimension> pfbz(geo);
    
    
    // INITIALIZE INT
    iMeanGS = 0;
    iCptGS = 0;
    
    
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
    std::fill(pfdx.begin(),pfdx.end(),0.0);
    std::fill(pfdy.begin(),pfdy.end(),0.0);
    std::fill(pfdz.begin(),pfdz.end(),0.0);
    
    std::fill(pfbx.begin(),pfbx.end(),0.0);
    std::fill(pfby.begin(),pfby.end(),0.0);
    std::fill(pfbz.begin(),pfbz.end(),0.0);
    
    std::fill(pfGb.begin(),pfGb.end(),1.0);
    std::fill(pfGb2.begin(),pfGb2.end(),1.0);
    
    
    // INITIALIZE FORCE TERM
    for (int iz=0; iz< iNz; iz++)
        for (int iy=0; iy< iNy; iy++)
            for (int ix=0; ix< iNx; ix++)
            {
                if(totalOut[X3(ix,iy,iz)]<pMin) totalOut[X3(ix,iy,iz)] = pMin;
                if(totalIn[X3(ix,iy,iz)]<pMin) totalIn[X3(ix,iy,iz)] = pMin;
                if(alpha[X3(ix,iy,iz)]<pMin) alpha[X3(ix,iy,iz)] = pMin;
                
                if((iInclusion_toWhich>=0) && (inclusion[X3(ix,iy,iz)]<0.1)){
                    totalOut[X3(ix,iy,iz)] = 0.9;
                    totalIn[X3(ix,iy,iz)] = 0.1;
                }
                
                pixel_type beta = 1.0-alpha[X3(ix,iy,iz)];
                if(beta<pMin) beta = pMin;
                
                pfHr[X3(ix,iy,iz)] = -std::log(totalIn[X3(ix,iy,iz)]) + std::log(totalOut[X3(ix,iy,iz)]) + pCoupling* (std::log(beta) - std::log(alpha[X3(ix,iy,iz)]));
            }
    
    

    // START SEG LOOP
    std::cout << "-------------------------SEG----------------------------" << std::endl;
    
    int k=0;
    while(k<iMaxIter){
        

        // STORE U_OLD
        for (int iz=0; iz< iNz; iz++)
            for (int iy=0; iy< iNy; iy++)
                for (int ix=0; ix< iNx; ix++)
                {
                    uOld[X3(ix,iy,iz)] = u[X3(ix,iy,iz)];
                }
    
        
        // UPDATE SEGMENTATION
        iGS=0;
        fError = 1e10;
        while ( fError>1e-2 && iGS<50 )
        {
            
            for (iz=1; iz< iNz-1; iz++)
                for (iy=1; iy< iNy-1; iy++)
                    for (ix=1; ix< iNx-1; ix++)
                    {
                        uOld2[X3(ix,iy,iz)] = u[X3(ix,iy,iz)];
                    }

            
            
            /* Center */ // *1center - each 6 neighbours*
            for (iz=1; iz< iNz-1; iz++)
                for (iy=1; iy< iNy-1; iy++)
                    for (ix=1; ix< iNx-1; ix++)
                    {
                        iX = X3(ix,iy,iz);
                        fG = u[X3(ix+1,iy,iz)] + u[X3(ix-1,iy,iz)] + u[X3(ix,iy+1,iz)] + u[X3(ix,iy-1,iz)] + u[X3(ix,iy,iz+1)] + u[X3(ix,iy,iz-1)];
                        fG += pfdx[X3(ix-1,iy,iz)] - pfdx[iX] - pfbx[X3(ix-1,iy,iz)] + pfbx[iX];
                        fG += pfdy[X3(ix,iy-1,iz)] - pfdy[iX] - pfby[X3(ix,iy-1,iz)] + pfby[iX];
                        fG += pfdz[X3(ix,iy,iz-1)] - pfdz[iX] - pfbz[X3(ix,iy,iz-1)] + pfbz[iX];
                        fG *= fct1;
                        fG -= fct2* pfHr[iX];
                        if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                        u[iX] = fG;
                    }
            
            /* Borders */ // *6 borders - each 5 neighbours*
            ix=0;
            for (iz=1; iz< iNz-1; iz++)
                for (iy=1; iy< iNy-1; iy++)
                {
                    iX = X3(ix,iy,iz);
                    fG = u[X3(ix+1,iy,iz)] + u[X3(ix,iy+1,iz)] + u[X3(ix,iy-1,iz)] + u[X3(ix,iy,iz+1)] + u[X3(ix,iy,iz-1)];
                    fG += - pfdx[iX] + pfbx[iX];
                    fG += pfdy[X3(ix,iy-1,iz)] - pfdy[iX] - pfby[X3(ix,iy-1,iz)] + pfby[iX];
                    fG += pfdz[X3(ix,iy,iz-1)] - pfdz[iX] - pfbz[X3(ix,iy,iz-1)] + pfbz[iX];
                    fG *= fct1b;
                    fG -= fct2b* pfHr[iX];
                    if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                    u[iX] = fG;
                }
            
            ix=iNx-1;
            for (iz=1; iz< iNz-1; iz++)
                for (iy=1; iy< iNy-1; iy++)
                {
                    iX = X3(ix,iy,iz);
                    fG = u[X3(ix-1,iy,iz)] + u[X3(ix,iy+1,iz)] + u[X3(ix,iy-1,iz)] + u[X3(ix,iy,iz+1)] + u[X3(ix,iy,iz-1)];
                    fG += pfdx[X3(ix-1,iy,iz)] - pfdx[iX] - pfbx[X3(ix-1,iy,iz)] + pfbx[iX];
                    fG += pfdy[X3(ix,iy-1,iz)] - pfdy[iX] - pfby[X3(ix,iy-1,iz)] + pfby[iX];
                    fG += pfdz[X3(ix,iy,iz-1)] - pfdz[iX] - pfbz[X3(ix,iy,iz-1)] + pfbz[iX];
                    fG *= fct1b;
                    fG -= fct2b* pfHr[iX];
                    if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                    u[iX] = fG;
                }
            
            iy=0;
            for (iz=1; iz< iNz-1; iz++)
                for (ix=1; ix< iNx-1; ix++)
                {
                    iX = X3(ix,iy,iz);
                    fG = u[X3(ix+1,iy,iz)] + u[X3(ix-1,iy,iz)] + u[X3(ix,iy+1,iz)] + u[X3(ix,iy,iz+1)] + u[X3(ix,iy,iz-1)];
                    fG += pfdx[X3(ix-1,iy,iz)] - pfdx[iX] - pfbx[X3(ix-1,iy,iz)] + pfbx[iX];
                    fG += - pfdy[iX] + pfby[iX];
                    fG += pfdz[X3(ix,iy,iz-1)] - pfdz[iX] - pfbz[X3(ix,iy,iz-1)] + pfbz[iX];
                    fG *= fct1b;
                    fG -= fct2b* pfHr[iX];
                    if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                    u[iX] = fG;
                }
            
            iy=iNy-1;
            for (iz=1; iz< iNz-1; iz++)
                for (ix=1; ix< iNx-1; ix++)
                {
                    iX = X3(ix,iy,iz);
                    
                    fG = u[X3(ix+1,iy,iz)] + u[X3(ix-1,iy,iz)] + u[X3(ix,iy-1,iz)] + u[X3(ix,iy,iz+1)] + u[X3(ix,iy,iz-1)];
                    fG += pfdx[X3(ix-1,iy,iz)] - pfdx[iX] - pfbx[X3(ix-1,iy,iz)] + pfbx[iX];
                    fG += pfdy[X3(ix,iy-1,iz)] - pfdy[iX] - pfby[X3(ix,iy-1,iz)] + pfby[iX];
                    fG += pfdz[X3(ix,iy,iz-1)] - pfdz[iX] - pfbz[X3(ix,iy,iz-1)] + pfbz[iX];
                    fG *= fct1b;
                    fG -= fct2b* pfHr[iX];
                    if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                    u[iX] = fG;
                }
            
            iz=0;
            for (iy=1; iy< iNy-1; iy++)
                for (ix=1; ix< iNx-1; ix++)
                {
                    iX = X3(ix,iy,iz);
                    fG = u[X3(ix+1,iy,iz)] + u[X3(ix-1,iy,iz)] + u[X3(ix,iy+1,iz)] + u[X3(ix,iy-1,iz)] + u[X3(ix,iy,iz+1)];
                    fG += pfdx[X3(ix-1,iy,iz)] - pfdx[iX] - pfbx[X3(ix-1,iy,iz)] + pfbx[iX];
                    fG += pfdy[X3(ix,iy-1,iz)] - pfdy[iX] - pfby[X3(ix,iy-1,iz)] + pfby[iX];
                    fG += - pfdz[iX] + pfbz[iX];
                    fG *= fct1b;
                    fG -= fct2b* pfHr[iX];
                    if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                    u[iX] = fG;
                }
            
            iz=iNz-1;
            for (iy=1; iy< iNy-1; iy++)
                for (ix=1; ix< iNx-1; ix++)
                {
                    iX = X3(ix,iy,iz);
                    fG = u[X3(ix+1,iy,iz)] + u[X3(ix-1,iy,iz)] + u[X3(ix,iy-1,iz)] + u[X3(ix,iy+1,iz)] + u[X3(ix,iy,iz-1)] ;
                    fG += pfdx[X3(ix-1,iy,iz)] - pfdx[iX] - pfbx[X3(ix-1,iy,iz)] + pfbx[iX];
                    fG += pfdy[X3(ix,iy-1,iz)] - pfdy[iX] - pfby[X3(ix,iy-1,iz)] + pfby[iX];
                    fG += pfdz[X3(ix,iy,iz-1)] - pfdz[iX] - pfbz[X3(ix,iy,iz-1)] + pfbz[iX];
                    fG *= fct1b;
                    fG -= fct2b* pfHr[iX];
                    if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                    u[iX] = fG;
                }
            
            /* edges */ // *12 edges - each 4 neighbours*
            ix=0; iy=0;
            for (iz=1; iz< iNz-1; iz++)
            {
                iX = X3(ix,iy,iz);
                
                fG = u[X3(ix+1,iy,iz)] + u[X3(ix,iy+1,iz)] + u[X3(ix,iy,iz+1)] + u[X3(ix,iy,iz-1)];
                fG += - pfdx[iX] + pfbx[iX];
                fG += - pfdy[iX] + pfby[iX];
                fG += pfdz[X3(ix,iy,iz-1)] - pfdz[iX] - pfbz[X3(ix,iy,iz-1)] + pfbz[iX];
                fG *= fct1e;
                fG -= fct2e* pfHr[iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                u[iX] = fG;
            }
            
            ix=iNx-1; iy=0;
            for (iz=1; iz< iNz-1; iz++)
            {
                iX = X3(ix,iy,iz);
                
                fG = u[X3(ix-1,iy,iz)] + u[X3(ix,iy+1,iz)] + u[X3(ix,iy,iz+1)] + u[X3(ix,iy,iz-1)];
                fG += pfdx[X3(ix-1,iy,iz)] - pfdx[iX] - pfbx[X3(ix-1,iy,iz)] + pfbx[iX];
                fG += - pfdy[iX] + pfby[iX];
                fG += pfdz[X3(ix,iy,iz-1)] - pfdz[iX] - pfbz[X3(ix,iy,iz-1)] + pfbz[iX];
                fG *= fct1e;
                fG -= fct2e* pfHr[iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                u[iX] = fG;
            }
            
            ix=0; iy=iNy-1;
            for (iz=1; iz< iNz-1; iz++)
            {
                iX = X3(ix,iy,iz);
                
                fG = u[X3(ix+1,iy,iz)] + u[X3(ix,iy-1,iz)] + u[X3(ix,iy,iz+1)] + u[X3(ix,iy,iz-1)];
                fG += - pfdx[iX] + pfbx[iX];
                fG += pfdy[X3(ix,iy-1,iz)] - pfdy[iX] - pfby[X3(ix,iy-1,iz)] + pfby[iX];
                fG += pfdz[X3(ix,iy,iz-1)] - pfdz[iX] - pfbz[X3(ix,iy,iz-1)] + pfbz[iX];
                fG *= fct1e;
                fG -= fct2e* pfHr[iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                u[iX] = fG;
            }
            
            ix=iNx-1; iy=iNy-1;
            for (iz=1; iz< iNz-1; iz++)
            {
                iX = X3(ix,iy,iz);
                
                fG = u[X3(ix-1,iy,iz)] + u[X3(ix,iy-1,iz)] + u[X3(ix,iy,iz+1)] + u[X3(ix,iy,iz-1)];
                fG += pfdx[X3(ix-1,iy,iz)] - pfdx[iX] - pfbx[X3(ix-1,iy,iz)] + pfbx[iX];
                fG += pfdy[X3(ix,iy-1,iz)] - pfdy[iX] - pfby[X3(ix,iy-1,iz)] + pfby[iX];
                fG += pfdz[X3(ix,iy,iz-1)] - pfdz[iX] - pfbz[X3(ix,iy,iz-1)] + pfbz[iX];
                fG *= fct1e;
                fG -= fct2e* pfHr[iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                u[iX] = fG;
            }
            
            iy=0; iz=0;
            for (ix=1; ix< iNx-1; ix++)
            {
                iX = X3(ix,iy,iz);
                
                fG = u[X3(ix-1,iy,iz)] + u[X3(ix+1,iy,iz)] + u[X3(ix,iy+1,iz)] + u[X3(ix,iy,iz+1)];
                fG += pfdx[X3(ix-1,iy,iz)] - pfdx[iX] - pfbx[X3(ix-1,iy,iz)] + pfbx[iX];
                fG += - pfdy[iX] + pfby[iX];
                fG += - pfdz[iX] + pfbz[iX];
                fG *= fct1e;
                fG -= fct2e* pfHr[iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                u[iX] = fG;
            }
            
            iy=0; iz=iNz-1;
            for (ix=1; ix< iNx-1; ix++)
            {
                iX = X3(ix,iy,iz);
                
                fG = u[X3(ix-1,iy,iz)] + u[X3(ix+1,iy,iz)] + u[X3(ix,iy+1,iz)] + u[X3(ix,iy,iz-1)];
                fG += pfdx[X3(ix-1,iy,iz)] - pfdx[iX] - pfbx[X3(ix-1,iy,iz)] + pfbx[iX];
                fG += - pfdy[iX] + pfby[iX];
                fG += pfdz[X3(ix,iy,iz-1)] - pfdz[iX] - pfbz[X3(ix,iy,iz-1)] + pfbz[iX];
                fG *= fct1e;
                fG -= fct2e* pfHr[iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                u[iX] = fG;
            }
            
            iy=iNy-1; iz=0;
            for (ix=1; ix< iNx-1; ix++)
            {
                iX = X3(ix,iy,iz);
                
                fG = u[X3(ix-1,iy,iz)] + u[X3(ix+1,iy,iz)] + u[X3(ix,iy-1,iz)] + u[X3(ix,iy,iz+1)];
                fG += pfdx[X3(ix-1,iy,iz)] - pfdx[iX] - pfbx[X3(ix-1,iy,iz)] + pfbx[iX];
                fG += pfdy[X3(ix,iy-1,iz)] - pfdy[iX] - pfby[X3(ix,iy-1,iz)] + pfby[iX];
                fG += - pfdz[iX] + pfbz[iX];
                fG *= fct1e;
                fG -= fct2e* pfHr[iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                u[iX] = fG;
            }
            
            iy=iNy-1; iz=iNz-1;
            for (ix=1; ix< iNx-1; ix++)
            {
                iX = X3(ix,iy,iz);
                
                fG = u[X3(ix-1,iy,iz)] + u[X3(ix+1,iy,iz)] + u[X3(ix,iy-1,iz)] + u[X3(ix,iy,iz-1)];
                fG += pfdx[X3(ix-1,iy,iz)] - pfdx[iX] - pfbx[X3(ix-1,iy,iz)] + pfbx[iX];
                fG += pfdy[X3(ix,iy-1,iz)] - pfdy[iX] - pfby[X3(ix,iy-1,iz)] + pfby[iX];
                fG += pfdz[X3(ix,iy,iz-1)] - pfdz[iX] - pfbz[X3(ix,iy,iz-1)] + pfbz[iX];
                fG *= fct1e;
                fG -= fct2e* pfHr[iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                u[iX] = fG;
            }
            
            ix=0; iz=0;
            for (iy=1; iy< iNy-1; iy++)
            {
                iX = X3(ix,iy,iz);
                
                fG = u[X3(ix+1,iy,iz)] + u[X3(ix,iy-1,iz)] + u[X3(ix,iy+1,iz)] + u[X3(ix,iy,iz+1)];
                fG += - pfdx[iX] + pfbx[iX];
                fG += pfdy[X3(ix,iy-1,iz)] - pfdy[iX] - pfby[X3(ix,iy-1,iz)] + pfby[iX];
                fG += - pfdz[iX] + pfbz[iX];
                fG *= fct1e;
                fG -= fct2e* pfHr[iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                u[iX] = fG;
            }
            
            ix=iNx-1; iz=0;
            for (iy=1; iy< iNy-1; iy++)
            {
                iX = X3(ix,iy,iz);
                
                fG = u[X3(ix-1,iy,iz)] + u[X3(ix,iy-1,iz)] + u[X3(ix,iy+1,iz)] + u[X3(ix,iy,iz+1)];
                fG += pfdx[X3(ix-1,iy,iz)] - pfdx[iX] - pfbx[X3(ix-1,iy,iz)] + pfbx[iX];
                fG += pfdy[X3(ix,iy-1,iz)] - pfdy[iX] - pfby[X3(ix,iy-1,iz)] + pfby[iX];
                fG += - pfdz[iX] + pfbz[iX];
                fG *= fct1e;
                fG -= fct2e* pfHr[iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                u[iX] = fG;
            }
            
            ix=0; iz=iNz-1;
            for (iy=1; iy< iNy-1; iy++)
            {
                iX = X3(ix,iy,iz);
                
                fG = u[X3(ix+1,iy,iz)] + u[X3(ix,iy-1,iz)] + u[X3(ix,iy+1,iz)] + u[X3(ix,iy,iz-1)];
                fG += - pfdx[iX] + pfbx[iX];
                fG += pfdy[X3(ix,iy-1,iz)] - pfdy[iX] - pfby[X3(ix,iy-1,iz)] + pfby[iX];
                fG += pfdz[X3(ix,iy,iz-1)] - pfdz[iX] - pfbz[X3(ix,iy,iz-1)] + pfbz[iX];
                fG *= fct1e;
                fG -= fct2e* pfHr[iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                u[iX] = fG;
            }
            
            ix=iNx-1; iz=iNz-1;
            for (iy=1; iy< iNy-1; iy++)
            {
                iX = X3(ix,iy,iz);
                
                fG = u[X3(ix-1,iy,iz)] + u[X3(ix,iy-1,iz)] + u[X3(ix,iy+1,iz)] + u[X3(ix,iy,iz-1)];
                fG += pfdx[X3(ix-1,iy,iz)] - pfdx[iX] - pfbx[X3(ix-1,iy,iz)] + pfbx[iX];
                fG += pfdy[X3(ix,iy-1,iz)] - pfdy[iX] - pfby[X3(ix,iy-1,iz)] + pfby[iX];
                fG += pfdz[X3(ix,iy,iz-1)] - pfdz[iX] - pfbz[X3(ix,iy,iz-1)] + pfbz[iX];
                fG *= fct1e;
                fG -= fct2e* pfHr[iX];
                if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
                u[iX] = fG;
            }
            
            /* Corners */ // *8 corners - each 3 neighbours*
            ix=0; iy=0; iz=0;
            iX = X3(ix,iy,iz);
            fG = u[X3(ix+1,iy,iz)] + u[X3(ix,iy+1,iz)] + u[X3(ix,iy,iz+1)];
            fG += - pfdx[iX] + pfbx[iX];
            fG += - pfdy[iX] + pfby[iX];
            fG += - pfdz[iX] + pfbz[iX];
            fG *= fct1c;
            fG -= fct2c* pfHr[iX];
            if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
            u[iX] = fG;
            
            ix=iNx-1; iy=0; iz=0;
            iX = X3(ix,iy,iz);
            fG = u[X3(ix-1,iy,iz)] + u[X3(ix,iy+1,iz)] + u[X3(ix,iy,iz+1)];
            fG += pfdx[X3(ix-1,iy,iz)] - pfdx[iX] - pfbx[X3(ix-1,iy,iz)] + pfbx[iX];
            fG += - pfdy[iX] + pfby[iX];
            fG += - pfdz[iX] + pfbz[iX];
            fG *= fct1c;
            fG -= fct2c* pfHr[iX];
            if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
            u[iX] = fG;
            
            ix=0; iy=iNy-1; iz=0;
            iX = X3(ix, iy, iz);
            fG = u[X3(ix+1, iy, iz)] + u[X3(ix, iy-1, iz)] + u[X3(ix, iy, iz+1)];
            fG += - pfdx[iX] + pfbx[iX];
            fG += pfdy[X3(ix,iy-1,iz)] - pfdy[iX] - pfby[X3(ix, iy-1,iz)] + pfby[iX];
            fG += - pfdz[iX] + pfbz[iX];
            fG *= fct1c;
            fG -= fct2c* pfHr[iX];
            if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
            u[iX] = fG;
            
            ix=iNx-1; iy=iNy-1; iz=0;
            iX = X3(ix,iy,iz);
            fG = u[X3(ix-1,iy,iz)] + u[X3(ix,iy-1,iz)]+ u[X3(ix, iy, iz+1)];
            fG += pfdx[X3(ix-1,iy,iz)] - pfdx[iX] - pfbx[X3(ix-1,iy,iz)] + pfbx[iX];
            fG += pfdy[X3(ix,iy-1,iz)] - pfdy[iX] - pfby[X3(ix,iy-1,iz)] + pfby[iX];
            fG += - pfdz[iX] + pfbz[iX];
            fG *= fct1c;
            fG -= fct2c* pfHr[iX];
            if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
            u[iX] = fG;
            
            ix=0; iy=0; iz=iNz-1;
            iX = X3(ix,iy,iz);
            fG = u[X3(ix+1,iy,iz)] + u[X3(ix,iy+1,iz)] + u[X3(ix,iy,iz-1)];
            fG += - pfdx[iX] + pfbx[iX];
            fG += - pfdy[iX] + pfby[iX];
            fG += pfdz[X3(ix,iy,iz-1)] - pfdz[iX] - pfbz[X3(ix,iy,iz-1)] + pfbz[iX];
            fG *= fct1c;
            fG -= fct2c* pfHr[iX];
            if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
            u[iX] = fG;
            
            ix=iNx-1; iy=0; iz=iNz-1;
            iX = X3(ix,iy,iz);
            fG = u[X3(ix-1,iy,iz)] + u[X3(ix,iy+1,iz)] + u[X3(ix,iy,iz-1)];
            fG += pfdx[X3(ix-1,iy,iz)] - pfdx[iX] - pfbx[X3(ix-1,iy,iz)] + pfbx[iX];
            fG += - pfdy[iX] + pfby[iX];
            fG += pfdz[X3(ix,iy,iz-1)] - pfdz[iX] - pfbz[X3(ix,iy,iz-1)] + pfbz[iX];
            fG *= fct1c;
            fG -= fct2c* pfHr[iX];
            if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
            u[iX] = fG;
            
            ix=0; iy=iNy-1; iz=iNz-1;
            iX = X3(ix, iy, iz);
            fG = u[X3(ix+1, iy, iz)] + u[X3(ix, iy-1, iz)] + u[X3(ix, iy, iz-1)];
            fG += - pfdx[iX] + pfbx[iX];
            fG += pfdy[X3(ix,iy-1,iz)] - pfdy[iX] - pfby[X3(ix, iy-1,iz)] + pfby[iX];
            fG += pfdz[X3(ix,iy,iz-1)] - pfdz[iX] - pfbz[X3(ix,iy,iz-1)] + pfbz[iX];
            fG *= fct1c;
            fG -= fct2c* pfHr[iX];
            if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
            u[iX] = fG;
            
            ix=iNx-1; iy=iNy-1; iz=iNz-1;
            iX = X3(ix,iy,iz);
            fG = u[X3(ix-1,iy,iz)] + u[X3(ix,iy-1,iz)]+ u[X3(ix, iy, iz-1)];
            fG += pfdx[X3(ix-1,iy,iz)] - pfdx[iX] - pfbx[X3(ix-1,iy,iz)] + pfbx[iX];
            fG += pfdy[X3(ix,iy-1,iz)] - pfdy[iX] - pfby[X3(ix,iy-1,iz)] + pfby[iX];
            fG += pfdz[X3(ix,iy,iz-1)] - pfdz[iX] - pfbz[X3(ix,iy,iz-1)] + pfbz[iX];
            fG *= fct1c;
            fG -= fct2c* pfHr[iX];
            if(fG>1.0) fG=1.0; else if(fG<0.0) fG=0.0;
            u[iX] = fG;
            /* end Borders */
            
            
            /* Compute diff ( u - uold ) */
            fSumDiff = 0.0;
            fSumU = 0.0;
            fSumUold = 0.0;
            for (iz=0; iz< iNz; iz++)
                for (iy=0; iy< iNy; iy++)
                    for (ix=0; ix< iNx; ix++)
                        fSumDiff += SQR(u[X3(ix,iy,iz)]-uOld2[X3(ix,iy,iz)]);
            fDiffNew2 = fSumDiff/ fNyxz;
            if ( iGS==0 )
            {
                fDiffFirst2 = fDiffNew2;
                fDiffNew2 = 1e10;
                fError = 1e10;
            }
            else
                fError = 1.0 - std::abs(fDiffNew2-fDiffFirst2)/fDiffFirst2;
            iGS++;
        } //End gauss seidel
        
        iMeanGS += iGS;
        iCptGS++;
        
        /* Center */
        for (iz=0; iz< iNz-1; iz++)
            for (iy=0; iy< iNy-1; iy++)
                for (ix=0; ix< iNx-1; ix++)
                {
                    iX = X3(ix,iy,iz);
                    /* d */
                    fDxu = u[X3(ix+1,iy,iz)] - u[iX];
                    fDyu = u[X3(ix,iy+1,iz)] - u[iX];
                    fDzu = u[X3(ix,iy,iz+1)] - u[iX];
                    f1 = fDxu+pfbx[iX];
                    f2 = fDyu+pfby[iX];
                    f3 = fDzu+pfbz[iX];
                    fs = SQR(f1)+SQR(f2)+SQR(f3);
                    fctST = fInvMu2* pfGb2[iX];
                    if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; pfdz[iX]=0.0; }
                    else {
                        fs = std::sqrt(fs);
                        fctST = std::sqrt(fctST);
                        fTemp = fs-fctST; fTemp /= fs;
                        pfdx[iX] = fTemp* f1;
                        pfdy[iX] = fTemp* f2;
                        pfdz[iX] = fTemp* f3;}
                    /* b */
                    pfbx[iX] += fDxu - pfdx[iX];
                    pfby[iX] += fDyu - pfdy[iX];
                    pfbz[iX] += fDzu - pfdz[iX];
                }
        
        /* Borders */
        ix=iNx-1;
        for (iz=1; iz< iNz-1; iz++)
            for (iy=1; iy< iNy-1; iy++)
            {
                iX = X3(ix,iy,iz);
                fDxu = 0.0;
                fDyu = u[X3(ix,iy+1,iz)] - u[iX];
                fDzu = u[X3(ix,iy,iz+1)] - u[iX];
                f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX]; f3 = fDzu+pfbz[iX];
                fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
                if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; pfdz[iX]=0.0; }
                else {
                    fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                    fTemp = fs-fctST; fTemp /= fs;
                    pfdx[iX] = fTemp* f1;
                    pfdy[iX] = fTemp* f2;
                    pfdz[iX] = fTemp* f3;}
                pfbx[iX] += fDxu - pfdx[iX];
                pfby[iX] += fDyu - pfdy[iX];
                pfbz[iX] += fDzu - pfdz[iX];
            }
        
        iy=iNy-1;
        for (iz=1; iz< iNz-1; iz++)
            for (ix=1; ix< iNx-1; ix++)
            {
                iX = X3(ix,iy,iz);
                fDxu = u[X3(ix+1,iy,iz)] - u[iX];
                fDyu = 0.0;
                fDzu = u[X3(ix,iy,iz+1)] - u[iX];
                f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX]; f3 = fDzu+pfbz[iX];
                fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
                if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; pfdz[iX]=0.0; }
                else {
                    fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                    fTemp = fs-fctST; fTemp /= fs;
                    pfdx[iX] = fTemp* f1;
                    pfdy[iX] = fTemp* f2;
                    pfdz[iX] = fTemp* f3;}
                pfbx[iX] += fDxu - pfdx[iX];
                pfby[iX] += fDyu - pfdy[iX];
                pfbz[iX] += fDzu - pfdz[iX];
            }
        
        iz=iNz-1;
        for (iy=1; iy< iNy-1; iy++)
            for (ix=1; ix< iNx-1; ix++)
            {
                iX = X3(ix,iy,iz);
                fDxu = u[X3(ix+1,iy,iz)] - u[iX];
                fDyu = u[X3(ix,iy+1,iz)] - u[iX];
                fDzu = 0.0;
                f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX]; f3 = fDzu+pfbz[iX];
                fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
                if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; pfdz[iX]=0.0; }
                else {
                    fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                    fTemp = fs-fctST; fTemp /= fs;
                    pfdx[iX] = fTemp* f1;
                    pfdy[iX] = fTemp* f2;
                    pfdz[iX] = fTemp* f3;}
                pfbx[iX] += fDxu - pfdx[iX];
                pfby[iX] += fDyu - pfdy[iX];
                pfbz[iX] += fDzu - pfdz[iX];
            }
        
        /*edges*/
        ix=iNx-1; iz=0;
        for (iy=1; iy< iNy-1; iy++)
        {
            iX = X3(ix,iy,iz);
            fDxu = 0.0;
            fDyu = u[X3(ix,iy+1,iz)] - u[iX];
            fDzu = u[X3(ix,iy,iz+1)] - u[iX];
            f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX]; f3 = fDzu+pfbz[iX];
            fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
            if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; pfdz[iX]=0.0; }
            else {
                fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                fTemp = fs-fctST; fTemp /= fs;
                pfdx[iX] = fTemp* f1;
                pfdy[iX] = fTemp* f2;
                pfdz[iX] = fTemp* f3;}
            pfbx[iX] += fDxu - pfdx[iX];
            pfby[iX] += fDyu - pfdy[iX];
            pfbz[iX] += fDzu - pfdz[iX];
        }
        
        ix=0; iz=iNz-1;
        for (iy=1; iy< iNy-1; iy++)
        {
            iX = X3(ix,iy,iz);
            fDxu = u[X3(ix+1,iy,iz)] - u[iX];
            fDyu = u[X3(ix,iy+1,iz)] - u[iX];
            fDzu = 0.0;
            f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX]; f3 = fDzu+pfbz[iX];
            fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
            if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; pfdz[iX]=0.0; }
            else {
                fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                fTemp = fs-fctST; fTemp /= fs;
                pfdx[iX] = fTemp* f1;
                pfdy[iX] = fTemp* f2;
                pfdz[iX] = fTemp* f3;}
            pfbx[iX] += fDxu - pfdx[iX];
            pfby[iX] += fDyu - pfdy[iX];
            pfbz[iX] += fDzu - pfdz[iX];
        }
        
        ix=iNx-1; iz=iNz-1;
        for (iy=1; iy< iNy-1; iy++)
        {
            iX = X3(ix,iy,iz);
            fDxu = 0.0;
            fDyu = u[X3(ix,iy+1,iz)] - u[iX];
            fDzu = 0.0;
            f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX]; f3 = fDzu+pfbz[iX];
            fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
            if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; pfdz[iX]=0.0; }
            else {
                fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                fTemp = fs-fctST; fTemp /= fs;
                pfdx[iX] = fTemp* f1;
                pfdy[iX] = fTemp* f2;
                pfdz[iX] = fTemp* f3;}
            pfbx[iX] += fDxu - pfdx[iX];
            pfby[iX] += fDyu - pfdy[iX];
            pfbz[iX] += fDzu - pfdz[iX];
        }
        
        
        iy=0; iz=iNz-1;
        for (ix=1; ix< iNx-1; ix++)
        {
            iX = X3(ix,iy,iz);
            fDxu = u[X3(ix+1,iy,iz)] - u[iX];
            fDyu = u[X3(ix,iy+1,iz)] - u[iX];
            fDzu = 0.0;
            f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX]; f3 = fDzu+pfbz[iX];
            fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
            if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; pfdz[iX]=0.0; }
            else {
                fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                fTemp = fs-fctST; fTemp /= fs;
                pfdx[iX] = fTemp* f1;
                pfdy[iX] = fTemp* f2;
                pfdz[iX] = fTemp* f3;}
            pfbx[iX] += fDxu - pfdx[iX];
            pfby[iX] += fDyu - pfdy[iX];
            pfbz[iX] += fDzu - pfdz[iX];
        }
        
        iy=iNy-1; iz=0;
        for (ix=1; ix< iNx-1; ix++)
        {
            iX = X3(ix,iy,iz);
            fDxu = u[X3(ix+1,iy,iz)] - u[iX];
            fDyu = 0.0;
            fDzu = u[X3(ix,iy,iz+1)] - u[iX];
            f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX]; f3 = fDzu+pfbz[iX];
            fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
            if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; pfdz[iX]=0.0; }
            else {
                fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                fTemp = fs-fctST; fTemp /= fs;
                pfdx[iX] = fTemp* f1;
                pfdy[iX] = fTemp* f2;
                pfdz[iX] = fTemp* f3;}
            pfbx[iX] += fDxu - pfdx[iX];
            pfby[iX] += fDyu - pfdy[iX];
            pfbz[iX] += fDzu - pfdz[iX];
        }
        
        iy=iNy-1; iz=iNz-1;
        for (ix=1; ix< iNx-1; ix++)
        {
            iX = X3(ix,iy,iz);
            fDxu = u[X3(ix+1,iy,iz)] - u[iX];
            fDyu = 0.0;
            fDzu = 0.0;
            f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX]; f3 = fDzu+pfbz[iX];
            fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
            if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; pfdz[iX]=0.0; }
            else {
                fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                fTemp = fs-fctST; fTemp /= fs;
                pfdx[iX] = fTemp* f1;
                pfdy[iX] = fTemp* f2;
                pfdz[iX] = fTemp* f3;}
            pfbx[iX] += fDxu - pfdx[iX];
            pfby[iX] += fDyu - pfdy[iX];
            pfbz[iX] += fDzu - pfdz[iX];
        }
        
        ix=iNx-1; iy=0;
        for (iz=1; iz< iNz-1; iz++)
        {
            iX = X3(ix,iy,iz);
            fDxu = 0.0;
            fDyu = u[X3(ix,iy+1,iz)] - u[iX];
            fDzu = u[X3(ix,iy,iz+1)] - u[iX];
            f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX]; f3 = fDzu+pfbz[iX];
            fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
            if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; pfdz[iX]=0.0; }
            else {
                fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                fTemp = fs-fctST; fTemp /= fs;
                pfdx[iX] = fTemp* f1;
                pfdy[iX] = fTemp* f2;
                pfdz[iX] = fTemp* f3;}
            pfbx[iX] += fDxu - pfdx[iX];
            pfby[iX] += fDyu - pfdy[iX];
            pfbz[iX] += fDzu - pfdz[iX];
        }
        
        ix=0; iy=iNy-1;
        for (iz=1; iz< iNz-1; iz++)
        {
            iX = X3(ix,iy,iz);
            fDxu = u[X3(ix+1,iy,iz)] - u[iX];
            fDyu = 0.0;
            fDzu = u[X3(ix,iy,iz+1)] - u[iX];
            f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX]; f3 = fDzu+pfbz[iX];
            fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
            if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; pfdz[iX]=0.0; }
            else {
                fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                fTemp = fs-fctST; fTemp /= fs;
                pfdx[iX] = fTemp* f1;
                pfdy[iX] = fTemp* f2;
                pfdz[iX] = fTemp* f3;}
            pfbx[iX] += fDxu - pfdx[iX];
            pfby[iX] += fDyu - pfdy[iX];
            pfbz[iX] += fDzu - pfdz[iX];
        }
        
        ix=iNx-1; iy=iNy-1;
        for (iz=1; iz< iNz-1; iz++)
        {
            iX = X3(ix,iy,iz);
            fDxu = 0.0;
            fDyu = 0.0;
            fDzu = u[X3(ix,iy,iz+1)] - u[iX];
            f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX]; f3 = fDzu+pfbz[iX];
            fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
            if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; pfdz[iX]=0.0; }
            else {
                fs = std::sqrt(fs); fctST = std::sqrt(fctST);
                fTemp = fs-fctST; fTemp /= fs;
                pfdx[iX] = fTemp* f1;
                pfdy[iX] = fTemp* f2;
                pfdz[iX] = fTemp* f3;}
            pfbx[iX] += fDxu - pfdx[iX];
            pfby[iX] += fDyu - pfdy[iX];
            pfbz[iX] += fDzu - pfdz[iX];
        }
        
        /*corners*/
        ix=iNx-1; iy=0; iz=0;
        iX = X3(ix,iy,iz);
        fDxu = 0.0;
        fDyu = u[X3(ix,iy+1,iz)] - u[iX];
        fDzu = u[X3(ix,iy,iz+1)] - u[iX];
        f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX]; f3 = fDzu+pfbz[iX];
        fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
        if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; pfdz[iX]=0.0; }
        else {
            fs = std::sqrt(fs); fctST = std::sqrt(fctST);
            fTemp = fs-fctST; fTemp /= fs;
            pfdx[iX] = fTemp* f1;
            pfdy[iX] = fTemp* f2;
            pfdz[iX] = fTemp* f3;}
        pfbx[iX] += fDxu - pfdx[iX];
        pfby[iX] += fDyu - pfdy[iX];
        pfbz[iX] += fDzu - pfdz[iX];
        
        ix=0; iy=iNy-1; iz=0;
        iX = X3(ix,iy,iz);
        fDxu = u[X3(ix+1,iy,iz)] - u[iX];
        fDyu = 0.0;
        fDzu = u[X3(ix,iy,iz+1)] - u[iX];
        f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX]; f3 = fDzu+pfbz[iX];
        fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
        if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; pfdz[iX]=0.0; }
        else {
            fs = std::sqrt(fs); fctST = std::sqrt(fctST);
            fTemp = fs-fctST; fTemp /= fs;
            pfdx[iX] = fTemp* f1;
            pfdy[iX] = fTemp* f2;
            pfdz[iX] = fTemp* f3;}
        pfbx[iX] += fDxu - pfdx[iX];
        pfby[iX] += fDyu - pfdy[iX];
        pfbz[iX] += fDzu - pfdz[iX];
        
        ix=iNx-1; iy=iNy-1; iz=0;
        iX = X3(ix,iy,iz);
        fDxu = 0.0;
        fDyu = 0.0;
        fDzu = u[X3(ix,iy,iz+1)] - u[iX];
        f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX]; f3 = fDzu+pfbz[iX];
        fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
        if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; pfdz[iX]=0.0;}
        else {
            fs = std::sqrt(fs); fctST = std::sqrt(fctST);
            fTemp = fs-fctST; fTemp /= fs;
            pfdx[iX] = fTemp* f1;
            pfdy[iX] = fTemp* f2;
            pfdz[iX] = fTemp* f3;}
        pfbx[iX] += fDxu - pfdx[iX];
        pfby[iX] += fDyu - pfdy[iX];
        pfbz[iX] += fDzu - pfdz[iX];
        
        ix=iNx-1; iy=0; iz=iNz-1;
        iX = X3(ix,iy,iz);
        fDxu = 0.0;
        fDyu = u[X3(ix,iy+1,iz)] - u[iX];
        fDzu = 0.0;
        f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX]; f3 = fDzu+pfbz[iX];
        fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
        if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; pfdz[iX]=0.0; }
        else {
            fs = std::sqrt(fs); fctST = std::sqrt(fctST);
            fTemp = fs-fctST; fTemp /= fs;
            pfdx[iX] = fTemp* f1;
            pfdy[iX] = fTemp* f2;
            pfdz[iX] = fTemp* f3;}
        pfbx[iX] += fDxu - pfdx[iX];
        pfby[iX] += fDyu - pfdy[iX];
        pfbz[iX] += fDzu - pfdz[iX];
        
        ix=0; iy=iNy-1; iz=iNz-1;
        iX = X3(ix,iy,iz);
        fDxu = u[X3(ix+1,iy,iz)] - u[iX];
        fDyu = 0.0;
        fDzu = 0.0;
        f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX]; f3 = fDzu+pfbz[iX];
        fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
        if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; pfdz[iX]=0.0; }
        else {
            fs = std::sqrt(fs); fctST = std::sqrt(fctST);
            fTemp = fs-fctST; fTemp /= fs;
            pfdx[iX] = fTemp* f1;
            pfdy[iX] = fTemp* f2;
            pfdz[iX] = fTemp* f3;}
        pfbx[iX] += fDxu - pfdx[iX];
        pfby[iX] += fDyu - pfdy[iX];
        pfbz[iX] += fDzu - pfdz[iX];
        
        ix=iNx-1; iy=iNy-1; iz=iNz-1;
        iX = X3(ix,iy,iz);
        fDxu = 0.0;
        fDyu = 0.0;
        fDzu = 0.0;
        f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX]; f3 = fDzu+pfbz[iX];
        fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
        if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; pfdz[iX]=0.0;}
        else {
            fs = std::sqrt(fs); fctST = std::sqrt(fctST);
            fTemp = fs-fctST; fTemp /= fs;
            pfdx[iX] = fTemp* f1;
            pfdy[iX] = fTemp* f2;
            pfdz[iX] = fTemp* f3;}
        pfbx[iX] += fDxu - pfdx[iX];
        pfby[iX] += fDyu - pfdy[iX];
        pfbz[iX] += fDzu - pfdz[iX];
        
        ix=0; iy=0; iz=iNz-1;
        iX = X3(ix,iy,iz);
        fDxu = u[X3(ix+1,iy,iz)] - u[iX];
        fDyu = u[X3(ix,iy+1,iz)] - u[iX];
        fDzu = 0.0;
        f1 = fDxu+pfbx[iX]; f2 = fDyu+pfby[iX]; f3 = fDzu+pfbz[iX];
        fs = SQR(f1)+SQR(f2)+SQR(f3); fctST = fInvMu2* pfGb2[iX];
        if ( fs<fctST ) { pfdx[iX]=0.0; pfdy[iX]=0.0; pfdz[iX]=0.0;}
        else {
            fs = std::sqrt(fs); fctST = std::sqrt(fctST);
            fTemp = fs-fctST; fTemp /= fs;
            pfdx[iX] = fTemp* f1;
            pfdy[iX] = fTemp* f2;
            pfdz[iX] = fTemp* f3;}
        pfbx[iX] += fDxu - pfdx[iX];
        pfby[iX] += fDyu - pfdy[iX];
        pfbz[iX] += fDzu - pfdz[iX];
        /* Borders */
        
        
        // CHECK CONVERGENCE
         pDeltaSumU = 0.0;
         pSumUold = 0.0;
         pSumU = 0.0;

         for (int iz=0; iz< iNz; iz++)
             for (int iy=0; iy< iNy; iy++)
                 for (int ix=0; ix< iNx; ix++)
                 {
                     pDeltaSumU += SQR(u[X3(ix,iy,iz)]-uOld[X3(ix,iy,iz)]);
                     pSumU += SQR(u[X3(ix,iy,iz)]);
                     pSumUold += SQR(uOld[X3(ix,iy,iz)]);
                 }

         pDeltaSumU = pDeltaSumU/ (pSumU*pSumUold);


         // SCREEN FEEDBACK
         std::stringstream ss; ss << k;
         std::cout << " End SEG iteration: " << ss.str() << " | Relative change u: " << pDeltaSumU << std::endl;


         if(std::abs(pDeltaSumU)<pEps) break;
         k++;
        
    }
    
    
    // UPDATE ALPHA
    if(bUpdateAlpha){
        
        pixel_type sum1 = 0.0;
        pixel_type sum2 = 0.0;
        for (int iz=0; iz< iNz; iz++)
        for (int iy=0; iy< iNy; iy++)
        for (int ix=0; ix< iNx; ix++){
            if(brainMask[X3(ix,iy,iz)]){
                sum1 += u[X3(ix,iy,iz)];
                sum2 += 1.0;
            }
        }
        std::fill(alpha.begin(),alpha.end(),sum1/sum2);
    }
    
    return;

}

}
}
#endif //SEG_JOINT_HPP_INCLUDED
