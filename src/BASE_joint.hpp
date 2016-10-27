#ifndef BASE_JOINT_HPP_INCLUDED
#define BASE_JOINT_HPP_INCLUDED

#include "EM_joint.hpp"
#include "SEG_joint.hpp"

#include "_help.hpp"

#define X3(ix,iy,iz) ((iz)*iNy+(iy))*iNx+(ix)

namespace image
{
namespace segmentation
{
    
template<typename pixel_type,unsigned int dimension>
void BASE_joint(std::vector<basic_image<pixel_type,dimension> >& im,
                    basic_image<pixel_type,dimension>& brainMask,
                    basic_image<pixel_type,dimension>& init,
                    basic_image<pixel_type,dimension>& u,
          
                    std::vector<basic_image<pixel_type,dimension> >& priorIn,
                    std::vector<basic_image<pixel_type,dimension> >& priorOut,
          
                    std::vector<basic_image<pixel_type,dimension> >& posteriorIn,
                    std::vector<basic_image<pixel_type,dimension> >& posteriorOut,
          
                    std::vector<std::vector<pixel_type> >& pMeanIn,
                    std::vector<std::vector<pixel_type> >& pSigmaIn,
                    std::vector<std::vector<pixel_type> >& pMeanOut,
                    std::vector<std::vector<pixel_type> >& pSigmaOut,
                    std::vector<pixel_type>& pPriorIn,
          
                    int iNrOfOtherClasses = 0,
          
                    pixel_type pStdMin = 1,
                    pixel_type pMin = 1e-40,
                    pixel_type pIOffset = 1.0,
                    pixel_type pIMax = 256.0,
          
                    int iMaxIter = 20,
                    int iMaxIterEM = 1,
                    int iMaxIterSeg = 1,
          
                    pixel_type pEpsEM = 1e-3,
                    pixel_type pEpsSeg = 1e-9,
             
                    pixel_type pLambda = 1e1,
                    pixel_type pMu = 1e1,
             
                    pixel_type pInitAlpha = 0.5,
                    bool bUpdateAlpha = true,
             
                    pixel_type pCoupling = 1.0,
             
                    bool bEMreinitialize = true,
             
                    std::string basename = "",
                    bool verbose = false,
             
                    bool bInit = false,
                
                    std::vector<int> iHyperHypoIntens_toWhich = std::vector<int>(),
                    std::vector<pixel_type> pHyperHypoIntens_howMuch = std::vector<pixel_type>(),
                    int iInclusion_toWhich = -1){
    
    
    
    // NORMALIZE
    normalizeV(im,brainMask,pIMax);

    
    // CREATE OTHER CLASSES
    create_other_classes(priorOut,brainMask,iNrOfOtherClasses);
    

    // DECLARE INT
    int iNrOfChannels = im.size();
    
    int iNrOfInClasses  = priorIn.size();
    int iNrOfOutClasses = priorOut.size();
    
    geometry<dimension> geo = im[0].geometry();
    int iNx = geo[0];
    int iNy = geo[1];
    int iNz = geo[2];
    
    
    // DECLARE PIXEL_TYPE
    pixel_type pSumU,pSumUold,pDeltaSumU;
    
    
    // DECLARE BASIC_IMAGE
    basic_image<pixel_type,dimension> alpha(geo);
    basic_image<pixel_type,dimension> uOld(geo);
    basic_image<pixel_type,dimension> totalIn(geo);
    basic_image<pixel_type,dimension> totalOut(geo);
    
    basic_image<pixel_type,dimension> pfHr(geo);
  
    
    // POSSIBLE INITIALIZATION OF U
    for (int iz=0; iz< iNz; iz++)
        for (int iy=0; iy< iNy; iy++)
            for (int ix=0; ix< iNx; ix++){
                
                if(bInit){
                    if(init[X3(ix,iy,iz)]>0.9) init[X3(ix,iy,iz)] = 0.9;
                    if(init[X3(ix,iy,iz)]<0.1) init[X3(ix,iy,iz)] = 0.1;
                    u[X3(ix,iy,iz)] = init[X3(ix,iy,iz)];
                }
                else{
                    if(brainMask[X3(ix,iy,iz)]) u[X3(ix,iy,iz)] = pInitAlpha;
                    else u[X3(ix,iy,iz)] = 0.0;
                }
            }
    
    
    // POSSIBLE INITIALIZATION OF ALPHA
    for (int iz=0; iz< iNz; iz++)
        for (int iy=0; iy< iNy; iy++)
            for (int ix=0; ix< iNx; ix++){
                
                if(bInit){
                    if(init[X3(ix,iy,iz)]>0.9) init[X3(ix,iy,iz)] = 0.9;
                    if(init[X3(ix,iy,iz)]<0.1) init[X3(ix,iy,iz)] = 0.1;
                    alpha[X3(ix,iy,iz)] = init[X3(ix,iy,iz)];
                }
                else{
                    //if(brainMask[X3(ix,iy,iz)]) alpha[X3(ix,iy,iz)] = pInitAlpha;
                    //else alpha[X3(ix,iy,iz)] = 0.0;
                    alpha[X3(ix,iy,iz)] = pInitAlpha;
                }
            }
    
    
    // INITIALIZATION OF LAMBDA WITH RESPECT TO RANGE OF HR
    EM_joint(im,
          brainMask,
          u,
          
          priorIn,
          priorOut,
          
          iNrOfOtherClasses,
          
          pStdMin,
          pMin,
          pIOffset,
          
          iMaxIterEM,
          pEpsEM,
          
          posteriorIn,
          posteriorOut,
          
          pMeanIn,
          pSigmaIn,
          pMeanOut,
          pSigmaOut,
          pPriorIn,
          
          totalIn,
          totalOut,
          
          verbose,
          0,
             
          bEMreinitialize,
             
          iHyperHypoIntens_toWhich,
          pHyperHypoIntens_howMuch);
    
    pixel_type pMinHr = 1e10;
    pixel_type pMaxHr = -1e10;
    for (int iz=0; iz< iNz; iz++)
        for (int iy=0; iy< iNy; iy++)
            for (int ix=0; ix< iNx; ix++)
            {
                if(brainMask[X3(ix,iy,iz)]){
                    if(totalOut[X3(ix,iy,iz)]<pMin) totalOut[X3(ix,iy,iz)] = pMin;
                    if(totalIn[X3(ix,iy,iz)]<pMin) totalIn[X3(ix,iy,iz)] = pMin;
                    if(alpha[X3(ix,iy,iz)]<pMin) alpha[X3(ix,iy,iz)] = pMin;
                    
                    pixel_type beta = 1.0-alpha[X3(ix,iy,iz)];
                    if(beta<pMin) beta = pMin;
                    
                    pfHr[X3(ix,iy,iz)] = -std::log(totalIn[X3(ix,iy,iz)]) + std::log(totalOut[X3(ix,iy,iz)]) + pCoupling* (std::log(beta) - std::log(alpha[X3(ix,iy,iz)]));
                        
                    if (pfHr[X3(ix,iy,iz)] > pMaxHr) pMaxHr = pfHr[X3(ix,iy,iz)];
                    if (pfHr[X3(ix,iy,iz)] < pMinHr) pMinHr = pfHr[X3(ix,iy,iz)];
                }
            }
    pixel_type pRangeHr = pMaxHr-pMinHr;
    pLambda /= pRangeHr;
    

    // ALGORITHM'S MAIN LOOP
    int k=0;
    while(k<iMaxIter){
        
        
        // STORE U_OLD
        for (int iz=0; iz< iNz; iz++)
            for (int iy=0; iy< iNy; iy++)
                for (int ix=0; ix< iNx; ix++)
                {
                    uOld[X3(ix,iy,iz)] = u[X3(ix,iy,iz)];
                }
        
        
        // EM
        EM_joint(im,
           brainMask,
           u,
           
           priorIn,
           priorOut,
           
           iNrOfOtherClasses,
           
           pStdMin,
           pMin,
           pIOffset,
           
           iMaxIterEM,
           pEpsEM,
           
           posteriorIn,
           posteriorOut,
           
           pMeanIn,
           pSigmaIn,
           pMeanOut,
           pSigmaOut,
           pPriorIn,
           
           totalIn,
           totalOut,
           
           verbose,
           k,
               
           bEMreinitialize,
                 
           iHyperHypoIntens_toWhich,
           pHyperHypoIntens_howMuch);
        
        
        // INCLUSION IMAGE
        basic_image<pixel_type,dimension> inclusion(geo);
        if(iInclusion_toWhich>=0){
            basic_image<pixel_type,dimension> total(geo);
            std::fill(total.begin(),total.end(),0.0);
            
                for(int j=0; j<iNrOfOutClasses; j++){
                    
                    for (int iz=0; iz< iNz; iz++)
                        for (int iy=0; iy< iNy; iy++)
                            for (int ix=0; ix< iNx; ix++)
                            {
                                if(brainMask[X3(ix,iy,iz)]){
                                    total[X3(ix,iy,iz)] += (1.0-u[X3(ix,iy,iz)])*posteriorOut[j][X3(ix,iy,iz)];
                                }
                            }
                    
                }
            
            for (int iz=0; iz< iNz; iz++)
                for (int iy=0; iy< iNy; iy++)
                    for (int ix=0; ix< iNx; ix++)
                    {
                        if(brainMask[X3(ix,iy,iz)]){
                            if(total[X3(ix,iy,iz)]>pMin){
                                inclusion[X3(ix,iy,iz)] = (1.0-u[X3(ix,iy,iz)])*posteriorOut[iInclusion_toWhich][X3(ix,iy,iz)]/total[X3(ix,iy,iz)];
                            }
                            else
                            {
                                inclusion[X3(ix,iy,iz)] += 1.0/(pixel_type) iNrOfOutClasses;
                            }
                        }
                    }
        }
        
        
        // SEG
        SEG_joint(totalIn,
            totalOut,
            alpha,
                             
            brainMask,
            u,
                             
            pMin,
                             
            iMaxIterSeg,
            pEpsSeg,
            
            pLambda,
            pMu,
               
            bUpdateAlpha,
               
            pCoupling,
                             
            basename,
            verbose,
                 
            iInclusion_toWhich,
            inclusion);
        
        
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
        std::cout << "-------------------------OUTER--------------------------" << std::endl;
        std::cout << " End Outer iteration: " << ss.str() << " | Relative change u: " << pDeltaSumU << std::endl;
        std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
        
        if(std::abs(pDeltaSumU)<pEpsSeg) break;
        k++;
        
    }
    
    return;
}

}
}
#endif //BASE_JOINT_HPP_INCLUDED