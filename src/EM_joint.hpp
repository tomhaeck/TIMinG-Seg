#ifndef EM_JOINT_HPP_INCLUDED
#define EM_JOINT_HPP_INCLUDED

#include "_gaussian.hpp"

#define SQR(x) (x)*(x)
#define X3(ix,iy,iz) ((iz)*iNy+(iy))*iNx+(ix)

namespace image
{
namespace segmentation
{
    
template<typename pixel_type,unsigned int dimension>
void EM_joint(std::vector<image::basic_image<pixel_type,dimension> >& im,
                  image::basic_image<pixel_type,dimension>& brainMask,
                  image::basic_image<pixel_type,dimension>& u,
                  
                  std::vector<basic_image<pixel_type,dimension> >& priorIn,
                  std::vector<basic_image<pixel_type,dimension> >& priorOut,

                  int iNrOfOtherClasses,
                  
                  pixel_type pStdMin,
                  pixel_type pMin,
                  pixel_type pIOffset,
                  
                  int iMaxIter,
                  pixel_type pEps,
                  
                  std::vector<image::basic_image<pixel_type,dimension> >& posteriorIn,
                  std::vector<image::basic_image<pixel_type,dimension> >& posteriorOut,
                  
                  std::vector<std::vector<pixel_type> >& pMeanIn,
                  std::vector<std::vector<pixel_type> >& pSigmaIn,
                  std::vector<std::vector<pixel_type> >& pMeanOut,
                  std::vector<std::vector<pixel_type> >& pSigmaOut,
                  std::vector<pixel_type>& pPriorIn,

                  basic_image<pixel_type,dimension>& totalIn,
                  basic_image<pixel_type,dimension>& totalOut,

                  bool verbose,
                  int iIter,

                  bool bEMreinitialize,
              
                  std::vector<int> iHyperHypoIntens_toWhich,
                  std::vector<pixel_type> pHyperHypoIntens_howMuch){
    
    
    // DECLARE INT
    int iNrOfChannels = im.size();
    
    int iNrOfInClasses  = priorIn.size();
    int iNrOfOutClasses = priorOut.size();
    
    geometry<dimension> geo = im[0].geometry();
    int iNx = geo[0];
    int iNy = geo[1];
    int iNz = geo[2];
    
    
    // DECLARE PIXEL_TYPE
    std::vector<pixel_type> Ln(2),L(2),deltaL(2);
    
    
    // DECLARE STD::VECTOR<PIXEL_TYPE>
    std::vector<std::vector<pixel_type> > pNormalizationIn(iNrOfInClasses,std::vector<pixel_type>(iNrOfChannels));
    std::vector<pixel_type> pNormalizationOut(iNrOfOutClasses);
    
    
    // DECLARE BASIC_IMAGE
    basic_image<pixel_type,dimension> total(geo);
    
    
    // DECLARE STD:VECTOR<BASIC_IMAGE>
    std::vector<basic_image<pixel_type,dimension> > likelihood_foreground(iNrOfInClasses, basic_image<pixel_type,dimension>(geo));
    std::vector<basic_image<pixel_type,dimension> > likelihood_background(iNrOfOutClasses, basic_image<pixel_type,dimension>(geo));
    
    std::vector<basic_image<pixel_type,dimension> > new_priorIn(iNrOfInClasses,basic_image<pixel_type,dimension>(geo));
    std::vector<basic_image<pixel_type,dimension> > new_priorOut(iNrOfOutClasses,basic_image<pixel_type,dimension>(geo));

    std::vector<gmm::dense_matrix<pixel_type> > covarIn(iNrOfInClasses,gmm::dense_matrix<pixel_type>(iNrOfChannels,iNrOfChannels));
    std::vector<gmm::dense_matrix<pixel_type> > covarOut(iNrOfOutClasses,gmm::dense_matrix<pixel_type>(iNrOfChannels,iNrOfChannels));
    
    
    // INITIALIZE STD::VECTOR<BASIC_IMAGE>
    if((iIter==0) | bEMreinitialize){
        for(int j=0; j<iNrOfInClasses; j++) posteriorIn[j] = priorIn[j];
        for(int j=0; j<iNrOfOutClasses; j++) posteriorOut[j] = priorOut[j];
    }
    
    for(int j=0; j<iNrOfInClasses; j++) new_priorIn[j] = priorIn[j];
    for(int j=0; j<iNrOfOutClasses; j++) new_priorOut[j] = priorOut[j];
    
    
    // INITIALIZE LOGLIKELIHOOD
    L[0] = std::numeric_limits<pixel_type>::max();
    L[1] = std::numeric_limits<pixel_type>::max();
    
    
    // START EM LOOP
    std::cout << "-------------------------EM-----------------------------" << std::endl;
    
    int k=0;
    while(k<iMaxIter){
        
        
        // BACKGROUND PARAMETER ESTIMATION (MEAN)
        for(int j=0; j<iNrOfOutClasses; j++){
            pixel_type normalizationOut = 0.0;
            for (int iz=0; iz< iNz; iz++)
                for (int iy=0; iy< iNy; iy++)
                    for (int ix=0; ix< iNx; ix++)
                    {
                        if(brainMask[X3(ix,iy,iz)]){
                            normalizationOut += posteriorOut[j][X3(ix,iy,iz)]*(1.0-u[X3(ix,iy,iz)]);
                        }
                    }
            
            for(int l=0; l<iNrOfChannels; l++ ){
                pixel_type meanOut = 0.0;
                for (int iz=0; iz< iNz; iz++)
                    for (int iy=0; iy< iNy; iy++)
                        for (int ix=0; ix< iNx; ix++)
                        {
                            if(brainMask[X3(ix,iy,iz)]){
                                meanOut += im[l][X3(ix,iy,iz)]*posteriorOut[j][X3(ix,iy,iz)]*(1.0-u[X3(ix,iy,iz)]);
                            }
                        }
                if (normalizationOut>pMin) meanOut /= normalizationOut;
                else throw std::string("error. NormalizationOut<=0.0.  Everything is lesion");
                pMeanOut[j][l]=meanOut;
            }
            pNormalizationOut[j]=normalizationOut;
        }
        
        
        // BACKGROUND PARAMETER OFFSET
        if(k==0)
            for(int i=0; i<iNrOfOtherClasses; i++)
                for(int l=0; l<iNrOfChannels; l++)
                    pMeanOut[i+iNrOfOutClasses-iNrOfOtherClasses][l] = (pIOffset*i+1.0)*pMeanOut[iNrOfOutClasses-iNrOfOtherClasses][l];
        
        
        // BACKGROUND PARAMETER ESTIMATION (COVARIANCE)
        for(int j=0; j<iNrOfOutClasses; j++){
            for(int l=0; l<iNrOfChannels; l++){
                for(int m=l; m<iNrOfChannels; m++){
                    pixel_type sigmaOut = 0.0;
                    for (int iz=0; iz< iNz; iz++)
                        for (int iy=0; iy< iNy; iy++)
                            for (int ix=0; ix< iNx; ix++)
                            {
                                if(brainMask[X3(ix,iy,iz)]){
                                    sigmaOut += (im[l][X3(ix,iy,iz)]-pMeanOut[j][l])*(im[m][X3(ix,iy,iz)]-pMeanOut[j][m])*posteriorOut[j][X3(ix,iy,iz)]*(1.0-u[X3(ix,iy,iz)]);
                                }
                            }
                    if (pNormalizationOut[j]>pMin) sigmaOut /= pNormalizationOut[j];
                    if (l==m && std::abs(sigmaOut)<pStdMin){
                        sigmaOut = pStdMin;
                    }
                    covarOut[j](l,m)=sigmaOut;
                    covarOut[j](m,l)=sigmaOut;
                }
            }
        }
    
        
        // FOREGROUND PARAMETER ESTIMATION (MEAN)
        for(int j=0; j<iNrOfInClasses; j++){
           for(int l=0; l<iNrOfChannels; l++){
               
               pixel_type meanIn = 0.0;
               pixel_type normalizationIn = 0.0;
               for (int iz=0; iz< iNz; iz++)
                   for (int iy=0; iy< iNy; iy++)
                       for (int ix=0; ix< iNx; ix++)
                       {
                           if(brainMask[X3(ix,iy,iz)]){
                               if(iHyperHypoIntens_toWhich[l]<0){
                                   meanIn += im[l][X3(ix,iy,iz)]*posteriorIn[j][X3(ix,iy,iz)]*u[X3(ix,iy,iz)];
                                   normalizationIn += posteriorIn[j][X3(ix,iy,iz)]*u[X3(ix,iy,iz)];
                               }
                               else{
                                   if(pHyperHypoIntens_howMuch[l]>=0.0 && (im[l][X3(ix,iy,iz)]>(pMeanOut[iHyperHypoIntens_toWhich[l]][l]+pHyperHypoIntens_howMuch[l]*std::sqrt(pSigmaOut[iHyperHypoIntens_toWhich[l]][l])))){
                                       meanIn += im[l][X3(ix,iy,iz)]*posteriorIn[j][X3(ix,iy,iz)]*u[X3(ix,iy,iz)];
                                       normalizationIn += posteriorIn[j][X3(ix,iy,iz)]*u[X3(ix,iy,iz)];
                                   }
                                   else if(pHyperHypoIntens_howMuch[l]<0.0 && (im[l][X3(ix,iy,iz)]<(pMeanOut[iHyperHypoIntens_toWhich[l]][l]+pHyperHypoIntens_howMuch[l]*std::sqrt(pSigmaOut[iHyperHypoIntens_toWhich[l]][l])))){
                                       meanIn += im[l][X3(ix,iy,iz)]*posteriorIn[j][X3(ix,iy,iz)]*u[X3(ix,iy,iz)];
                                       normalizationIn += posteriorIn[j][X3(ix,iy,iz)]*u[X3(ix,iy,iz)];
                                   }
                               }
                           }
                        }
               if (normalizationIn>pMin) meanIn /= normalizationIn;
               else throw std::string("error. NormalizationIn<=0.0.  No lesion found in at least one of the channels");
               pMeanIn[j][l]=meanIn;
               pNormalizationIn[j][l]=normalizationIn;
           }
        }
        
        // FOREGROUND PARAMETER OFFSET
        if(k==0)
            for(int i=0; i<iNrOfInClasses-1; i++)
                for(int l=0; l<iNrOfChannels; l++)
                    pMeanIn[i+1][l] = (pIOffset*(i+1)+1.0)*pMeanIn[0][l];
        
        
        // FOREGROUND PARAMETER ESTIMATION (COVARIANCE)
        for(int j=0; j<iNrOfInClasses; j++){
            for(int l=0; l<iNrOfChannels; l++){
                for(int m=l; m<iNrOfChannels; m++){
                    pixel_type sigmaIn = 0.0;
                    for (int iz=0; iz< iNz; iz++)
                        for (int iy=0; iy< iNy; iy++)
                            for (int ix=0; ix< iNx; ix++)
                            {
                                if(brainMask[X3(ix,iy,iz)]){
                                    if(iHyperHypoIntens_toWhich[l]<0){
                                        sigmaIn += (im[l][X3(ix,iy,iz)]-pMeanIn[j][l])*(im[m][X3(ix,iy,iz)]-pMeanIn[j][m])*posteriorIn[j][X3(ix,iy,iz)]*u[X3(ix,iy,iz)];
                                    }
                                    else{
                                        if(pHyperHypoIntens_howMuch[l]>=0.0 && (im[l][X3(ix,iy,iz)]>(pMeanOut[iHyperHypoIntens_toWhich[l]][l]+pHyperHypoIntens_howMuch[l]*std::sqrt(pSigmaOut[iHyperHypoIntens_toWhich[l]][l])))){
                                                sigmaIn += (im[l][X3(ix,iy,iz)]-pMeanIn[j][l])*(im[m][X3(ix,iy,iz)]-pMeanIn[j][m])*posteriorIn[j][X3(ix,iy,iz)]*u[X3(ix,iy,iz)];
                                        }
                                        else if(pHyperHypoIntens_howMuch[l]<0.0 && (im[l][X3(ix,iy,iz)]<(pMeanOut[iHyperHypoIntens_toWhich[l]][l]+pHyperHypoIntens_howMuch[l]*std::sqrt(pSigmaOut[iHyperHypoIntens_toWhich[l]][l])))){
                                                sigmaIn += (im[l][X3(ix,iy,iz)]-pMeanIn[j][l])*(im[m][X3(ix,iy,iz)]-pMeanIn[j][m])*posteriorIn[j][X3(ix,iy,iz)]*u[X3(ix,iy,iz)];
                                        }
                                    }
                                }
                            }
                    if (pNormalizationIn[j][l]>pMin) sigmaIn /= pNormalizationIn[j][l];
                    else throw std::string("NormalizationIn<=0.0.  No lesion found");
                    if (l==m && std::abs(sigmaIn)<pStdMin){
                        sigmaIn = pStdMin;
                    }
                    covarIn[j](l,m)=sigmaIn;
                    covarIn[j](m,l)=sigmaIn;
                }
            }
        }
        
        
        // FOREGROUND PARAMETER ESTIMATION (GAUSSIAN WEIGHTING)
        for(int j=0; j<iNrOfInClasses; j++){
                pixel_type priorIn = 0.0;
                pixel_type priorNormalizationIn = 0.0;
                for (int iz=0; iz<iNz; iz++)
                    for (int ix=0; ix< iNx; ix++)
                        for (int iy=0; iy< iNy; iy++)
                        {
                            if(brainMask[X3(ix,iy,iz)]){
                                priorIn += posteriorIn[j][X3(ix,iy,iz)]*u[X3(ix,iy,iz)];
                                priorNormalizationIn += u[X3(ix,iy,iz)];
                            }
                        }
                
                if (priorNormalizationIn>pMin) priorIn /= priorNormalizationIn;
                else throw std::string("priorNormalizationIn<=0.0.  One of the Gaussians has weight zero");
                pPriorIn[j] = priorIn;
        }
        
        
        // UPDATE THE LIKELIHOODS
        for(int j=0; j<iNrOfOutClasses; j++) image::evaluate_in_gaussian_mult(im,likelihood_background[j], pMeanOut[j], covarOut[j]);
        for(int j=0; j<iNrOfInClasses; j++)  image::evaluate_in_gaussian_mult(im,likelihood_foreground[j], pMeanIn[j], covarIn[j]);
        

        // FOREGROUND POSTERIORS
        std::fill(totalIn.begin(),totalIn.end(),0.0);
        for(int j=0; j<iNrOfInClasses; j++){
            for (int iz=0; iz< iNz; iz++)
                for (int iy=0; iy< iNy; iy++)
                    for (int ix=0; ix< iNx; ix++)
                    {
                        if(brainMask[X3(ix,iy,iz)]){
                            totalIn[X3(ix,iy,iz)] += likelihood_foreground[j][X3(ix,iy,iz)]*pPriorIn[j];
                        }
                    }
        }
        
        for (int iz=0; iz< iNz; iz++)
            for (int iy=0; iy< iNy; iy++)
                for (int ix=0; ix< iNx; ix++)
                {
                    if(brainMask[X3(ix,iy,iz)]){
                        if(totalIn[X3(ix,iy,iz)]>pMin){
                            for(int j=0; j<iNrOfInClasses; j++){
                                posteriorIn[j][X3(ix,iy,iz)] = likelihood_foreground[j][X3(ix,iy,iz)]*pPriorIn[j]/totalIn[X3(ix,iy,iz)];
                            }
                        }
                        else
                        {
                            for(int j=0; j<iNrOfInClasses; j++){
                                posteriorIn[j][X3(ix,iy,iz)] = 1.0/(pixel_type) iNrOfInClasses;
                            }
                        }
                    }
                }


        // BACKGROUND POSTERIORS
        std::fill(totalOut.begin(),totalOut.end(),0.0);
        for(int j=0; j<iNrOfOutClasses; j++){

            for (int iz=0; iz< iNz; iz++)
                for (int iy=0; iy< iNy; iy++)
                    for (int ix=0; ix< iNx; ix++)
                    {
                        if(brainMask[X3(ix,iy,iz)]){
                            totalOut[X3(ix,iy,iz)] += likelihood_background[j][X3(ix,iy,iz)]*new_priorOut[j][X3(ix,iy,iz)];
                        }
                    }
        }
        
        for (int iz=0; iz< iNz; iz++)
            for (int iy=0; iy< iNy; iy++)
                for (int ix=0; ix< iNx; ix++)
                {
                    if(brainMask[X3(ix,iy,iz)]){
                        if(totalOut[X3(ix,iy,iz)]>pMin){
                            for(int j=0; j<iNrOfOutClasses; j++){
                                posteriorOut[j][X3(ix,iy,iz)] = likelihood_background[j][X3(ix,iy,iz)]*new_priorOut[j][X3(ix,iy,iz)]/totalOut[X3(ix,iy,iz)];
                            }
                        }
                        else
                        {
                            for(int j=0; j<iNrOfOutClasses; j++){
                                posteriorOut[j][X3(ix,iy,iz)] = 1.0/(pixel_type) iNrOfOutClasses;
                            }
                        }
                    }
                }
        
        
        for (int iz=0; iz< iNz; iz++)
            for (int iy=0; iy< iNy; iy++)
                for (int ix=0; ix< iNx; ix++)
                {
                    if(!brainMask[X3(ix,iy,iz)]){
                        totalOut[X3(ix,iy,iz)] = 0.9;
                        totalIn[X3(ix,iy,iz)] = 0.1;
                    }
                }
        

        // CHECK CONVERGENCE CRITERIUM
        Ln[0] = 0.0;
        for (int iz=0; iz< iNz; iz++)
            for (int iy=0; iy< iNy; iy++)
                for (int ix=0; ix< iNx; ix++)
                {
                    if(brainMask[X3(ix,iy,iz)]){
                        pixel_type total = totalIn[X3(ix,iy,iz)];
                        if(total > pMin)  Ln[0] += -std::log(total);
                    }
                }
        deltaL[0] = (L[0]-Ln[0])/Ln[0]*100;
        L[0] = Ln[0];
        
        Ln[1] = 0.0;
        for (int iz=0; iz< iNz; iz++)
            for (int iy=0; iy< iNy; iy++)
                for (int ix=0; ix< iNx; ix++)
                {
                    if(brainMask[X3(ix,iy,iz)]){
                        pixel_type total = totalOut[X3(ix,iy,iz)];
                        if(total > pMin)  Ln[1] += -std::log(total);
                    }
                }
        deltaL[1] = (L[1]-Ln[1])/Ln[1]*100;
        L[1] = Ln[1];
        

        // REWRITE SIGMA
        for(int i=0; i<iNrOfOutClasses; i++)
            for(int k=0; k<iNrOfChannels; k++)
                for(int l=0; l<iNrOfChannels; l++)
                    pSigmaOut[i][k*iNrOfChannels+l] = covarOut[i](k,l);
        for(int i=0; i<iNrOfInClasses; i++)
            for(int k=0; k<iNrOfChannels; k++)
                for(int l=0; l<iNrOfChannels; l++)
                    pSigmaIn[i][k*iNrOfChannels+l] = covarIn[i](k,l);
        
        
        // SCREEN FEEDBACK
        std::stringstream ss; ss << k;
        std::cout << " End EM iteration: " << ss.str() << " | Relative change loglikelihood: " << deltaL << std::endl;
        
        if(verbose){
        std::cout << " * Mu foreground: " << pMeanIn << std::endl;
        std::cout << " * Sigma foreground: " << pSigmaIn << std::endl;
        
        std::cout << " * Mu background: " << pMeanOut << std::endl;
        std::cout << " * Sigma background: " << pSigmaOut << std::endl;
        }
      
        bool converged = true;
        if(std::abs(deltaL[0])>=pEps) converged = false;
        if(std::abs(deltaL[1])>=pEps) converged = false;
        if(converged) break;
        
        k++;
        
    }    
    return;
}
    
}
}
#endif //EM_JOINT_HPP_INCLUDED