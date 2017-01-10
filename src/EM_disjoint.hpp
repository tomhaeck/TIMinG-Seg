#ifndef EM_DISJOINT_HPP_INCLUDED
#define EM_DISJOINT_HPP_INCLUDED

#include "_gaussian.hpp"
#include "_help.hpp"

#define SQR(x) (x)*(x)
#define X3(ix,iy,iz) ((iz)*iNy+(iy))*iNx+(ix)

namespace image
{
namespace segmentation
{
    
template<typename pixel_type,unsigned int dimension>
void EM_disjoint(std::vector<image::basic_image<pixel_type,dimension> >& im,
                      image::basic_image<pixel_type,dimension>& brainMask,
                      std::vector<image::basic_image<pixel_type,dimension> >& u,
                      
                      std::vector<basic_image<pixel_type,dimension> >& priorIn,
                      std::vector<basic_image<pixel_type,dimension> >& priorOut,
        
                      int iNrOfOtherClasses,
                      
                      pixel_type pStdMin,
                      pixel_type pMin,
                      pixel_type pIOffset,
                      
                      int iMaxIter,
                      pixel_type pEps,
                      
                      std::vector<std::vector<image::basic_image<pixel_type,dimension> > >& posteriorIn,
                      std::vector<std::vector<image::basic_image<pixel_type,dimension> > >& posteriorOut,
                      
                      std::vector<std::vector<pixel_type> >& pMeanIn,
                      std::vector<std::vector<pixel_type> >& pSigmaIn,
                      std::vector<std::vector<pixel_type> >& pMeanOut,
                      std::vector<std::vector<pixel_type> >& pSigmaOut,
                      std::vector<std::vector<pixel_type> >& pPriorIn,
        
                      std::vector<basic_image<pixel_type,dimension> >& totalIn,
                      std::vector<basic_image<pixel_type,dimension> >& totalOut,
        
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
    
    std::vector<std::vector<pixel_type> > lesionIndicator;
    int iNrOfLesionStates = image::fillLesionIndicator(lesionIndicator, iNrOfChannels);
    
    
    // DECLARE PIXEL_TYPE
    std::vector<pixel_type> Ln(iNrOfLesionStates), L(iNrOfLesionStates), deltaL(iNrOfLesionStates);
    

    // DECLARE STD::VECTOR<PIXEL_TYPE>
    std::vector<std::vector<pixel_type> > pNormalizationIn(iNrOfInClasses,std::vector<pixel_type>(iNrOfChannels));
    std::vector<std::vector<pixel_type> > pNormalizationOut(iNrOfOutClasses,std::vector<pixel_type>(iNrOfChannels));
    
    
    // DECLARE STD:VECTOR<BASIC_IMAGE>
    std::vector<std::vector<basic_image<pixel_type,dimension> > > likelihood_foreground(iNrOfInClasses, std::vector<basic_image<pixel_type,dimension> >(iNrOfChannels, basic_image<pixel_type,dimension>(geo)));
    std::vector<std::vector<basic_image<pixel_type,dimension> > > likelihood_background(iNrOfOutClasses, std::vector<basic_image<pixel_type,dimension> >(iNrOfChannels, basic_image<pixel_type,dimension>(geo)));
    
    std::vector<basic_image<pixel_type,dimension> > new_priorIn(iNrOfInClasses,basic_image<pixel_type,dimension>(geo));
    std::vector<basic_image<pixel_type,dimension> > new_priorOut(iNrOfOutClasses,basic_image<pixel_type,dimension>(geo));
    
    std::vector<basic_image<pixel_type,dimension> > posterior(iNrOfLesionStates,basic_image<pixel_type,dimension>(geo));
    
    
    // INITIALIZE STD::VECTOR<BASIC_IMAGE>
    if((iIter==0) | bEMreinitialize){
        for(int t=0; t<iNrOfLesionStates; t++) for(int j=0; j<iNrOfInClasses; j++) posteriorIn[t][j] = priorIn[j];
        for(int t=0; t<iNrOfLesionStates; t++) for(int j=0; j<iNrOfOutClasses; j++) posteriorOut[t][j] = priorOut[j];
    }
    
    for(int j=0; j<iNrOfInClasses; j++) new_priorIn[j] = priorIn[j];
    for(int j=0; j<iNrOfOutClasses; j++) new_priorOut[j] = priorOut[j];
    
    for(int t=0; t<iNrOfLesionStates; t++){
        std::fill(posterior[t].begin(),posterior[t].end(),1.0);
        for(int l=0; l<iNrOfChannels; l++){
            for (int iz=0; iz< iNz; iz++)
            for (int iy=0; iy< iNy; iy++)
            for (int ix=0; ix< iNx; ix++)
            {
                if(brainMask[X3(ix,iy,iz)]){
                    if(lesionIndicator[t][l]==0) posterior[t][X3(ix,iy,iz)] *= (1.0 - u[l][X3(ix,iy,iz)]);
                    if(lesionIndicator[t][l]==1) posterior[t][X3(ix,iy,iz)] *= u[l][X3(ix,iy,iz)];
                }
                else
                posterior[t][X3(ix,iy,iz)] = 0.0;
                
            }
        }
    }
    
    
    // INITIALIZE LOGLIKELIHOOD
    for(int t=0; t<iNrOfLesionStates; t++){
        L[t] = std::numeric_limits<pixel_type>::max();
    }
    
    
    // START EM LOOP
    std::cout << "-------------------------EM-----------------------------" << std::endl;
    
    int k=0;
    while(k<iMaxIter){
        
        
        // BACKGROUND PARAMETER ESTIMATION (MEAN)
        for(int j=0; j<iNrOfOutClasses; j++){
            for(int l=0; l<iNrOfChannels; l++ ){
                
                pixel_type meanOut = 0.0;
                pixel_type normalizationOut = 0.0;
                for (int iz=0; iz< iNz; iz++)
                    for (int iy=0; iy< iNy; iy++)
                        for (int ix=0; ix< iNx; ix++)
                        {
                            if(brainMask[X3(ix,iy,iz)]){
                                pixel_type sum = 0.0;
                                for(int t=0; t<iNrOfLesionStates; t++){
                                    sum+=posteriorOut[t][j][X3(ix,iy,iz)]*posterior[t][X3(ix,iy,iz)]*(1.0-lesionIndicator[t][l]);
                                }
                                meanOut += im[l][X3(ix,iy,iz)]*sum;
                                normalizationOut += sum;
                            }
                        }
                if (normalizationOut>pMin) meanOut /= normalizationOut;
                else throw std::string( "error. NormalizationOut<=0.0.  Everything is lesion");
                pMeanOut[j][l]=meanOut;
                pNormalizationOut[j][l]=normalizationOut;
            }
        }
        
        
        // BACKGROUND PARAMETER OFFSET
        if(k==0)
            for(int i=0; i<iNrOfOtherClasses; i++)
                for(int l=0; l<iNrOfChannels; l++)
                    pMeanOut[i+iNrOfOutClasses-iNrOfOtherClasses][l] = (pIOffset*i+1.0)*pMeanOut[iNrOfOutClasses-iNrOfOtherClasses][l];
        
        
        // BACKGROUND PARAMETER ESTIMATION (VARIANCE)
        for(int j=0; j<iNrOfOutClasses; j++){
            for(int l=0; l<iNrOfChannels; l++){
                
                pixel_type sigmaOut = 0.0;
                for (int iz=0; iz< iNz; iz++)
                    for (int iy=0; iy< iNy; iy++)
                        for (int ix=0; ix< iNx; ix++)
                        {
                            if(brainMask[X3(ix,iy,iz)]){
                                pixel_type sum = 0.0;
                                for(int t=0; t<iNrOfLesionStates; t++){
                                    sum+=posteriorOut[t][j][X3(ix,iy,iz)]*posterior[t][X3(ix,iy,iz)]*(1.0-lesionIndicator[t][l]);
                                }
                                sigmaOut += (im[l][X3(ix,iy,iz)]-pMeanOut[j][l])*(im[l][X3(ix,iy,iz)]-pMeanOut[j][l])*sum;
                            }
                        }
                if (pNormalizationOut[j][l]>pMin) sigmaOut = std::sqrt(sigmaOut/pNormalizationOut[j][l]);
                if (std::abs(sigmaOut)<pStdMin){
                    sigmaOut = pStdMin;
                }
                pSigmaOut[j][l]=sigmaOut;
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
                                pixel_type sum = 0.0;
                                for(int t=0; t<iNrOfLesionStates; t++){
                                    if(iHyperHypoIntens_toWhich[l]<0 || iIter<=5){
                                        sum+=posteriorIn[t][j][X3(ix,iy,iz)]*posterior[t][X3(ix,iy,iz)]*lesionIndicator[t][l];
                                    }
                                    else{
                                        if(pHyperHypoIntens_howMuch[l]>=0.0 && (im[l][X3(ix,iy,iz)]>(pMeanOut[iHyperHypoIntens_toWhich[l]][l]+pHyperHypoIntens_howMuch[l]*pSigmaOut[iHyperHypoIntens_toWhich[l]][l]))){
                                            sum+=posteriorIn[t][j][X3(ix,iy,iz)]*posterior[t][X3(ix,iy,iz)]*lesionIndicator[t][l];
                                        }
                                        else if(pHyperHypoIntens_howMuch[l]<0.0 && (im[l][X3(ix,iy,iz)]<(pMeanOut[iHyperHypoIntens_toWhich[l]][l]+pHyperHypoIntens_howMuch[l]*pSigmaOut[iHyperHypoIntens_toWhich[l]][l]))){
                                            sum+=posteriorIn[t][j][X3(ix,iy,iz)]*posterior[t][X3(ix,iy,iz)]*lesionIndicator[t][l];
                                        }
                                    }
                                }
                                meanIn += im[l][X3(ix,iy,iz)]*sum;
                                normalizationIn += sum;
                            }
                        }
                if (normalizationIn>pMin) meanIn /= normalizationIn;
                else throw std::string("error.  NormalizationIn<=0.0.  No lesion found in at least one of the channels.");
                pMeanIn[j][l]=meanIn;
                pNormalizationIn[j][l]=normalizationIn;
            }
        }
        
        
        // FOREGROUND PARAMETER OFFSET
        if(k==0)
            for(int i=0; i<iNrOfInClasses-1; i++)
                for(int l=0; l<iNrOfChannels; l++)
                    pMeanIn[i+1][l] = (pIOffset*(i+1)+1.0)*pMeanIn[0][l];
        
        
        // FOREGROUND PARAMETER ESTIMATION (VARIANCE)
        for(int j=0; j<iNrOfInClasses; j++){
            for(int l=0; l<iNrOfChannels; l++){
                
                pixel_type sigmaIn = 0.0;
                for (int iz=0; iz< iNz; iz++)
                    for (int iy=0; iy< iNy; iy++)
                        for (int ix=0; ix< iNx; ix++)
                        {
                            if(brainMask[X3(ix,iy,iz)]){
                                pixel_type sum = 0.0;
                                for(int t=0; t<iNrOfLesionStates; t++){
                                    if(iHyperHypoIntens_toWhich[l]<0 || iIter<=5){
                                        sum+=posteriorIn[t][j][X3(ix,iy,iz)]*posterior[t][X3(ix,iy,iz)]*lesionIndicator[t][l];
                                    }
                                    else{
                                        if(pHyperHypoIntens_howMuch[l]>=0.0 && (im[l][X3(ix,iy,iz)]>(pMeanOut[iHyperHypoIntens_toWhich[l]][l]+pHyperHypoIntens_howMuch[l]*pSigmaOut[iHyperHypoIntens_toWhich[l]][l]))){
                                            sum+=posteriorIn[t][j][X3(ix,iy,iz)]*posterior[t][X3(ix,iy,iz)]*lesionIndicator[t][l];
                                        }
                                        else if(pHyperHypoIntens_howMuch[l]<0.0 && (im[l][X3(ix,iy,iz)]<(pMeanOut[iHyperHypoIntens_toWhich[l]][l]+pHyperHypoIntens_howMuch[l]*pSigmaOut[iHyperHypoIntens_toWhich[l]][l]))){
                                            sum+=posteriorIn[t][j][X3(ix,iy,iz)]*posterior[t][X3(ix,iy,iz)]*lesionIndicator[t][l];
                                        }
                                    }
                                }
                                sigmaIn += (im[l][X3(ix,iy,iz)]-pMeanIn[j][l])*(im[l][X3(ix,iy,iz)]-pMeanIn[j][l])*sum;
                            }
                        }
                if (pNormalizationIn[j][l]>pMin) sigmaIn = std::sqrt(sigmaIn/pNormalizationIn[j][l]);
                else throw std::string("NormalizationIn<=0.0.  No lesion found in at least one of the channels");
                if (std::abs(sigmaIn)<pStdMin){
                    sigmaIn = pStdMin;
                }
                pSigmaIn[j][l]=sigmaIn;
            }
        }
        
        
        // FOREGROUND PARAMETER ESTIMATION (GAUSSIAN WEIGHTING)
        for(int j=0; j<iNrOfInClasses; j++){
            for(int l=0; l<iNrOfChannels; l++){
                
                pixel_type priorIn = 0.0;
                pixel_type normalizationIn = 0.0;
                for (int iz=0; iz< iNz; iz++)
                    for (int iy=0; iy< iNy; iy++)
                        for (int ix=0; ix< iNx; ix++)
                        {
                            if(brainMask[X3(ix,iy,iz)]){
                                pixel_type sum1 = 0.0;
                                pixel_type sum2 = 0.0;
                                for(int t=0; t<iNrOfLesionStates; t++){
                                    if(iHyperHypoIntens_toWhich[l]<0 || iIter<=5){
                                        sum1+=posteriorIn[t][j][X3(ix,iy,iz)]*posterior[t][X3(ix,iy,iz)]*lesionIndicator[t][l];
                                        sum2+=posterior[t][X3(ix,iy,iz)]*lesionIndicator[t][l];
                                    }
                                    else{
                                        if(pHyperHypoIntens_howMuch[l]>=0.0 && (im[l][X3(ix,iy,iz)]>(pMeanOut[iHyperHypoIntens_toWhich[l]][l]+pHyperHypoIntens_howMuch[l]*pSigmaOut[iHyperHypoIntens_toWhich[l]][l]))){
                                            sum1+=posteriorIn[t][j][X3(ix,iy,iz)]*posterior[t][X3(ix,iy,iz)]*lesionIndicator[t][l];
                                            sum2+=posterior[t][X3(ix,iy,iz)]*lesionIndicator[t][l];
                                        }
                                        else if(pHyperHypoIntens_howMuch[l]<0.0 && (im[l][X3(ix,iy,iz)]<(pMeanOut[iHyperHypoIntens_toWhich[l]][l]+pHyperHypoIntens_howMuch[l]*pSigmaOut[iHyperHypoIntens_toWhich[l]][l]))){
                                            sum1+=posteriorIn[t][j][X3(ix,iy,iz)]*posterior[t][X3(ix,iy,iz)]*lesionIndicator[t][l];
                                            sum2+=posterior[t][X3(ix,iy,iz)]*lesionIndicator[t][l];
                                        }
                                    }
                                }
                                priorIn += sum1;
                                normalizationIn += sum2;
                            }
                        }
                if (normalizationIn>pMin) priorIn /= normalizationIn;
                else throw std::string("error.  priorNormalizationIn<=0.0.  One of the Gaussians has weight zero in at least one of the channels");
                pPriorIn[j][l]=priorIn;
            }
        }
        
    
        // UPDATE THE LIKELIHOODS
        for(int j=0; j<iNrOfOutClasses; j++)
            for(int l=0; l<iNrOfChannels; l++)
                image::evaluate_in_gaussian(im[l],likelihood_background[j][l], pMeanOut[j][l], pSigmaOut[j][l]);
        for(int j=0; j<iNrOfInClasses; j++)
            for(int l=0; l<iNrOfChannels; l++)
                image::evaluate_in_gaussian(im[l],likelihood_foreground[j][l], pMeanIn[j][l], pSigmaIn[j][l]);
    
        
        // FOREGROUND POSTERIORS
        for(int t=0; t<iNrOfLesionStates; t++){
            
            std::fill(totalIn[t].begin(),totalIn[t].end(),0.0);
            for(int j=0; j<iNrOfInClasses; j++){
                for (int iz=0; iz< iNz; iz++)
                    for (int iy=0; iy< iNy; iy++)
                        for (int ix=0; ix< iNx; ix++)
                        {
                            if(brainMask[X3(ix,iy,iz)]){
                                pixel_type prod = 1.0;
                                for(int i=0; i<iNrOfChannels; i++){
                                    if(lesionIndicator[t][i]==1) prod *= likelihood_foreground[j][i][X3(ix,iy,iz)]*pPriorIn[j][i];
                                }
                                totalIn[t][X3(ix,iy,iz)] += prod;
                            }
                        }
            }
            
            for (int iz=0; iz< iNz; iz++)
                for (int iy=0; iy< iNy; iy++)
                    for (int ix=0; ix< iNx; ix++)
                    {
                        if(brainMask[X3(ix,iy,iz)]){
                            if(totalIn[t][X3(ix,iy,iz)]>pMin){
                                for(int j=0; j<iNrOfInClasses; j++){
                                    pixel_type prod = 1.0;
                                    for(int i=0; i<iNrOfChannels; i++){
                                        if(lesionIndicator[t][i]==1) prod *= likelihood_foreground[j][i][X3(ix,iy,iz)]*pPriorIn[j][i];
                                    }
                                    posteriorIn[t][j][X3(ix,iy,iz)] = prod/totalIn[t][X3(ix,iy,iz)];
                                }
                            }
                            else
                            {
                                for(int j=0; j<iNrOfInClasses; j++){
                                    posteriorIn[t][j][X3(ix,iy,iz)] = 1.0/(pixel_type) iNrOfInClasses;
                                }
                            }
                        }
                    }
        }
    
        
        // BACKGROUND POSTERIORS
        for(int t=0; t<iNrOfLesionStates; t++){
            
            std::fill(totalOut[t].begin(),totalOut[t].end(),0.0);
            for(int j=0; j<iNrOfOutClasses; j++){

                for (int iz=0; iz< iNz; iz++)
                    for (int iy=0; iy< iNy; iy++)
                        for (int ix=0; ix< iNx; ix++)
                        {
                            if(brainMask[X3(ix,iy,iz)]){
                                pixel_type prod = 1.0;
                                for(int i=0; i<iNrOfChannels; i++){
                                    if(lesionIndicator[t][i]==0) prod *= likelihood_background[j][i][X3(ix,iy,iz)];
                                }
                                totalOut[t][X3(ix,iy,iz)] += prod * new_priorOut[j][X3(ix,iy,iz)];
                            }
                        }
                
            }
            
            
            for (int iz=0; iz< iNz; iz++)
                for (int iy=0; iy< iNy; iy++)
                    for (int ix=0; ix< iNx; ix++)
                    {
                        if(brainMask[X3(ix,iy,iz)]){
                            if(totalOut[t][X3(ix,iy,iz)]>pMin){
                                for(int j=0; j<iNrOfOutClasses; j++){
                                    pixel_type prod = 1.0;
                                    for(int i=0; i<iNrOfChannels; i++){
                                        if(lesionIndicator[t][i]==0) prod *= likelihood_background[j][i][X3(ix,iy,iz)];
                                    }
                                    posteriorOut[t][j][X3(ix,iy,iz)] = prod*new_priorOut[j][X3(ix,iy,iz)]/totalOut[t][X3(ix,iy,iz)];
                                }
                            }
                            else
                            {
                                for(int j=0; j<iNrOfOutClasses; j++){
                                    posteriorOut[t][j][X3(ix,iy,iz)] = 1.0/(pixel_type) iNrOfOutClasses;
                                }
                            }
                        }
                    }
        }        

        
        for(int t=0; t<iNrOfLesionStates; t++){
            
            for (int iz=0; iz< iNz; iz++)
                for (int iy=0; iy< iNy; iy++)
                    for (int ix=0; ix< iNx; ix++)
                    {
                        if(!brainMask[X3(ix,iy,iz)]){
                            if(t==0){
                                totalOut[t][X3(ix,iy,iz)] = 0.9;
                                totalIn[t][X3(ix,iy,iz)] = 1.0;
                            }
                            else if(t==iNrOfLesionStates-1){
                                totalOut[t][X3(ix,iy,iz)] = 1.0;
                                totalIn[t][X3(ix,iy,iz)] = 0.1;
                            }
                            else{
                                totalOut[t][X3(ix,iy,iz)] = 0.0;
                                totalIn[t][X3(ix,iy,iz)] = 0.0;
                            }
                        }
                    }
        }
    
        // CHECK CONVERGENCE CRITERIUM
        for(int t=0; t<iNrOfLesionStates; t++){
            Ln[t] = 0.0;
            for (int iz=0; iz< iNz; iz++)
                for (int iy=0; iy< iNy; iy++)
                    for (int ix=0; ix< iNx; ix++)
                    {
                        if(brainMask[X3(ix,iy,iz)]){
                            pixel_type total = totalIn[t][X3(ix,iy,iz)]*totalOut[t][X3(ix,iy,iz)];
                            if(total > pMin)  Ln[t] += -std::log(total);
                        }
                    }
            deltaL[t] = (L[t]-Ln[t])/Ln[0]*100;
            L[t] = Ln[t];
        }
        
        
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
        for(int t=0; t<iNrOfLesionStates; t++) if(std::abs(deltaL[t])>=pEps) converged = false;
        if(converged) break;
        
        k++;
        
    }
    return;
}
    
}
}
#endif //EM_DSIJOINT_HPP_INCLUDED