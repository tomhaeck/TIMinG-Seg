#include <sys/stat.h>

#ifdef _WIN32
    #include <windows.h>
#endif

#include <tclap/CmdLine.h>

#include "BASE_joint.hpp"
#include "BASE_disjoint.hpp"

int main(int argc,char *argv[])
{
    
    // DECLARE INPUT
    std::vector<std::string> inputname;
    int iNrOfChannels;
    
    // DECLARE PRIORS
    std::vector<std::string> priorname;
    int iNrOfPriors;
    
    // DECLARE OUTPUT
    std::string output;
    
    // DECLARE PARAMETERS
    int iNrOfOtherClasses;
    int iNrOfInsideClasses;
        double pIMax;
        double pStdMin;
        double pMin;
        double pIOffset;
    int iMaxIter;
    int iMaxIterEM;
        int iMaxIterSeg;
        double pEpsEM;
        double pEpsSeg;
        double pLambda;
    double pMu;
    double pInitAlpha;
    bool bUpdateAlpha;
        double pCoupling;
        bool bEMreinitialize;
    bool bVerbose;
    std::string maskname;
        bool bInit;
        std::string initname;
    std::vector<int> iHyperHypoIntens_toWhich;
    std::vector<double> pHyperHypoIntens_howMuch;
    int iInclusion_toWhich;
    bool bDisjoint;
    

    // INIT PARAMETERS
    try {
        
        TCLAP::CmdLine cmd("Implements the lesion segmentation method as described by Haeck et al.", ' ', "0.0");
        
        TCLAP::MultiArg<std::string> inputArg("i","channel","Input Channel.  Multiple input channels may be provided by repeating this flag.  Typical input channels are Flair.nii, T2.nii, etc. ",true,"string");
        TCLAP::MultiArg<std::string> priorArg("p","prior","Prior.  Multiple priors may be provided by repeating this flag.  Typical priors are white.nii, grey.nii and csf.nii.",false,"string");
        TCLAP::ValueArg<std::string> outputArg("o","output","Output name.  Example: seg.nii",false,"seg.nii","string");
    
        TCLAP::ValueArg<int> iNrOfOtherClassesArg("","nrOfAdditionalPriors","Number of additional spatially flat priors that are added to the healthy model.",false,0,"int");
        TCLAP::ValueArg<int> iNrOfInsideClassesArg("","nrOfLesionGaussians","Number of Gaussians to model the lesion.",false,1,"int");
        TCLAP::ValueArg<int> iMaxIterArg("","maxIter","Maximum number of iterations of the outer loop.",false,15,"int");
        TCLAP::ValueArg<int> iMaxIterEMArg("","maxIterEM","Maximum number of iterations of the EM loop.",false,10,"int");
        TCLAP::ValueArg<double> pMuArg("","mu","Spatial regularization parameter Mu.",false,1e1,"double");
        TCLAP::ValueArg<double> pInitAlphaArg("","alpha","Initial value for the lesion-healthy mixture coefficient alpha.",false,0.5,"double");
        
        TCLAP::ValueArg<std::string> maskArg("","mask","Name of nifti-image that contains a brain mask.",false,"mask.nii","string");
        
        TCLAP::MultiArg<int> iHyperHypoIntens_toWhichArg("","hyp","",false,"int");
        TCLAP::MultiArg<double> pHyperHypoIntens_howMuchArg("","hypamount","",false,"double");
        TCLAP::ValueArg<int> iInclusion_toWhichArg("","inclusion","This integer indicates the number of the prior in which the lesion should be included.  Indexing starts from zero.  Example: if 'white.nii' is the first prior that is provided on the command line and if you want the lesion to be included in the white matter, set this parameter to '0'.  Default value of this parameter is '-1' and no inclusion constraint will be enforced.",false,-1,"int");
        
        TCLAP::SwitchArg bVerboseArg("","verbose","If set to TRUE, the verbose level is increased.", cmd, false);
        TCLAP::SwitchArg bDisjointArg("","disjoint","If set to TRUE, the channels are handled independently.  Set this to TRUE if the channels provide complementary information on the lesion.", cmd, false);
        TCLAP::SwitchArg bUpdateAlphaArg("","updateAlpha","If set to TRUE, alpha is updated for each iteration according to the current volume fraction of the lesion.", cmd, false);
        cmd.add(iInclusion_toWhichArg);
        cmd.add(pHyperHypoIntens_howMuchArg);
        cmd.add(iHyperHypoIntens_toWhichArg);
        cmd.add(maskArg);
        cmd.add(pInitAlphaArg);
        cmd.add(pMuArg);
        cmd.add(iMaxIterEMArg);
        cmd.add(iMaxIterArg);
        cmd.add(iNrOfInsideClassesArg);
        cmd.add(iNrOfOtherClassesArg);
        cmd.add(outputArg);
        cmd.add(priorArg);
        cmd.add(inputArg);
        
        // Parse the argv array.
        cmd.parse( argc, argv );
        
        // INIT INPUT
        inputname = inputArg.getValue();
        iNrOfChannels = inputname.size();
        
        // INIT PRIORS
        priorname = priorArg.getValue();
        iNrOfPriors = priorname.size();
        
        // INIT OUTPUT
        output = outputArg.getValue();
        
        // INIT PARAMETERS
        iNrOfOtherClasses = iNrOfOtherClassesArg.getValue();
        iNrOfInsideClasses = iNrOfInsideClassesArg.getValue();
            pIMax = 256.0;
            pStdMin = 1.0;
            pMin = 1e-20;
            pIOffset = 1.0;
        iMaxIter = iMaxIterArg.getValue();
        iMaxIterEM = iMaxIterEMArg.getValue();
            iMaxIterSeg = 1;
            pEpsEM = 1e-3;
            pEpsSeg = 1e-12;
            pLambda = 1e1;
        pMu = pMuArg.getValue();
        pInitAlpha = pInitAlphaArg.getValue();
        bUpdateAlpha = bUpdateAlphaArg.getValue();
            pCoupling = 1.0;
            bEMreinitialize = false;
        bVerbose = bVerboseArg.getValue();
        maskname = maskArg.getValue();
            bInit = false;
            initname = std::string("") ;
        iHyperHypoIntens_toWhich = iHyperHypoIntens_toWhichArg.getValue();
        pHyperHypoIntens_howMuch = pHyperHypoIntens_howMuchArg.getValue();
        if(iHyperHypoIntens_toWhich.size()<1){
            iHyperHypoIntens_toWhich.resize(iNrOfChannels);
            pHyperHypoIntens_howMuch.resize(iNrOfChannels);
            std::fill(iHyperHypoIntens_toWhich.begin(),iHyperHypoIntens_toWhich.end(),-1);
            std::fill(pHyperHypoIntens_howMuch.begin(),pHyperHypoIntens_howMuch.end(),0);
        }
        iInclusion_toWhich = iInclusion_toWhichArg.getValue();
        bDisjoint = bDisjointArg.getValue();
        
    } catch (TCLAP::ArgException &e)
    { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }


    // REPORT PARAMETERS:
    std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "* Number of Channels: " << iNrOfChannels << std::endl;
    for(int i=0; i<iNrOfChannels; i++)
    std::cout << "Channel " << i << ": " << inputname[i] << std::endl;
    std::cout << "* Number of Classes: " << iNrOfPriors << std::endl;
    for(int i=0; i<priorname.size(); i++)
    std::cout << "Prior " << i << ": " << priorname[i] << std::endl;
    if(priorname.size()>=1)
    std::cout << " (+ " << iNrOfOtherClasses << " rest classes)" << std::endl;
  
    std::cout << "iNrOfOtherClasses: " << iNrOfOtherClasses << std::endl;
    std::cout << "iNrOfInsideClasses: " << iNrOfInsideClasses << std::endl;
    std::cout << "iNrOfIter: " << iMaxIter << std::endl;
    std::cout << "pMu: " << pMu << std::endl;
    std::cout << "pInitAlpha: " << pInitAlpha << std::endl;
    std::cout << "bUpdateAlpha: " << bUpdateAlpha << std::endl;
    std::cout << "basename: " << output << std::endl;
    std::cout << "bVerbose: " << bVerbose << std::endl;
    std::cout << "maskname: " << maskname << std::endl;
    std::cout << "HyperHypoIntens_toWhich: " << iHyperHypoIntens_toWhich << std::endl;
    std::cout << "HyperHypoIntens_howMuch: " << pHyperHypoIntens_howMuch << std::endl;
    std::cout << "Inclusion_toWhich: " << iInclusion_toWhich << std::endl;
    std::cout << "bDisjoint: " << bDisjoint << std::endl;
    
    std::cout << "--------------------------------------------------------------------------------------------" << std::endl;
   
    
    //////////
    //LOAD
    /////////
    
    image::io::nifti nifti0;
    
    // MAKE DIRECTORY
    std::string basename = output.substr(0,output.size()-4);
    #ifdef _WIN32
        CreateDirectory(basename.c_str(),NULL);
    #elif __APPLE__
        mkdir(basename.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    #else
        #error "Unknown Compiler"
    #endif
    
    // LOAD INPUT
    std::vector<image::basic_image<double,3> > input(iNrOfChannels);
    for(int i=0; i<iNrOfChannels; i++)
        image::serialization::nii_in<double,3>(input[i],inputname[i],nifti0);
    
    int iNrOfOutsideClasses = iNrOfPriors + iNrOfOtherClasses;
    
    // LOAD MASK
    image::basic_image<double,3> mask;
    image::serialization::nii_in<double,3>(mask,maskname,nifti0);
    
    // LOAD INIT
    image::basic_image<double,3> init;
    if(bInit) image::serialization::nii_in<double,3>(init,initname,nifti0);
    
    // LOAD PRIORS OUT
    std::vector<image::basic_image<double,3> > priorsOut(priorname.size());
    for(int i=0; i<iNrOfPriors; i++){
        image::serialization::nii_in<double,3>(priorsOut[i],priorname[i],nifti0);
    }
    
    // ELIMINATE AREAS IN MASK WHERE SUM OF PRIORS OUT IS NEARLY ZERO
    image::geometry<3> geo = input[0].geometry();
    int iNx = geo[0];
    int iNy = geo[1];
    int iNz = geo[2];
    for (int iz=0; iz< iNz; iz++)
        for (int iy=0; iy< iNy; iy++)
            for (int ix=0; ix< iNx; ix++){
                double sum = 0.0;
                for(int i=0; i<iNrOfPriors; i++) sum += priorsOut[i][X3(ix,iy,iz)];
                if(sum<1e-5) mask[X3(ix,iy,iz)] = 0.0;
    }
    
    
    // LOAD PRIORS IN
    std::vector<image::basic_image<double,3> > priorsIn(iNrOfInsideClasses,image::basic_image<double,3>(input[0].geometry()));
    for(int i=0; i<iNrOfInsideClasses; i++)
        std::fill(priorsIn[i].begin(),priorsIn[i].end(),1.0/(double) iNrOfInsideClasses);
    
    // RENORMALIZE ATLAS WITHIN MASK
    for (int iz=0; iz< iNz; iz++)
        for (int iy=0; iy< iNy; iy++)
            for (int ix=0; ix< iNx; ix++)
            {
                if(mask[X3(ix,iy,iz)]==1.0){
                    double temp = 0.0;
                    for(int l=0; l<priorsIn.size(); l++) temp += priorsIn[l][X3(ix,iy,iz)];
                    for(int l=0; l<priorsIn.size(); l++) priorsIn[l][X3(ix,iy,iz)] /= temp;
                }
            }

if(!bDisjoint){
    
    //////////
    //DO
    //////////
    
    std::vector<image::basic_image<double,3> > posteriorIn(iNrOfInsideClasses,image::basic_image<double,3>(input[0].geometry()));
    std::vector<image::basic_image<double,3> > posteriorOut(iNrOfOutsideClasses,image::basic_image<double,3>(input[0].geometry()));
    
    image::basic_image<double,3> u(input[0].geometry());
    
    std::vector<std::vector<double> > pMeanIn(iNrOfInsideClasses,std::vector<double>(iNrOfChannels));
    std::vector<std::vector<double> > pSigmaIn(iNrOfInsideClasses,std::vector<double>(iNrOfChannels*iNrOfChannels));
    std::vector<std::vector<double> > pMeanOut(iNrOfOutsideClasses,std::vector<double>(iNrOfChannels));
    std::vector<std::vector<double> > pSigmaOut(iNrOfOutsideClasses,std::vector<double>(iNrOfChannels*iNrOfChannels));
    
    std::vector<double> pPriorIn(iNrOfInsideClasses);

    image::segmentation::BASE_joint(input, mask, init, u, priorsIn, priorsOut, posteriorIn, posteriorOut, pMeanIn, pSigmaIn, pMeanOut, pSigmaOut, pPriorIn,
                              
                              iNrOfOtherClasses,

                              pStdMin,
                              pMin,
                              pIOffset,
                              pIMax,
                              
                              iMaxIter,
                              iMaxIterEM,
                              iMaxIterSeg,
                              
                              pEpsEM,
                              pEpsSeg,
                              
                              pLambda,
                              pMu,
                                 
                              pInitAlpha,
                              bUpdateAlpha,
                                 
                              pCoupling,
                                 
                              bEMreinitialize,

                              basename,
                              bVerbose,
                                 
                              bInit,
                                    
                              iHyperHypoIntens_toWhich,
                              pHyperHypoIntens_howMuch,
                              iInclusion_toWhich);

    //////////
    //WRITE
    //////////

    // POSTERIORS
    for(int i=0; i<posteriorIn.size(); i++){
        std::stringstream ss; ss << i;
        image::serialization::nii_out<double,3>(posteriorIn[i],basename + "/classIn_" + ss.str()  + ".nii",nifti0);}

    for(int i=0; i<posteriorOut.size(); i++){
        std::stringstream ss; ss << i;
        image::serialization::nii_out<double,3>(posteriorOut[i],basename + "/classOut_" + ss.str()  + ".nii",nifti0);}

    // SEGS
    image::serialization::nii_out<double,3>(u,basename+"/seg.nii",nifti0);
    for (int iz=0; iz< iNz; iz++)
    for (int iy=0; iy< iNy; iy++)
    for (int ix=0; ix< iNx; ix++)
    {
        if(u[X3(ix,iy,iz)]>=0.5) u[X3(ix,iy,iz)]=1.0; else u[X3(ix,iy,iz)]=0.0;
    }
    image::serialization::nii_out<double,3>(u,basename+"/seg_bin.nii",nifti0);
    
}
else{
    
    //////////
    //DO
    //////////
    
    std::vector<std::vector<double> > lesionIndicator;
    int iNrOfLesionStates = image::fillLesionIndicator(lesionIndicator, iNrOfChannels);
    
    std::vector<std::vector<image::basic_image<double,3> > > posteriorIn(iNrOfLesionStates,std::vector<image::basic_image<double,3> >(iNrOfInsideClasses,image::basic_image<double,3>(input[0].geometry())));
    std::vector<std::vector<image::basic_image<double,3> > > posteriorOut(iNrOfLesionStates,std::vector<image::basic_image<double,3> >(iNrOfOutsideClasses,image::basic_image<double,3>(input[0].geometry())));
    
    std::vector<image::basic_image<double,3> > u(iNrOfChannels,image::basic_image<double,3>(input[0].geometry()));
    
    std::vector<std::vector<double> > pMeanIn(iNrOfInsideClasses,std::vector<double>(iNrOfChannels));
    std::vector<std::vector<double> > pSigmaIn(iNrOfInsideClasses,std::vector<double>(iNrOfChannels));
    std::vector<std::vector<double> > pMeanOut(iNrOfOutsideClasses,std::vector<double>(iNrOfChannels));
    std::vector<std::vector<double> > pSigmaOut(iNrOfOutsideClasses,std::vector<double>(iNrOfChannels));
    
    std::vector<std::vector<double> > pPriorIn(iNrOfInsideClasses,std::vector<double>(iNrOfChannels));
    
    image::segmentation::BASE_disjoint(input, mask, init, u, priorsIn, priorsOut, posteriorIn, posteriorOut, pMeanIn, pSigmaIn, pMeanOut, pSigmaOut, pPriorIn,
                                       
                                       iNrOfOtherClasses,
                                       
                                       pStdMin,
                                       pMin,
                                       pIOffset,
                                       pIMax,
                                       
                                       iMaxIter,
                                       iMaxIterEM,
                                       iMaxIterSeg,
                                       
                                       pEpsEM,
                                       pEpsSeg,
                                       
                                       pLambda,
                                       pMu,
                                       
                                       pInitAlpha,
                                       bUpdateAlpha,
                                       
                                       pCoupling,
                                       
                                       bEMreinitialize,
                                       
                                       basename,
                                       bVerbose,
                                       
                                       bInit,
                                       
                                       iHyperHypoIntens_toWhich,
                                       pHyperHypoIntens_howMuch,
                                       iInclusion_toWhich);
    
    
    //////////
    //WRITE
    //////////
    
    //SEGS
    for(int l=0; l<iNrOfChannels; l++){
        std::stringstream ss; ss << l;
        image::serialization::nii_out<double,3>(u[l],basename+"/seg_" + ss.str() + ".nii",nifti0);
    }
    for(int l=0; l<iNrOfChannels; l++){
        std::stringstream ss; ss << l;
        for (int iz=0; iz< iNz; iz++)
            for (int iy=0; iy< iNy; iy++)
                for (int ix=0; ix< iNx; ix++)
                {
                    if(u[l][X3(ix,iy,iz)]>=0.5) u[l][X3(ix,iy,iz)]=1.0; else u[l][X3(ix,iy,iz)]=0.0;
                }
        image::serialization::nii_out<double,3>(u[l],basename+"/seg_bin_"+ss.str()+".nii",nifti0);
    }    
}

    return 0;
}