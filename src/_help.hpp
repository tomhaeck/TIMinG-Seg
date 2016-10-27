#ifndef HELP_HPP_INCLUDED
#define HELP_HPP_INCLUDED

namespace image{
//----------------------------------------------------------------------------------------------------------------------------
//normalize a vector of images
template<typename pixel_type, unsigned int dimension>
void normalizeV(std::vector<basic_image<pixel_type,dimension> >& src, basic_image<pixel_type,dimension>& brainmask, pixel_type Imax){
    
    int iNrOfChannels = src.size();
    for(int i=0; i<iNrOfChannels;i++){
        pixel_type min_el = *std::min_element(src[i].begin(),src[i].end());
        src[i]-=min_el;
        image::normalize(src[i], Imax);
    }
}
//----------------------------------------------------------------------------------------------------------------------------
//create other classes
template<typename pixel_type, unsigned int dimension>
void create_other_classes(std::vector<basic_image<pixel_type,dimension> >& priorOut, basic_image<pixel_type,dimension>& brainmask, int iNrOfOtherClasses){
    
    geometry<dimension> geo = priorOut[0].geometry();
    basic_image<pixel_type,dimension> im(geo);
    
    if(iNrOfOtherClasses>0){
        for(image::pixel_index<dimension> index; index.valid(geo); index.next(geo)){
            pixel_type tmp = 1.0;
            for(int i=0; i<priorOut.size(); i++)
                tmp -= priorOut[i][index.index()];
                tmp /= (pixel_type) iNrOfOtherClasses;
                if(tmp>=0.0 && brainmask[index.index()])
                    im[index.index()] = tmp;
                    else
                        im[index.index()] = 0.0;
                        }
    for(int i=0; i<iNrOfOtherClasses; i++)
        priorOut.push_back(im);
    }
}
//----------------------------------------------------------------------------------------------------------------------------
//lesion states (i.e. binary counter)
template<typename pixel_type>
int fillLesionIndicator(std::vector<std::vector<pixel_type> >& lesionIndicator, int iNrOfChannels)
{
    int iNrOfLesionStates = std::pow((pixel_type) 2, (pixel_type) iNrOfChannels);
    lesionIndicator.resize(iNrOfLesionStates);
    
    for(int t=0; t<iNrOfLesionStates;t++){
        
        lesionIndicator[t].resize(iNrOfChannels);
        
        int dividend = t;
        int bin_num;
        int counter = iNrOfChannels-1;
        while ( dividend >= 1 )
        {
            bin_num = dividend % 2;
            lesionIndicator[t][counter]=bin_num;
            dividend /= 2;
            counter--;
        }
    }
    return iNrOfLesionStates;
}
}
#endif //HELP_HPP_INCLUDED