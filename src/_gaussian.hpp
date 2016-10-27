#ifndef GAUSSIAN_HPP_INCLUDED
#define GAUSSIAN_HPP_INCLUDED

#include "gmm/gmm.h"
#include "tipl/_nii.hpp"

namespace image{
//----------------------------------------------------------------------------------------------------------------------------
// Evaluates a scalar image in a 1D gaussian
template<typename pixel_type, unsigned int dimension>
inline void evaluate_in_gaussian(const basic_image<pixel_type,dimension>& src, basic_image<pixel_type,dimension>& dest, pixel_type mu, pixel_type std)
{
    typename basic_image<pixel_type,dimension>::const_iterator iter = src.begin();
    typename basic_image<pixel_type,dimension>::const_iterator end = src.end();
    typename basic_image<pixel_type,dimension>::iterator out = dest.begin();
    
    if(std>std::numeric_limits<pixel_type>::epsilon()){
        pixel_type twoStdSquaredRecip = 1.0 / (2.0 * std * std);
        pixel_type sqrtTwoPiTimesStdRecip = 1.0 / (sqrt(2.0 * M_PI) * std);
        
        while(iter!=end){
            pixel_type x = *iter++-mu;
            x *= x;
            x = (pixel_type) sqrtTwoPiTimesStdRecip * std::exp(-x * twoStdSquaredRecip);
            if(x>1.0) x = 1.0;
            *out++ = x;
        }
    }
    else{
        while(iter!=end){
            if(*iter++==mu) *out++ = 1.0;
            else *out++ = 0.0;
        }
    }
    
}
//----------------------------------------------------------------------------------------------------------------------------
// Evaluates a vector valued image in a multidimensional gaussian
template<typename pixel_type, unsigned int dimension>
void evaluate_in_gaussian_mult(const std::vector<basic_image<pixel_type,dimension> >& src, basic_image<pixel_type,dimension>& dest, std::vector<pixel_type>& mu, gmm::dense_matrix<pixel_type>& covar)
{
    geometry<dimension> geo = src[0].geometry();
    
    std::vector<pixel_type> y(src.size()), prod(src.size());
    typename std::vector<pixel_type>::iterator iter;
    typename std::vector<pixel_type>::const_iterator end = y.end();
    typename std::vector<pixel_type>::const_iterator itermu;
    typename std::vector<basic_image<pixel_type,dimension> >::const_iterator iterim;
    typename basic_image<pixel_type,dimension>::iterator out = dest.begin();
    
    pixel_type det = lu_det(covar);
    if(std::abs(det)>std::numeric_limits<pixel_type>::epsilon()){
        gmm::dense_matrix<pixel_type> covarinv(src.size(),src.size());
        gmm::copy(covar,covarinv);
        gmm::lu_inverse(covarinv);
        pixel_type sqrtTwoPiTimesStdRecip = 1.0 / std::sqrt(std::pow(2.0*M_PI,(pixel_type) src.size())*det);
        for(pixel_index<dimension> index; index.valid(geo); index.next(geo)){
            iter = y.begin();
            iterim = src.begin();
            itermu = mu.begin();
            while(iter!=end){
                *iter++ = (*iterim++)[index.index()]-(*itermu++);
            }
            gmm::mult(covarinv,y,prod);
            pixel_type x = (pixel_type) sqrtTwoPiTimesStdRecip * std::exp(- 0.5 * gmm::vect_sp(prod,y));
            if(x>1.0) x = 1.0;
            *out++ = x;
        }
    }
    else{
        for(pixel_index<dimension> index; index.valid(geo); index.next(geo)){
            iter = y.begin();
            iterim = src.begin();
            itermu = mu.begin();
            pixel_type sum = 0.0;
            while(iter!=end){
                sum += (*iterim++)[index.index()]-(*itermu++);
                iter++;
            }
            if(sum==0)
                *out++ = 1.0;
            else
                *out++ = 0.0;	    	
        }
    }
}
//----------------------------------------------------------------------------------------------------------------------------
}
#endif //GAUSSIAN_HPP_INCLUDED