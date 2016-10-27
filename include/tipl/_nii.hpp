#ifndef NII_HPP
#define NII_HPP

#include "_nifti.hpp"

namespace image{

namespace serialization{

//---------------------------------------------------------------------------------------------------------
//LOAD NII-FILES (single)

// create basic_image<pixel_type,dimension> from nii-file
template<typename pixel_type, unsigned int dimension>
void nii_in(basic_image<pixel_type,dimension>& I0, std::string filename, image::io::nifti& myniifile){
	myniifile.load_from_file(filename.c_str());
	myniifile >> I0;
}

//---------------------------------------------------------------------------------------------------------
//WRITE NII-FILES (single)

// write basic_image<pixel_type,dimension> to nii-file
template<typename pixel_type, unsigned int dimension>
void nii_out(const basic_image<pixel_type,dimension>& I0, std::string filename, image::io::nifti& myniifile){
	myniifile << I0;
	myniifile.save_to_file(filename.c_str());
}

}
}
#endif //NII_HPP
