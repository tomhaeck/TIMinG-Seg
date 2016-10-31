# Patient Nr1 from the MICCAI Ischemic Stroke Lesion Segmentation (ISLES) 2015 Challenge SPES data set [1].  Patient images were flipped along y-axis.

# STEP 1: BIAS FIELD CORRECTION WITH METHOD AND SOFTWARE OF YOUR CHOICE
# recommended:
#    ANTS N4BiasFieldCorrection [2]:
#    >> N4BiasFieldCorrection -i t1c.nii -o t1c_bias.nii -b [500,3]
#    >> N4BiasFieldCorrection -i t2.nii -o t2_bias.nii -b [500,3]


# STEP 2: REGISTRATION OF ATLAS (MOVING) TO PATIENT IMAGE (FIXED) WITH METHOD AND SOFTWARE OF YOUR CHOICE
# recommended:
#    ANTS:
#    >> ANTS 3 -m CC[t1c_bias.nii,../_atlas/T1.nii,1,4] -i 50x30x10 -o ab.nii
#    >> WarpImageMultiTransform 3 ../_atlas/WM.nii _WM_warped.nii -R t1c_bias.nii abWarp.nii abAffine.txt
#    >> WarpImageMultiTransform 3 ../_atlas/GM.nii _GM_warped.nii -R t1c_bias.nii abWarp.nii abAffine.txt
#    >> WarpImageMultiTransform 3 ../_atlas/CSF.nii _CSF_warped.nii -R t1c_bias.nii abWarp.nii abAffine.txt

# STEP 3: LESION SEGMENTATION WITH TIMinG-Seg
TIMinG-Seg -i t2_bias.nii -i ttp.nii -p _WM_warped.nii -p _GM_warped.nii -p _CSF_warped.nii -o timingseg.nii --mu 1e1 --alpha 0.5 --mask brainmask.nii --maxIterEM 1






# [1] O. Maier et al.. “ISLES 2015-A public evaluation benchmark for ischemic stroke lesion segmentation from multispectral MRI”. In: Media 35 (2017), pp. 250–269.
# [2] B. B. Avants, C. L. Epstein, M. Grossman, and J. C. Gee. “Symmetric diffeomorphic image registration with cross-correlation: evaluating automated labeling of elderly and neurodegenerative brain”. In: Media 12.1 (2008), pp. 26–41.
# [2] http://stnava.github.io/ANTs/ (on October 2016)