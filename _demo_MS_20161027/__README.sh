#!/bin/bash

# Patient training01_01 from the ISBI Longitudinal MS Lesion Segmentation 2015 Challenge [1]

# STEP 1: BIAS FIELD CORRECTION WITH METHOD AND SOFTWARE OF YOUR CHOICE
# recommended:
#    ANTS N4BiasFieldCorrection [2]:
#    >> N4BiasFieldCorrection -i t1.nii -o t1_bias.nii
#    >> N4BiasFieldCorrection -i flair.nii -o flair_bias.nii

# STEP 2: REGISTRATION OF ATLAS (MOVING) TO PATIENT IMAGE (FIXED) WITH METHOD AND SOFTWARE OF YOUR CHOICE
# recommended:
#    ANTS:
#    >> ANTS 3 -m CC[flair_bias.nii,../_atlas/T1.nii,0.75,4] -m CC[t1_bias.nii,../_atlas/T1.nii,0.25,4] -i 50x30x10 -o ab.nii
#    >> WarpImageMultiTransform 3 ../_atlas/WM.nii _WM_warped.nii -R t1_bias.nii abWarp.nii abAffine.txt
#    >> WarpImageMultiTransform 3 ../_atlas/GM.nii _GM_warped.nii -R t1_bias.nii abWarp.nii abAffine.txt
#    >> WarpImageMultiTransform 3 ../_atlas/CSF.nii _CSF_warped.nii -R t1_bias.nii abWarp.nii abAffine.txt

# STEP 3: LESION SEGMENTATION WITH TIMinG-Seg
TIMinG-Seg -i flair_bias.nii -i t1_bias.nii -p _WM_warped.nii -p _GM_warped.nii -p _CSF_warped.nii -o timingseg.nii --mu 1e0 --alpha 0.5 --mask brainmask.nii --disjoint --hyp 1 --hyp 0 --hypamount 2  --hypamount -2 --inclusion 0 --maxIterEM 10



# Remark: In the case of MS segmentation, it may be recommended to mask out areas close to air-filled cavities, i.e. sinuses or ear canals,  before running the TIMinG-Seg algorithm, as these areas may cause false positive lesion segmentations.

# [1] http://iacl.ece.jhu.edu/index.php/MSChallenge (on October 2016)
# [2] B. B. Avants, C. L. Epstein, M. Grossman, and J. C. Gee. “Symmetric diffeomorphic image registration with cross-correlation: evaluating automated labeling of elderly and neurodegenerative brain”. In: Media 12.1 (2008), pp. 26–41.
# [2] http://stnava.github.io/ANTs/ (on October 2016)