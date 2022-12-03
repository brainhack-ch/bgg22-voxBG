from copy import deepcopy
from pathlib import Path

import numpy as np
from nilearn.image import load_img, resample_img, resample_to_img, new_img_like
from scipy import ndimage

folder_hcp = Path(
    "/Users/scheltie/git/brainhack-ch/bgg22-voxBG/demos/data/HCP/dataset/"
)
idx = 100307
fname = folder_hcp / str(idx) / "T1w" / "ribbon.nii.gz"
img = load_img(fname)

# select part of the image
value = 3
data = img.get_fdata()
idx = np.where(data != value)
# set values different from target to 0
data[idx] = 0
img_sel = new_img_like(img, data=data, copy_header=True)
del img

# down-sampling, from 0.7 mm to 2 mm. The diagonal from affine is retrieved and
# the diagonal values are swapped with the new resolution.
new_res = 2
target_affine = deepcopy(img_sel.affine)
diagonal = [elt * new_res for elt in np.sign(target_affine.diagonal()[:-1])]
np.fill_diagonal(target_affine, diagonal + [1])
img_resampled = resample_img(
    img_sel, target_affine=target_affine, interpolation="nearest"
)

# to match another image, resample_to_img can be used.
# resample_to_img(source_img, target_img)

# remove isolated voxels
label, nb_obj = ndimage.label(img_resampled.get_fdata())
# find the largest connected components
sizes = [np.where(label.flatten() == k)[0].size for k in range(1, nb_obj + 1)]
keep = np.argmax(sizes) + 1
# set to 0 voxels to remove
data = img_resampled.get_fdata()
data[np.where(label != keep)] = 0
img_resampled = new_img_like(img_resampled, data=data, copy_header=True)

# compute adjacency of the voxels selected and resampled
