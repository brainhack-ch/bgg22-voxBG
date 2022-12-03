import os
from pathlib import Path

import nibabel as nib


folder = Path(
    "/Users/scheltie/git/brainhack-ch/bgg22-voxBG/demos/data/HCP/dataset/"
)
folder_save = Path(
    "/Users/scheltie/git/brainhack-ch/bgg22-voxBG/demos/data/HCP/dataset/"
)
idx = 100307
gtype = "gmlh.res2000.spaceT1w"
tissue = gtype.split(".res")[0]

# skip wm for now and also_get_files_related_to_fMRI

hemisphere = "left" if "lh" in tissue else "right"
surftype = None
res = int(gtype.split(".res")[1].split(".")[0])
space = "T1w" if "T1w" in gtype else "Diffusion"

# creates some tags
resTag = ".res2000"
spaceTag = ".spaceT1w"
rsTag = ".res2000.spaceT1w"
trsTag = "gmlh.res2000spaceT1w"

settingsTag = ""
neighb = 3
maskTypeTag = f".neighb{neighb}"

# G.f structure with all folders
# hcp_root, hcpsave_root are equal to folder.
folder_t1w = folder / str(idx) / "T1w"
folder_t1w_save = folder_save / str(idx) / "T1w"
folder_t1w_results_save = folder_t1w_save / "Results"
folder_mni = folder / str(idx) / "MNINonLinear"
folder_mni_results = folder_mni / "Results"
folder_mni_save = folder_save / str(idx) / "MNINonLinear"
folder_mni_results_save = folder_mni_save / "results"
graphmain = folder_t1w_save / "graph"
graph = graphmain / gtype.replace(".", "_")

# create folders
os.makedirs(folder_t1w_save, exist_ok=True)
os.makedirs(folder_mni_save, exist_ok=True)
os.makedirs(folder_mni_results_save, exist_ok=True)
os.makedirs(graph, exist_ok=True)

# graph mat file
file_G = graph / f"G.{gtype}.mat"
n = f"{tissue}{maskTypeTag}.res{res}"

# volumetric files
source_file = folder_t1w / "ribbon.nii"
mask = graph / f"{n}.spaceT1w.nii"

# surface files
folder_surf = folder_t1w_save / "surf"  # d_surf
os.makedirs(folder_surf, exist_ok=True)
d = tissue[2:]  # lh
file_surf_pial = folder_surf / f"{d}.pial.surf.gii"
file_surf_white = folder_surf / f"{d}.white.surf.gii"
file_surfmesh_pial = folder_surf / f"{d}.pial"
file_surfmesh_white = folder_surf / f"{d}.white"
file_surfmesh_pial_ascii = folder_surf / f"{d}.pial.asc"
file_surfmesh_white_ascii = folder_surf / f"{d}.white.asc"

if not file_surf_pial.exists():
    if (folder / str(idx) / "preproc" / "dir").exists():
        d1 = folder / str(idx) / "preproc" / str(idx) / "T1w" / "surf"
    else:
        d1 = folder / str(idx) / "T1w" / str(idx) / "surf"
    d2 = file_surf_pial.stem
    d3 = file_surf_pial.suffix

    # copy file
    # copy(d1,[d2, d3]) to file_surf_pial

# same with white
if not file_surf_white.exists():
    if (folder / str(idx) / "preproc" / "dir").exists():
        d1 = folder / str(idx) / "preproc" / str(idx) / "T1w" / "surf"
    else:
        d1 = folder / str(idx) / "T1w" / str(idx) / "surf"
    d2 = file_surf_white.stem
    d3 = file_surf_white.suffix

    # copy file
    # copy(d1,[d2, d3]) to file_surf_white

# G.fname is equal to file_G
