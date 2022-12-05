import nibabel as nib

def save_gifti_as_obj(filename_gifti, save_name):
    gifti_file = nib.load(filename_gifti)
    f = open(save_name, "w")
    f.write("o surface_mesh\n")
    vertices = gifti_file.agg_data('NIFTI_INTENT_POINTSET')
    triangles = gifti_file.agg_data('NIFTI_INTENT_TRIANGLE')
    for i in range(0, vertices.shape[0]):
        f.write(" ".join(["v",str(vertices[i, 0]), str(vertices[i, 1]), str(vertices[i, 2])]) + "\n")
    for i in range(0, triangles.shape[0]):
        f.write(" ".join(["f",str(triangles[i, 0] + 1), str(triangles[i, 1] + 1), str(triangles[i, 2] +1 )]) + "\n")
    f.close()
