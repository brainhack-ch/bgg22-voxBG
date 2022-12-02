function [v_o,n,f,v_rm] = hb_make_connected(v,conn)
% HB_MAKE_CONNECTED makes input bw array connected; the largest component
% will be kept and the remainder of pixels/voxels will be deleted.
% 
% Inputs:
%   v: bw 2D/3D array, or absolute address of a nifti, or spm_vol header. 
%
%   conn: 4 or 8 for 2D. 6, 18 or 26 for 3D; see bwconncomp.m for details.
%
% Outputs:
%   v_o: connected bw array. 
%
%   n: number of elements removed to make array connected.
%
%   f: fraction of voxels removed to make array connected.
%
%   v_rm: bw array showing removed elements: v == or(v_o,v_rm); 
%
% Hamid Behjat 

if ischar(v)
    v = hb_nii_load(v);
elseif isstruct(v)
    v = hb_nii_load(v.fname);
end

vdim = size(v);

CC = bwconncomp(v,conn);

[~,imax] = max(cellfun(@length,CC.PixelIdxList));

v_o = zeros(vdim);

v_o(CC.PixelIdxList{imax}) = 1;

v_rm = zeros(vdim);
for n=1:length(CC.PixelIdxList)
    if n==imax
        continue;
    end
    v_rm(CC.PixelIdxList{n}) = 1;
end

n = nnz(v_rm);
f = n/nnz(v);  
if f>0.05
    warning('A large number of elements were removed!'); 
end
end