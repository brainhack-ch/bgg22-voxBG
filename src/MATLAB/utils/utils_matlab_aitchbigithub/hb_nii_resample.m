function f_o = hb_nii_resample(f_i,res_o,varargin)
%HB_NII_RESAMPLE Resamples an input nifti volume to a new resolution. 
% The new resampled volume is written to the directory of the input volume,
% unless name of output file is specified.
%
% Uses: spm_sample_vol.m, spm_slice_vol.m
%
% Inputs:
%   f_i: input file name, full address. Or spm_vol header.
%
%   res_o: resolution to resample f_i, in mm. Scalar value for isotropic
%   resolution voxels, or 3x1 vector for non-isotropic resolution.
%
%   Name-Value Pair Arguments: 
%   'InterOrder': interpolation order; see spm_sample_vol.m or
%   spm_slice_vol.m for details. [default: 1]
%
%   'OutputFile': output file name, absolute address.
%
%   'Method': resampling approach; 'approach1' that works in 3D [sample by
%   sample] or 'approach2' that works in 2D [slice by slice]. Results from
%   both approaches should be identical.
%
%   'InterpOrder': interpolation order. 
% 
%   'MemroySafe': logical; applicable to 'approach1'. 
%
% Outputs:
%   f_o: output file name.
%
% NOTES:
% approach1 slower than approach1, if memroySafe=true. In general,
% {'approach1',memorySafe=true}, {'approach1',memorySafe=false],
% {'approach2'} are all pretty fast. Computation time for instance for
% 0.7mm^3 to 1.25mm^3 around 0.3 secs.
%
% Example:
% f_o = hb_nii_resample(f_i,1.25);
%
% See also: hb_nii_reslice.m
%
% Hamid Behjat

d = inputParser;
addParameter(d,'Method','approach1');
addParameter(d,'InterpOrder',1);
addParameter(d,'MemorySafe',true);
addParameter(d,'OutputFile',[]);
parse(d,varargin{:});
opts = d.Results;

assert(ischar('f_i'),'f_i: input file absolute address.')

try
    h_i = spm_vol(f_i);
catch
    f_igz = [f_i,'.gz'];
    gunzip(f_igz);
    h_i = spm_vol(f_i);
end

d = abs(diag(h_i.mat));
res_i = [d(1);d(2);d(3)];

if length(res_o)==1
    res_o = repmat(res_o,3,1);
end

scl = res_o./res_i;

if any(scl<1)
    % approach1/2 can only handle downsampling.
    assert(strcmp(opts.Method,'approach3'),...
        'HB: mising method needed for upsampling vol.');
end

h_o_dim = round(h_i.dim(:)'./scl(:)');
xx = 1:h_o_dim(1);
yy = 1:h_o_dim(2);
zz = 1:h_o_dim(3);

if isempty(opts.OutputFile)
    if ~isequal(res_o(1),res_o(2),res_o(3))
        tag = '_hbResampled';
    else
        d = sprintf('%04d',res_o(1)*1e3);
        tag = ['_res',d];
    end
    [p_i,n_i] = fileparts(f_i);
    f_o = fullfile(p_i,[n_i,tag,'.nii']);
else
   f_o = opts.OutputFile;
end

h_o = struct();
h_o.fname = f_o;
h_o.dim = h_o_dim;
d = [res_o(:);1].*sign(diag(h_i.mat));
h_o.mat = h_i.mat-diag(diag(h_i.mat))+diag(d);
h_o.dt = h_i.dt;
if any(h_i.pinfo(1:2)~=[1;0])
    error('fishy. See NOTE 1.'); 
end
h_o.pinfo = h_i.pinfo;
d = dbstack;
h_o.descrip = [h_i.descrip,...
    ' - resampled with: ',...
    d(1).file,' -',...
    opts.Method];

h_o = spm_create_vol(h_o);

switch opts.Method
    case 'approach1' % 3D interpolation [sample by sample]
        A = h_i.mat\h_o.mat;
        if opts.MemorySafe
            for iZ=1:h_o.dim(3)
                [X,Y,Z] = ndgrid(xx,yy,zz(iZ));
                XYZ = A*[X(:),Y(:),Z(:),ones(length(X(:)),1)]'; % see NOTE 2.
                X = reshape(XYZ(1,:),h_o.dim(1),h_o.dim(2));
                Y = reshape(XYZ(2,:),h_o.dim(1),h_o.dim(2));
                Z = reshape(XYZ(3,:),h_o.dim(1),h_o.dim(2));
                d = spm_sample_vol(h_i,X,Y,Z,opts.InterpOrder);
                spm_write_plane(h_o,d,iZ);
            end
        else
            [X,Y,Z] = ndgrid(xx,yy,zz);
            XYZ = A*[X(:),Y(:),Z(:),ones(length(X(:)),1)]';
            X = reshape(XYZ(1,:),h_o.dim(1),h_o.dim(2),h_o.dim(3));
            Y = reshape(XYZ(2,:),h_o.dim(1),h_o.dim(2),h_o.dim(3));
            Z = reshape(XYZ(3,:),h_o.dim(1),h_o.dim(2),h_o.dim(3));
            v = zeros(size(X));
            for iZ=1:size(X,3)
                v(:,:,iZ) = spm_sample_vol(h_i,...
                    X(:,:,iZ),Y(:,:,iZ),Z(:,:,iZ),...
                    opts.InterpOrder);
            end
            spm_write_vol(h_o,v);
        end
    case 'approach2' % 2D interpolation [slice by slice]
        for j=1:h_o.dim(3)
            d = spm_slice_vol(h_i,...
                spm_matrix(...
                [0 0 j*scl(3) 0 0 0 scl(:)' 0 0 0]),...
                h_o.dim(1:2),opts.InterpOrder);
            spm_write_plane(h_o,d,j);
        end
        
end

f_o = h_o.fname;

% NOTE 1 ------------------------------------------------------------------
% pinfo[1] is the scaling factor. pinfo[2] is the offset value. By default,
% these should be 1 and 0, respectively. If they are not, I am not sure how
% this affect the resulting resampled volume. Should these values also be
% modified?

% NOTE 2 ------------------------------------------------------------------
% The output image voxel coordinates are first transformed to the milimeter
% space using h_o.mat, and then from milimeter space to voxel cooridinates
% of the input image. The resulting voxel coordinates, which will be non
% integers, will be positions at which the input image should be sampled
% at, and the resulting values will be assigned to the output image.
% Mathwise:
%
% Let, 
% 1) M_i   = h_i.mat
% 2) M_o   = h_o.mat
% 3) VC_o  = 4xN matrix, where the first 3 elements of each column gives
% the coordinate of a single voxel (integer values) in the output space,
% and the last element is 1.
% 4) VC_i  = 4xN matrix, giving the corresponding coordinates (not
% necessarily integers) in the input image; the input image should be
% resampled at these values to obtain the voxel values of the output image.
% VC_mm = 4xN matrix, giving the correspoding coordinates in mm space.
%
% Then, 
% VC_mm = M_o * VC_o;
% VC_i = inv(M_i)*VC_mm
% which is implemented as: VC_i = (M_i\M_o)*VC_o
% -------------------------------------------------------------------------
