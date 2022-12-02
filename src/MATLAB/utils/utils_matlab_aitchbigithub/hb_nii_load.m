function [v,h] = hb_nii_load(f,varargin)
%HB_NII_LOAD load nifti volume and header using external software. 

d = inputParser;
addParameter(d,'HeaderType','spm');
addParameter(d,'IndicesToLoad',[]);
parse(d,varargin{:});
opts = d.Results;

assert(ischar(f),'Enter absolute adress of nifti file.');

if contains(f,'.gz')
    gunzip(f);
    f = f(1:end-3);
    DELNONZIP = true;
else
    if exist(f,'file')
        DELNONZIP = false;
    else
        fgz = [f,'.gz'];
        if exist(fgz,'file')
            gunzip(fgz);
            DELNONZIP = true;
        else
            error('File missing: %s',f);
        end
    end
end

if isempty(opts.IndicesToLoad)
    FASTLOAD = false;
else
    FASTLOAD = true;
end
    
switch opts.HeaderType
    case 'spm'
        %assert(exist('spm_vol.m','file'),'"spm" package not in path.')
        h = spm_vol(f);
        Nv = length(h);
        if Nv==1
            v = spm_read_vols(h);
        else
            if FASTLOAD
                h1 = h(1);
                I = opts.IndicesToLoad;
                [x,y,z] = ind2sub(h1.dim,I);
                v1d0 = zeros(prod(h1.dim),1);
                v = zeros([h1.dim,Nv]);
                for iV=1:Nv
                    v1d = v1d0;
                    v1d(I) = spm_sample_vol(h(iV),x,y,z,0);
                    v(:,:,:,iV) = reshape(v1d,h1.dim);
                    showprgs(iV,Nv);
                end
            else
                v = spm_read_vols(h);
            end
        end
        if DELNONZIP
            delete(f);
        end
    otherwise
        error('extend.')
end
end

function showprgs(n,N)
l = numel(num2str(N));
if n==1
    fprintf('\n..Loading 4D nifti.. ');
else
    fprintf(repmat('\b',1,2*l+1),n);
end
eval(['fprintf(''%-',num2str(l),'d/%-',num2str(l),'d'',n,N)'])
end