function out = ml_get_volume(in)
% ML_GET_VOLUME Load NIFTI volume
%   out = ML_GET_VOLUME(in) is a convenience function. Input can be a NIFTI
%   filename (*.nifti), a NIFTI header structure or a double volume. In any
%   case, the output will be a double volume.

    if ischar(in)
        hb_gunzip(in);
        out = spm_read_vols(spm_vol(in));
    elseif isstruct(in)
        out = spm_read_vols(in);
    else
        out = in;
    end
end
