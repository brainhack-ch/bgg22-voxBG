function ind = ml_parse_plane(plane)
% ML_PARSE_PLANE Convert char array to plane index
%	ind = ml_parse_plane(plane) returns a plane index [1,2,3] corresponding
%       to the inputs {1, 'x', 'sagittal'}, {2, 'y', 'coronal'} and 
%       {3, 'z', 'axial'}, respectively.
    
    switch plane
        case {1, 'x', 'sagittal'}
            ind = 1;
        case {2, 'y', 'coronal'}
            ind = 2;
        case {3, 'z', 'axial'}
            ind = 3;
        otherwise
            error('Unable to convert input to plane.');
    end
end

