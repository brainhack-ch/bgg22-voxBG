function ml_plot_odfs_box(G,node,plane,odfs,vertices,varargin)
% ML_PLOT_ODFS_BOX Plot slice of volume with ODFs
%   ML_PLOT_ODFS_BOX(G,node,plane,odfs,vertices) plots the ODFs surrounding
%       a node in a slice. G specifies the graph structure and node is a
%       linear index in the vertices. plane specifies which plane to slice
%       ('sagittal', 'coronal', 'axial' or index of normal axis). odfs
%       contains the ODF values corresponding to the direction in vertices.
%   ML_PLOT_ODFS_BOX(__,name,value) specifies plotting
%       properties using one or more name-value pairs. The following
%       properites are available:
%
%       BBWidth - the width of the (square) bounding box. Must be odd.
%       Underlay - NIFTI filename, header, or volume to use a an underlay.

    p = inputParser;
    addRequired(p,'G');
    addRequired(p,'node');
    addRequired(p,'plane');
    addParameter(p,'BBWidth',13,@(x) isscalar(x) && x >= 1 && mod(x,2) == 1);
    addParameter(p,'Underlay',[]);
    addParameter(p,'pixels_toplot',[]);
    
    parse(p,G,node,plane,varargin{:});
    opts = p.Results;
    
    plane = ml_parse_plane(opts.plane);
    bbwidth = opts.BBWidth;
    
    mask = zeros(G.dim);
    mask(G.indices) = 1;

    % Parameters for pre-processing.
    % TODO: Consider moving some of these to name-value pairs in opts.
    params.shDegree         = 8;
    params.minNorm          = true;
    params.maxNorm          = false;
    params.odfPow           = 1; 
    params.directColor      = false; % directional color-coding
    params.indices          = G.indices;
    params.mask             = mask;
    params.vertices         = vertices;
    params.odfs             = odfs;
    
    % Options for showODFs.
    hbopts.modifiedGrayscale     = true;
    hbopts.modifiedGrayscaleType = 2;
    hbopts.underlayShading       = 'flat'; % Underlay shading type: 'interp', 'flat', 'faceted'
    hbopts.plotOdfAtVoxelCenter  = true;
    hbopts.openFigure            = false;
    hbopts.pixels_toplot         = opts.pixels_toplot;
    
    % Define slice.
    bb = floor((bbwidth-1)/2)*[-1 1; -1 1; -1 1];
    [ax,ay,az] = ind2sub(G.dim,G.indices(node));
    abb = {
        ax+(bb(1,1):bb(1,2)+1) % TODO: +1 is a hack.
        ay+(bb(2,1):bb(2,2)+1) % TODO: +1 is a hack.
        az+(bb(3,1):bb(3,2)) % TODO: no +1 is a hack.
        };
    abb{plane} = round(mean(abb{plane}));
    assert(plane == 3,'The code about does only work for axial plane.'); % TODO: Allow other planes.

    if params.minNorm
        odfs = params.odfs-min(params.odfs);
    end

    if params.maxNorm
        odfs = odfs./max(odfs);
    end

    odfs = odfs.^params.odfPow;
    maxNormPostPower = true;
    scalePostMaxNormPostPower = 1;
    if params.odfPow>1&&maxNormPostPower
        odfs = odfs./max(odfs);
        odfs = odfs/scalePostMaxNormPostPower;
    end

    % Spherical harmonic (sh) representation of ODFs.
    anglesODF = cart2sph_phys(params.vertices');
    [~,pT] = makePT(params.shDegree,anglesODF);
    sh = pT*odfs;

    % Place sh coefficients in volume.
    Nb = size(sh,1);
    shinv = zeros([size(params.mask),Nb]);
    for k=1:Nb
        d = zeros(size(params.mask));
        d(params.indices) = sh(k,:);
        shinv(:,:,:,k) = d;
    end

    if isempty(opts.Underlay)
        uv = 0;
    else
        uv = ml_get_volume(opts.Underlay);
        uv = uv(abb{1},abb{2},abb{3});
        uv = (uv-min(uv(:)))./max(uv(:));
    end

    roiSize = [1 bbwidth; 1 bbwidth; 1 bbwidth];
    roiSize(plane,2) = 1;
    
    showODFs_hb(shinv(abb{1},abb{2},abb{3},:),roiSize,uv,...
        params.mask(abb{1},abb{2},abb{3}),60,1,1,params.directColor,[],...
        hbopts);
    
    if plane ~= 1
        set(gca,'XDir','reverse');
    end
    axis off
end

