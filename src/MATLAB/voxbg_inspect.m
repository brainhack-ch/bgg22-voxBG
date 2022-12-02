function voxbg_inspect(G,opts,nodes)
%VOXBG_INSPECT generates multiple plots for inspecting a voxBG graph
%around specified input nodes or nodes selected interactively.
%
% Inputs:
%   G: graph structure output from hb_get_G.m. If more than one graph is
%   to be input, G should be a cell array of structres.
%
%   opts: structure with multiple fields related to the plots that are to
%   be made. If G is a cell array, i.e., multiple graphs, then opts may
%   also be a cell array; if the same opts is to be used for both graphs,
%   then a single opts structure is sufficient to be input, thus no cell
%   array.
%
%   nodes: (optional, default: [])
%   A vector of node indices, i.e., values within [1,G.N]. These are
%   node(s) at which to inspect the input graph(s). If length(G)==2, then
%   nodes should be a two element cell array, each cell specifying the
%   nodes at which two inspect each graph. If nodes is not input, or input
%   as [], then a desired node can be interactively selected.
%
%   NOTE 1: If two graph are input, and if nodes is empty, both input
%   graphs hould be from the same subject since the node that will be
%   selected interactively is based on the same structural image (t1w, or
%   whatever other image indicated by opts.refImage). As such, if two opts
%   structures are input (which is ok), some othe fields should be
%   identical between the two opts, but not all. In particular, the
%   following fields should be the same:
%   -refImage
%   -plane
%   -slice
%
% Examples:
% (1) hb_inspect_graph(G,opts); % select a node interactively
% (2) hb_inspect_graph(G,opts,nodes); % input node(s)
% (3) hb_inspect_graph({G1,G2},opts); % same opts for both Gs, select a node
% (4) hb_inspect_graph({G1,G2},{opts1,opts2}); % G-specific opts, select a node
% (5) hb_inspect_graph({G1,G2},{opts1,opts2},{nodes1,nodes2}); % G-specific opts and nodes
%
% Hamid Behjat

if not(exist('nodes','var'))
    nodes = [];
end

%-Check inputs.
[G,opts,nodes] = checkInputs(G,opts,nodes);

%-Append inputs.
[G,opts] = appendInputs(G,opts);

%-Make plots.
if isempty(nodes) % interactively select node(s)
    
    I = permute(spm_read_vols(spm_vol(hb_gunzip(opts{1}.refImage))),[2 1 3]);
    
    slDir = zeros(1,3);
    slDir(opts{1}.plane) = 1;
    sv = sliceViewer(...
        I,...
        'SliceDirection',slDir,...
        'SliceNumber',opts{1}.slice,...
        'DisplayRangeInteraction','off');
    
    ax = gca;
    ax.XDir = 'reverse';
    ax.YDir = 'normal';
    set(ax.Children,'ButtonDownFcn',...
        @(src,evt) onButtonDownFcn(src,evt,sv,G,opts));
    
    title(['Click in image to show',...
        ' ODFs and graph around selected point',...
        '+ realization of an atom']);
    set(gcf,...
        'Position',[800,10,1000,1000]);

else % use given nodes
    buildplots(G,opts,nodes);
end
end

%==========================================================================
function [GG,OPTS,NODES] = checkInputs(GG,OPTS,NODES)

Ng = length(GG);

%assert(Ng<=2,'Either one or two graphs can be input.');

switch Ng
    case 1
        if not(iscell(GG))
            assert(isstruct(GG),...
                'G should be a structre: hb_get_G.m output.');
            GG = {GG};
        end
    otherwise
        assert(iscell(GG),'For multiple graphs, enter G as a cell array.');
end

if iscell(OPTS)
    assert(length(OPTS)==Ng,...
        ['opts can either be a single structure, ',...
        'or a cell array of structures of length GG']);
else
    assert(isstruct(OPTS),...
        ['opts should be a structre, ',...
        'or a cell array of structures of length GG']);
    d = OPTS;
    OPTS = cell(1,Ng);
    for iG=1:Ng
        OPTS{iG} = d;
    end
end

for iG=1:Ng
    switch OPTS{iG}.plane
        case {1, 'x', 'sagittal'}
            OPTS{iG}.plane = 1;
        case {2, 'y', 'coronal'}
            OPTS{iG}.plane = 2;
        case {3, 'z', 'axial'}
            OPTS{iG}.plane = 3;
    end
    
    % ref image
    switch OPTS{iG}.whichRefImage
        case 't1w'
            OPTS{iG}.refImage = GG{iG}.f.t1w_graphspace;
        case 'mask'
            OPTS{iG}.refImage = GG{iG}.f.mask;
    end
end

if isempty(NODES)
    
    if Ng>1
        % Given two graphs, for interactive selection of a ndoe, the code
        % enables to only select a single coordinate on a single image, thus,
        % we expect the specified opts.refImage on which coordinate is
        % selected, as wella s the opts.plane be the same for both graphs.
        % Also it i assumed that the G.dim field be the same for both
        % graphs since it used to extract the coordinates of the selected
        % pixel and the associated index of graph node for each graph.
        
        %assert(isequal(OPTS{1}.refImage,OPTS{2}.refImage),'opts.refImage should be the same for both graphs.');
        fprintf('\n..Assuming ref image the same for all graphs; using that of 1st graph.');
        disp(' ');
        disp(' ');
        assert(isequal(OPTS{1}.plane,OPTS{2}.plane),...
            'opts.plane should be the same for both graphs.');
        assert(isequal(OPTS{1}.slice,OPTS{2}.slice),...
            'opts.slice should be the same for both graphs.');
        assert(isequal(GG{1}.dim,GG{2}.dim),...
            'G.dim should be the same for both graphs.');
    end
else
    switch Ng
        case 1
            assert(isnumeric(NODES),...
                'A vector specifying node indices is expected.');
        otherwise
            d1 = iscell(NODES);
            d2 = length(NODES)==Ng;
            for iG=1:Ng
                d3 = isnumeric(NODES{iG});
            end
            assert(all([d1,d2,d3]),...
                ['A cell array of length equal to the number of graph, ',...
                'each cell being a vector ',...
                'specifying node indices is expected ']);
            
            % NOTE: same nodes cannot be specified for two graph, as it is
            % very lickely that the two graphs are represented by different
            % brain masks, even if they have the same nature and are from
            % the same subject; for e.g. the cerebrum-voxBG graph of the
            % same subject, when built with neighb3 and neigh5 option may
            % be represented by a slightly different mask since the voxels
            % that are removed post pruning differ for the two designs,
            % reulting in slightly varied masks, which in turn, completely
            % changes node indices for the two designs.
    end
end
end

%==========================================================================
function [G,opts] = appendInputs(G,opts)

Ng = length(G);

GG = G;
optss = opts;
clear G opts;

for iG=1:Ng
    
    G = GG{iG};
    opts = optss{iG};
    
    %-Append G.
    if not(isfield(G,'L'))
        G.L = sgwt_laplacian(G.A,'opt','normalized');
        G.lmax = sgwt_rough_lmax(G.L);
    end
    if not(isfield(G,'dim'))
        G.dim = spm_vol(hb_gunzip(G.f.mask)).dim;
    end
    
    %-Append opts.
    % command window print
    if not(isfield(opts,'PrintNodeIndex'))
        opts.PrintNodeIndex = false;
    end
    
    % figures
    if not(isfield(opts,'FigMenuBar'))
        %'figure' | 'none'
        % [set to 'none' to skip menu bar in the figures]
        opts.FigMenuBar = 'figure';
    end
    if not(isfield(opts,'FigNumberTitle'))
        opts.FigNumberTitle = 'off'; % 'on' | 'off'
    end
    
    % odfs
    if not(isfield(opts,'odfs'))
        d = load(G.f.odf.fib_mat,'odf').odf;
        if size(d,2)~=G.N
            d = d(:,G.pruning.ind_pre_pruning_A_remained_post_pruning);
            assert(size(d,2)==G.N);
        end
        opts.odfs = repmat(d,2,1);
        opts.vertices = load(G.f.odf.fib_mat,'odf_vertices').odf_vertices;
    end
    
    % underlay
    switch opts.whichUnderlay
        case 't1w'
            opts.underlay = G.f.t1w_graphspace;
        case 'none'
            opts.underlay = [];
        otherwise
            error('extend');
    end
    opts.underlay_bw = G.f.mask;
    
    % colormap
    switch opts.atom_colormap
        case 'RdBu'
            opts.atom_colormap = flipud(cbrewer('div','RdBu',101));
        case 'fslRender3'
            opts.atom_colormap = load('colormap_fsl_render3.mat').LUT;
        otherwise
            % good to go; potentially should have addpathed. 
    end

    % update
    GG{iG} = G;
    
    optss{iG} = opts;
end
clear G opts

G = GG;
opts = optss;
clear GG optss
end

%==========================================================================
function onButtonDownFcn(~,evt,sv,GG,OPTS)

Ng = length(GG);

if evt.Button == 1
    xy = round(evt.IntersectionPoint(1:2));
    
    switch OPTS{1}.plane
        case  1
            pos = [sv.SliceNumber xy];
        case 2
            pos = [xy(1) sv.SliceNumber xy(2)];
        case 3
            pos = [xy sv.SliceNumber];
    end
    
    NodeOk = true;
    NODES = cell(1,Ng);
    indinv = zeros(GG{1}.dim);
    for iG=1:Ng
        indinv(GG{iG}.indices) = 1:length(GG{iG}.indices);
        NODES{iG} = indinv(pos(1),pos(2),pos(3));
        if NODES{iG} == 0
            if Ng==1
                fprintf('No node at pressed location.\n');
            else
                fprintf('No node at pressed location for graph %d.\n',iG);
            end
            NodeOk = false;
        else
            if OPTS{iG}.PrintNodeIndex
                if Ng==1
                    fprintf('Pressed (%d,%d,%d) = node %d\n',...
                        pos(1),pos(2),pos(3),NODES{iG});
                else
                    fprintf('Pressed (%d,%d,%d) = node %d [graph %d]\n',...
                        pos(1),pos(2),pos(3),NODES{iG},iG);
                end
            end
        end
    end
    disp(' ');
    
    if NodeOk
        buildplots(GG,OPTS,NODES);
    end
end
end

%==========================================================================
function buildplots(GG,OPTS,NODES)

Ng = length(GG);

for iG=1:Ng
    
    G = GG{iG};
    
    opts = OPTS{iG};
    
    nodes = NODES{iG};

    Nn = length(nodes);
    
    bbw = opts.bbwidth;
    
    d_figs = opts.d_savefigs;
    
    FigNumberTitle = opts.FigNumberTitle;
    
    if Ng==1
        figtag = '';
    else
        figtag = sprintf('graph:%d, ',iG);
    end

    %-Plot reference image-------------------------------------------------
    if opts.plot_ref
        
        % all nodes in same opts.plane slice?
        AllNodesInSameSlice = true;
        for iN = 1:Nn
            node = nodes(iN);
            switch opts.plane
                case 1
                    error('extend');
                case 2
                    error('extend');
                case 3
                    [~,~,d] = ind2sub(G.dim,G.indices(node));
                    if iN==1
                        d0 = d;
                    else
                        if ~isequal(d0,d)
                            AllNodesInSameSlice = false;
                            break;
                        end
                    end
            end
        end
        refv = spm_read_vols(spm_vol(hb_gunzip(opts.refImage)));
        
        if AllNodesInSameSlice
            h_refim = figure;
            set(h_refim,...
                'Position',[10 1000 500 600],...
                'MenuBar',opts.FigMenuBar,...
                'NumberTitle',FigNumberTitle);
        else
            disp(' ');
            disp(' ');
            warning('[HB] Given nodes do not fall in same plane.');
            fprintf('A seperate ref image displayed for each node.\n');
        end
        
        if Nn>7
            warning(...
                'HB: only 7 color availbale for bounding-boxes.');
        end
        cs = lines();
        
        for iN = 1:Nn
            node = nodes(iN);
            
            [nx,ny,nz] = ind2sub(G.dim,G.indices(node));
            
            if AllNodesInSameSlice
                if iN==1
                    [Lx,Ly] = plotRefSlice(refv,nx,ny,nz,opts.plane); %#ok<ASGLU>
                    tag = sprintf('%06d',node);
                else
                    tag = sprintf('%s_%06d',tag,node);
                end
            else
                h_refim = figure;
                set(h_refim,...
                    'Position',[10 1000 500 600],...
                    'MenuBar',opts.FigMenuBar,...
                    'NumberTitle',FigNumberTitle);
                [Lx,Ly] = plotRefSlice(refv,nx,ny,nz,opts.plane); %#ok<ASGLU>
            end
            
            switch opts.plane
                case 1
                    error('extend');
                case 2
                    error('extend');
                case 3
                    d = bbw-1;
                    rectangle(...
                        'Position',[nx-(d)/2-0.5,ny-(d)/2-0.5,bbw,bbw],...
                        'EdgeColor',cs(iN,:),...
                        'LineWidth',3);
                    text(30,Lx-5*iN,sprintf('node %d',iN),...
                        'Color',cs(iN,:),...
                        'FontWeight','bold',...
                        'FontSize',15);
            end
            
            if opts.saveFigs
                if AllNodesInSameSlice
                    if iN==Nn
                        d = fullfile(d_figs,sprintf('reference_nodes%s',tag));
                        saveas(h_refim,d,'png');
                        set(h_refim,'renderer','Painters');
                        saveas(h_refim,d,'epsc');
                        close(h_refim);
                    end
                else
                    d = fullfile(d_figs,sprintf('reference_node%06d',node));
                    saveas(h_refim,d,'png');
                    set(h_refim,'renderer','Painters');
                    saveas(h_refim,d,'epsc');
                    close(h_refim);
                end
            end
        end
    end
    
    % Plot ODFs----------------------------------------------------------------
    if opts.plot_odf
        hfodf = cell(1,length(nodes));
        
        for iN = 1:Nn
            node = nodes(iN);
            
            hfodf{iN} = figure;
            if opts.largePlots
                d = [10 10 2*500 2*475];
            else
                d = [10 10 .6*500 .6*475];
            end
            set(hfodf{iN},...
                'Position',d,...
                'Name',sprintf('%snode:%d',figtag,iN),...
                'MenuBar',opts.FigMenuBar,...
                'NumberTitle',FigNumberTitle);
            clf(hfodf{iN});
            
            if opts.odfsWithUnderlay
                % with underlay.
                ml_plot_odfs_box(G,node,opts.plane,opts.odfs,opts.vertices,...
                    'BBWidth',bbw,...
                    'Underlay',opts.underlay_bw);
                
                if opts.saveFigs
                    print(hfodf{iN},fullfile(d_figs,...
                        sprintf('odf.node%06d_r300',node)),'-dpng','-r300');
                    print(hfodf{iN},fullfile(d_figs,...
                        sprintf('odf.node%06d_r600',node)),'-dpng','-r600');
                end
            else
                % without underlay.
                % [since the underlay apparently masks ODFs that are
                % perpendicular to the 2D slice]
                ml_plot_odfs_box(G,node,opts.plane,opts.odfs,opts.vertices,...
                    'BBWidth',opts.bbwidth,'Underlay',[]);
                
                if opts.saveFigs
                    print(hfodf{iN},fullfile(d_figs,...
                        sprintf('odf.node%06d_r300_noUnderlay',node)),...
                        '-dpng','-r300');
                    print(hfodf{iN},fullfile(d_figs,...
                        sprintf('odf.node%06d_r600_noUnderlay',node)),...
                        '-dpng','-r600');
                    close(hfodf{iN});
                end
            end
            axis equal
        end
    end
    
    % Plot graphs--------------------------------------------------------------
    if opts.plot_graph
        optsg = struct();
        optsg.EdgeColorbar = false;
        optsg.EdgeColormap = [0 0 0; 1 0 0; 1 1 0];
        optsg.EdgeCLim     = [0.8 1];
        optsg.EdgeAlphamap = [0.2 1];
        optsg.EdgeALim     = [0.8 1];
        optsg.EdgeWidthmap = [1 4];
        optsg.EdgeWLim     = [0.8 1];
        optsg.VertColor    = [0.25 0.25 0.25];
        optsg.VertAlpha    = 1;
        optsg.VertSize     = 3;
        optsg.Underlay     = opts.underlay_bw;
        
        hfgraph = cell(1,length(nodes));
        for iN = 1:Nn
            node = nodes(iN);
            [nx,ny,nz] = ind2sub(G.dim,G.indices(node));
            sliceee = nz;
            optsg.XLim = nx+[-(bbw-1)/2 (bbw-1)/2];
            optsg.YLim = ny+[-(bbw-1)/2 (bbw-1)/2];
            
            
            hfgraph{iN} = figure;
            set(hfgraph{iN},...
                'Position',[350 10 .6*500 .6*475],...
                'Name',sprintf('%snode:%d, neighb:%d',figtag,iN,G.neighb),...
                'MenuBar',opts.FigMenuBar,...
                'NumberTitle',FigNumberTitle);
            
            plot_graph_slice(G,opts.plane,sliceee,optsg);
            
            if opts.saveFigs
                set(hfgraph{iN},'renderer','Painters');
                saveas(hfgraph{iN},fullfile(d_figs,...
                    sprintf('graph.node%06d',node)),'epsc');
                close(hfgraph{iN});
            end
        end
        
        if opts.plotForColorbar
            hf_clb = figure;
            set(hf_clb,...
                'Position',[100 1000 300 300],...
                'MenuBar',opts.FigMenuBar,...
                'NumberTitle',FigNumberTitle);
            optsg.EdgeColorbar = true;
            plot_graph_slice(G,opts.plane,sliceee,optsg);
            
            if opts.saveFigs
                set(hf_clb,'renderer','Painters');
                saveas(hf_clb,fullfile(d_figs,'graph_colorbar'),'epsc');
                close(hf_clb);
            end
        end
        
        % Add colorbar to last graph plot.
        % NOTE: The colormap will affect the underlay as well.
        % colorbar(bax(3));
        % colormap(bax(3),interp1([0.8 0.9 1],opts.EdgeColormap,linspace(0.8,1)));
        % bax(3).CLim = opts.EdgeCLim;
        
    end
    
    % Plot atoms---------------------------------------------------------------
    if opts.plot_atom
        
        ABB = cell(1,length(nodes));
        
        %ploti = 1;
        hfs = cell(1,9);
        
        for iN = 1:Nn
            node = nodes(iN);
            if isempty(opts.GaussianFWHM)
                hfs{iN} = figure;
                if opts.largePlots
                    d = [10 10 2*500 2*475];
                else
                    d = [700 10 .6*500 .6*475];
                end
                if isempty(opts.kernel) && not(isempty(opts.tau))
                    n = sprintf('%snode:%d, neighb:%d, tau:%d',...
                        figtag,iN,G.neighb,opts.tau);
                else
                    n = sprintf('%snode:%d, neighb:%d',figtag,iN,G.neighb);
                end
                set(hfs{iN},...
                    'Position',d,...
                    'Name',n,...
                    'MenuBar',opts.FigMenuBar,...
                    'NumberTitle',FigNumberTitle);
                axis;
                opts.axis_atom = gca;
                opts.axis_gaussian = [];
                opts.axis_gaussian_masked = [];
            else
                if isfield(opts,'axis_atom')
                    assert(isfield(opts,'axis_gaussian'),...
                        ['If axis given for atom, ',...
                        'it shoud also be given for Gaussian kernel.']);
                    if opts.MaskedGaussian
                        assert(isfield(opts,'axis_gaussian_masked'),...
                            ['If axis given for atom, it shoud ',...
                            'also be given for masked Gaussian kernel.']);
                    end
                else
                    opts.axis_atom = [];
                    opts.axis_gaussian = [];
                    opts.axis_gaussian_masked = [];
                end
            end
            [~,~,ABB{iN}] = hb_plot_atom(G,node,opts.plane,...
                'BBWidth',bbw,...
                'Underlay',opts.underlay,...
                'Alpha',opts.Alpha,...
                'Tau',opts.tau,...
                'plotType',opts.plotType,...
                'Colormap',opts.atom_colormap,...
                'ChebOrd',opts.ChebOrd,...
                'kernel',opts.kernel,...
                'GaussianFWHM',opts.GaussianFWHM,...
                'AlsoPlotMaskedGaussian',opts.MaskedGaussian,...
                'flipColormap',opts.flipColormap,...
                'axis_atom',opts.axis_atom,...
                'axis_gaussian',opts.axis_gaussian,...
                'axis_gaussian_masked',opts.axis_gaussian_masked);
            axis equal
        end
        
        if opts.plotForColorbar
            % Add colorbar to last atom plot.
            % colorbar(axa);
            
            hfa_clb = figure;
            set(hfa_clb,...
                'Position',[10 10 300 300],...
                'MenuBar',opts.FigMenuBar,...
                'NumberTitle',FigNumberTitle);
            hb_plot_atom(G,node,opts.plane,...
                'BBWidth',bbw,...
                'Underlay',opts.underlay,...
                'Alpha',opts.Alpha,...
                'Tau',opts.tau,...
                'plotType',opts.plotType,...
                'Colormap',opts.atom_colormap,...
                'ChebOrd',opts.ChebOrd,...
                'kernel',opts.kernel,...
                'flipColormap',opts.flipColormap);
            colorbar
        end
        
        dir_atoms = fullfile(d_figs,['atoms_',opts.atom_colormap]);
        if not(exist(dir_atoms,'dir'))
            mkdir(dir_atoms);
        end
        
        %set(hfa_clb,'renderer','Painters');
        %saveas(hfa_clb,fullfile(dir_atoms,['atom_colorbar','_',aopts.Colormap]),'epsc');
        %close(hfa_clb);
        
        % *** Manually save the 4 atom figures: File> Save as... > *.png
        % [4 figures: 3 atoms + 1 colorbar]
        % Name them as 'atom*_colormap', where * is 1, 2 or 3 & save them in
        % related colormap atom folder, i.e., dir_atoms, see above.
        %
        % Why save like this? Because for some weired reason, when saving using
        % saveas or print, either as png or eps, the region outside mask gets shown
        % as white rather than black.
        
    end
    
    % Plot tissue mask(s)------------------------------------------------------
    if opts.plotTissueMasks
        for iN = 1:Nn
            node = nodes(iN);
            inode = G.indices(node);
            abb = ABB{iN};
            hf_mask = figure;
            set(hf_mask,...
                'Position',[700 430 300 300],...
                'Name',sprintf('%snode:%d, neighb:%d',figtag,iN,G.neighb),...
                'MenuBar',opts.FigMenuBar,...
                'NumberTitle',FigNumberTitle);
            switch G.tissue
                case 'cerebrum'
                    if iN==1
                        h_gm = spm_vol(G.f.mask_gm);
                        h_wm = spm_vol(G.f.mask_wm);
                        v_gm = spm_read_vols(h_gm);
                        v_wm = spm_read_vols(h_wm);
                    end
                    d1 = double(squeeze(v_gm(abb{1},abb{2},abb{3})))*0.5;
                    d2 = double(squeeze(v_wm(abb{1},abb{2},abb{3})));
                    chk = intersect(find(d1(:)),find(d2(:)));
                    assert(isempty(chk),'fishy: GM & WM masks intersect.');
                    v = d1+d2; % WM: 1 - GM: 0.5
                    if ismember(inode,find(v_gm))
                        nodeTissue = 'gm';
                    else
                        nodeTissue = 'wm';
                    end
                case {'wm','wmlh','wmrh','gm','gmlh','gmrh'}
                    v = spm_read_vols(spm_vol(G.f.mask));
                    v = double(squeeze(v(abb{1},abb{2},abb{3})));
                    switch G.tissue
                        case {'wm','wmlh','wmrh'}
                            nodeTissue = 'wm';
                        case {'gm','gmlh','gmrh'}
                            nodeTissue = 'gm';
                    end
                otherwise
                    error('extend.')
            end
            switch nodeTissue
                case 'wm'
                    cl = [0 0 1];
                case 'gm'
                    cl = [0 1 0];
            end
            switch opts.plane % note1, see below.
                case 1
                    error('extend.');
                case 2
                    error('extend.');
                case 3
                    v = rot90(v,2);
                    v = v';
            end
            imagesc(v); % if alpha needed, see hb_plot_atom.m: ,'AlphaData',?);
            colormap gray;
            hold on;
            d = bbw/2;
            rectangle('position',[d,d,1,1],...
                'EdgeColor','none',...
                'FaceColor',cl);
            axis image;
            xticks([]);
            yticks([]);
            caxis([0 1]); % note2, see below
            %
            % note1:
            % this is sloppy, and not necesarily correct for dofferent plane
            % settings; the way the atom is created and plotted (in
            % hb_plot_atom.m) is not informed by spm_vol(G.f.mask).dim, and
            % thus I am doing a messy fix here.
            %
            % note2: to ensure GM always shown as gray even if all voxel within
            % bounding-box are part of mask.
        end
    end
    
    % Plot subset of ODFs------------------------------------------------------
    % Ploting only a subset of ODFs within the 2D slice.
    if opts.plot_odf_subset
        roi = opts.whichODFroi;
        
        if isfield(opts,'pixels_toplot')
            pixels_toplot = opts.pixels_toplot;
            assert(size(pixels_toplot,1)==bbw,...
                'specified ROI does not match specified bbwidth.');
            assert(size(pixels_toplot,2)==bbw,...
                'specified ROI does not match specified bbwidth.');
        else
            switch roi
                case 1
                    pixels_toplot = [
                        0 0 0 0 0 0 0 0 0 0 0 0 0
                        0 0 0 0 0 0 0 0 0 0 0 0 0
                        0 0 0 0 0 0 1 0 1 0 0 0 0
                        0 0 0 0 0 1 1 1 1 1 0 0 0
                        0 0 0 0 0 0 1 2 1 0 0 0 0
                        0 0 0 0 0 1 1 1 1 1 0 0 0
                        0 0 0 0 0 0 1 0 1 0 0 0 0
                        0 0 0 0 0 0 0 0 0 0 0 0 0
                        0 0 0 0 0 0 0 0 0 0 0 0 0
                        0 0 0 0 0 0 0 0 0 0 0 0 0
                        0 0 0 0 0 0 0 0 0 0 0 0 0
                        0 0 0 0 0 0 0 0 0 0 0 0 0
                        0 0 0 0 0 0 0 0 0 0 0 0 0
                        ];
                    
                case 2
                    pixels_toplot = [
                        0 0 0 0 0 0 0 0 0 0 0 0 0
                        0 0 0 0 0 0 0 1 0 1 0 0 0
                        0 0 0 0 0 0 1 1 1 1 1 0 0
                        0 0 0 0 0 0 0 1 2 1 0 0 0
                        0 0 0 0 0 0 1 1 1 1 1 0 0
                        0 0 0 0 0 0 0 1 0 1 0 0 0
                        0 0 0 0 0 0 0 0 0 0 0 0 0
                        0 0 0 0 0 0 0 0 0 0 0 0 0
                        0 0 0 0 0 0 0 0 0 0 0 0 0
                        0 0 0 0 0 0 0 0 0 0 0 0 0
                        0 0 0 0 0 0 0 0 0 0 0 0 0
                        0 0 0 0 0 0 0 0 0 0 0 0 0
                        0 0 0 0 0 0 0 0 0 0 0 0 0
                        ];
            end
        end
        pixels_toplot = fliplr(flipud(pixels_toplot')); %#ok<FLUDLR>
        
        % x_plot: columns to plot 1xN
        % y_plot: matrix, specifying which rows to plot in each of the columns
        %xy_toplot = 4:8%[]; % if both x/y_toplot=[], plots all.
        
        whichNode = 1;
        node = nodes(whichNode);
        
        % with underlay.
        hf_subset = figure;
        set(hf_subset,...
            'Position',[100 1000 2*500 2*475],...
            'MenuBar',opts.FigMenuBar,...
            'NumberTitle',FigNumberTitle);
        ml_plot_odfs_box(G,node,opts.plane,opts.odfs,opts.vertices,...
            'BBWidth',bbw,...
            'Underlay',opts.underlay_bw,...
            'pixels_toplot',pixels_toplot);
        axis equal
        
        if opts.saveFigs
            n1 = sprintf('ODF_node%06d_subset_roi%d_r300',node,roi);
            n2 = sprintf('ODF_node%06d_subset_roi%d_r600',node,roi);
            print(hf_subset,fullfile(d_figs,n1),'-dpng','-r300');
            print(hf_subset,fullfile(d_figs,n2),'-dpng','-r600');
        end
        
        % without underlay.
        % [since the underlay apparently masks ODFs that are perpendicular to
        % the 2D slice]
        clf(hf_subset);
        ml_plot_odfs_box(G,node,opts.plane,opts.odfs,opts.vertices,...
            'BBWidth',bbw,...
            'pixels_toplot',pixels_toplot);
        axis equal
        
        if opts.saveFigs
            n1 = sprintf('ODF_node%06d_subset_roi%d_r300_noUnderlay',node,roi);
            n2 = sprintf('ODF_node%06d_subset_roi%d_r600_noUnderlay',node,roi);
            print(hf_subset,fullfile(d_figs,n1),'-dpng','-r300');
            print(hf_subset,fullfile(d_figs,n2),'-dpng','-r600');
            close(hf_subset);
        end
    end
end
end

%==========================================================================
function [Lx,Ly] = plotRefSlice(v,nx,ny,nz,plane) %#ok<INUSL>
switch plane
    case 1
        error('extend');
    case 2
        error('extend');
    case 3
        I = squeeze(v(:,:,nz));
        I = I';
        imagesc(I);
        [Lx,Ly] = size(I);
end
colormap(gca,'gray');
axis image
axis off
ax = gca;
ax.XDir = 'reverse';
ax.YDir = 'normal';
end

