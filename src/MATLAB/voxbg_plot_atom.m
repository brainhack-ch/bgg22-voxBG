% Based on ml_plot_atom_box.m, extensively extended on multiple fronts,
% e.g.:
%1. Added option to plot atom of any given kernel, not just heat kernel. 
%2. Added option to plot Gaussian kernel of a given fwhm (in mm) to compare graph atom against. 
%3. Added option to input axis handle to plot in.
%4. Added option yo plot atoms in a format adapted to FSLEYES.
%5. Added option to plot atoms in form of contours instead of as imagesc.
%6. Made standalone by nesting required helping functions. 

%--------------------------------------------------------------------------
% ML_PLOT_ATOM_BOX Plot a slice of an atom within a bounding box
%   [axa,axu] = ML_PLOT_ATOM_BOX(G,node,plane) plots a slice of an atom
%       produced by the heat kernel. node specifies a linear index in the
%       nodes in G at which the atom is localized. plane specifies which
%       plane to slice ('sagittal', 'coronal', 'axial' or index of normal
%       axis). The outputs axa and axu specifies the axes for the atom and
%       underlay, respectively.
%   [axa,axu] = ML_PLOT_ATOM_BOX(__,name,value) specifies plotting
%       properties using one or more name-value pairs. The following
%       properites are available:
%
%       Alpha - the transparency of the atom. 'mask' indicates that voxels
%           outside of the graph should be completely transparent. 'linear'
%           sets the transparency at each pixel based on the corresponding
%           atom value.
%       BBWidth - the width of the (square) bounding box. Must be odd.
%       ChebOrd - the order of the Chebyshev polynomial used to approximate
%           the heat kernel.
%       Colormap - the colormap of the atom.
%       Normalize - whether or not the atom should be normalize to the
%           range [0,1].
%       Tau - parameter controlling the bandwidth of the heat kernel.
%       Underlay - NIFTI filename, header, or volume to use a an underlay.
%       UnderlayColor - the background color of the plot.
%       UnderlayColormap - the colormap of the underlay.
%--------------------------------------------------------------------------

function [axa,axu,abb] = voxbg_plot_atom(G,node,plane,varargin)

    p = inputParser;
    addRequired(p,'G');
    addRequired(p,'node');
    addRequired(p,'plane');
    addParameter(p,'Alpha','mask',@(x) ischar(x) || isscalar(x));
    addParameter(p,'BBWidth',13,@(x) isscalar(x) && x >= 1 && mod(x,2) == 1);
    addParameter(p,'ChebOrd',50);
    addParameter(p,'Colormap','parula');
    addParameter(p,'Normalize',true,@(x) isscalar(x) && islogical(x));
    addParameter(p,'Tau',7,@(x) isscalar(x) && x > 0);
    addParameter(p,'Underlay',[]);
    addParameter(p,'UnderlayColor','k');
    addParameter(p,'UnderlayColormap','gray');
    addParameter(p,'flipColormap',false);
    addParameter(p,'plotType','imagesc'); % 'imagesc', 'contour'
    addParameter(p,'kernel',[]);
    addParameter(p,'GaussianFWHM',[]);
    addParameter(p,'AlsoPlotMaskedGaussian',false);
    addParameter(p,'axis_atom',[]);
    addParameter(p,'axis_gaussian',[]);
    addParameter(p,'axis_gaussian_masked',[]);
    addParameter(p,'viewType','standard'); % 'standard', 'fsleyes'(Left/Anterior shown on right)
    
    parse(p,G,node,plane,varargin{:});
    opts = p.Results;
    
    plane = getplane(opts.plane);
    if isempty(opts.kernel)
        kernel = @(x) exp(-x*opts.Tau);
        heatKernel = true;
    else
        kernel = opts.kernel;
        heatKernel = false;
    end

    uv = getvol(opts.Underlay);
    
    if ~isempty(opts.Underlay)
        assert(isequal(size(uv),G.dim),'Underlay size is not matching graph.');
    end
    
    % Construct atom.
    d = zeros(G.N, 1);
    d(node) = 1;
    if isempty(opts.ChebOrd)
        opts.ChebOrd=50;
    end
    f = getatom(G.L,G.lmax,kernel,d,opts.ChebOrd);
    if heatKernel
        assert(abs(min(f)) < 1e-6*max(f),...
            'Expected non-negative values when using the heat kernel.');
        f(f < 0) = 0;
    end
    
    %caxis(limits);
    
    % Normalized atom.
    if opts.Normalize
        if heatKernel
            f = f./max(f(:));
        else
            f = f./max(abs(f(:)));
        end
    end
    
    atomv = zeros(G.dim);
    atomv(G.indices) = f;
    
    % Gaussian kernel.
    if ~isempty(opts.GaussianFWHM)
        gaussianv = zeros(G.dim);
        gaussianv(G.indices(node)) = 1; % delta
        
        vox_dim = sqrt(sum(G.mat(1:3,1:3).^2));
        sigma = (opts.GaussianFWHM ./ vox_dim) / sqrt(8*log(2));
        gaussianv = imgaussfilt3(gaussianv,sigma,'padding',0,'FilterDomain','spatial');
        
        if opts.AlsoPlotMaskedGaussian
            gaussianv_masked = zeros(G.dim);
            gaussianv_masked(G.indices) = gaussianv(G.indices);
        end
    end
    
    % Initiate figure if missing.
    if isempty(opts.axis_atom)
        % if axis for atom is given, we assume that other axis(es) (if
        % applicable) are also given, and vice versa. 
       
        if isempty(opts.GaussianFWHM)
            Nplots = 1;
        else
            if opts.AlsoPlotMaskedGaussian
                Nplots = 3;
            else
                Nplots = 2;
            end
        end
        
        figure;
        for iP=1:Nplots
            switch iP
                case 1
                    subplot(1,Nplots,1);
                    %axis;
                    axis image
                    opts.axis_atom = gca;
                case 2
                    subplot(1,Nplots,2);
                    axis image
                    opts.axis_gaussian = gca;
                case 3
                    subplot(1,Nplots,3);
                    axis image
                    opts.axis_gaussian_masked = gca;
            end
        end
    end
    
    % Plot atom.
    bb = floor((opts.BBWidth-1)/2)*[-1 1; -1 1; -1 1];
    [axa,ay,az] = ind2sub(size(atomv),G.indices(node));
    abb = {
        axa+(bb(1,1):bb(1,2))
        ay+(bb(2,1):bb(2,2))
        az+(bb(3,1):bb(3,2))
        };
    abb{plane} = round(mean(abb{plane}));
    atomslice = squeeze(atomv(abb{1},abb{2},abb{3}));
    if isempty(opts.GaussianFWHM)
        slices{1} = atomslice;
        plotAx{1} = opts.axis_atom;
    else
        gaussianslice = squeeze(gaussianv(abb{1},abb{2},abb{3}));
        slices{1} = atomslice;
        slices{2} = gaussianslice;
        plotAx{1} = opts.axis_atom;
        plotAx{2} = opts.axis_gaussian;
        if opts.AlsoPlotMaskedGaussian
            gaussianslice_masked = squeeze(gaussianv_masked(abb{1},abb{2},abb{3}));
            slices{3} = gaussianslice_masked;
            plotAx{3} = opts.axis_gaussian_masked;
        end
    end
    Nplot = length(slices);
    
    alphas = cell(1,Nplot);
    
    for iP=1:Nplot
        
        if strcmp(opts.Alpha,'mask')
            alphas{iP} = slices{iP} ~= 0;
        elseif strcmp(opts.Alpha,'linear')
            alphas{iP} = (slices{iP}-min(slices{iP}(:)))./(max(slices{iP}(:))-min(slices{iP}(:)));
        else
            alphas{iP} = opts.Alpha;
        end
        
        if isempty(uv)
            
            axu = [];
            switch opts.plotType
                case 'imagesc'
                    imagesc(plotAx{iP},slices{iP}','AlphaData',alphas{iP}');
                case 'contour'
                    contour(plotAx{iP},slices{iP}');
            end
            
            axa = plotAx{iP};
            axa.Color = opts.UnderlayColor;
            axis(axa,'image')
            %         axa.CLim = [0 max(f)];
            
            % symmetrize colormap
            if ~heatKernel
                if any(slices{iP}(:)<0) % at least one negative voxel
                    d = max(abs(min(slices{iP}(:))),max(slices{iP}(:)));
                    axa.CLim = [-d d];
                else
                    d = length(opts.Colormap)/2;
                    opts.Colormap = opts.Colormap(d+1:end,:);
                end
            end
            
            if opts.flipColormap
                switch opts.Colormap
                    case 'pink'
                        d = pink;
                    case 'gray'
                        d = gray;
                    otherwise
                        error('extend switch.')
                end
                colormap(axa,flipud(d));
            else
                colormap(axa,opts.Colormap);
            end
            xticks(plotAx{iP},[]);
            yticks(plotAx{iP},[]);
        else
            uslice = squeeze(uv(abb{1},abb{2},abb{3}));
            imagesc(uslice',[min(uv(:)) max(uv(:))]);
            axu = gca;
            axu.Color = opts.UnderlayColor;
            colormap(axu,opts.UnderlayColormap);
            axis image
            xticks([]);
            yticks([]);
            
            switch opts.viewType
                case 'standard'
                    if plane == 1
                        axu.XDir = 'normal';
                    else
                        axu.XDir = 'reverse';
                    end
                    axu.YDir = 'normal';
                case 'fsleyes'
                    axu.XDir = 'normal';
                    axu.YDir = 'normal';
            end
            
            axa = axes();
            imagesc(slices{iP}','AlphaData',alphas{iP}');
            %         axa.CLim = [0 max(f)];
            axa.Color = 'none';
            axis image
            xticks([]);
            yticks([]);
            
            % Attempt to keep the axes align. TODO: Make 100 % reliable...
            axa.Position = axu.Position;
            f = gcf;
            f.SizeChangedFcn = @(~,~) set(axa,'Position',get(axu,'Position'));
            
            if isa(axu.Parent,'matlab.graphics.layout.TiledChartLayout')
                axes(axu);
            end
        end
        
        switch opts.viewType
            case 'standard'
                % Flip X-axis if in plane.
                if plane == 1
                    axa.XDir = 'normal';
                else
                    axa.XDir = 'reverse';
                end
                axa.YDir = 'normal';
            case 'fsleyes'
                axa.XDir = 'normal';
                axa.YDir = 'normal';
        end
        axis(plotAx{iP},'equal');
    end
end

function ind = getplane(plane)
    switch plane
        case {1,'x','sagittal'}
            ind = 1;
        case {2,'y','coronal'}
            ind = 2;
        case {3,'z','axial'}
            ind = 3;
    end
end

function out = getvol(in)
    if ischar(in)
        hb_gunzip(in);
        out = spm_read_vols(spm_vol(in));
    elseif isstruct(in)
        out = spm_read_vols(in);
    else
        out = in;
    end
end

function atom = getatom(L,lmax,kernel,signal,cheb_ord)
    g = {kernel, @(x) sqrt(1-kernel(x).^2)}; % note 
    c = cell(1,length(g));
    for k = 1:length(g)
        c{k} = sgwt_cheby_coeff(g{k},cheb_ord,cheb_ord+1,[0,lmax]);
    end
    d = sgwt_cheby_op(signal(:),L,c,[0,lmax]);
    atom = d{1};
    % note: additional dummy kernel for sgwt_inverse to converge.
end
