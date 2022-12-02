% [06.04.2022] extended "switch plane" to allow numerical plane
% identifiers.
%
function plot_graph_slice(G,plane,slice,varargin)
% PLOT_GRAPH_SLICE Plot a slice of a brain graph
%   PLOT_GRAPH_SLICE(G,plane,slice) plot a slice of the graph G in the
%       given plane ('sagittal', 'coronal', 'axial').
%   PLOT_GRAPH_SLICE(__,name,value) specifies plotting properties using one
%       or more name-value pairs. The following properites are available:
%
%       XLim, YLim - the plotting limits with the slice.
%       EdgeThreshold - edges with weights smaller than this a not plotted.
%       EdgeColor - 1 x 3 vector specifying the color of the edges or [] if
%           no edges should be drawn.
%       EdgeColorbar - whether to show a colorbar
%       EdgeColormap - m x 3 array of colors.
%       EdgeCLim - the weight limits within the colormap is applied.
%       EdgeAlpha - the transparency of the edges.
%       EdgeAlphamap - vector of alpha values.
%       EdgeALim - the weight limits within the alphamap is applied.
%       EdgeWidth - the line width used to draw the edges.
%       EdgeWidthmap - vector of values specifying the edge widths.
%       EdgeWLim - the weight limits within the widthmap is applied.
%       VertColor - 1 x 3 vector specifying the color of the vertices or []
%           if no vertices should be drawn.
%       VertEdgeColor, VertFaceColor - 1 x 3 vectors specifying the edge
%           and face color of the veritices.
%       VertAlpha - the transparency of the vertices.
%       VertSize - the marker size of the vertices.
%       VertLineWidth - the line width of the edge color for the vertices.
%       Underlay - a filename specifying a NIFTI file to use as underlay.
%       UnderlayAlpha - the transparency of the underlay. By setting the
%           Color property of the current axis to, e.g., white or black the
%           underlay can be made brighter or darker.

% Set default parameters.
switch plane
    case {'sagittal',1}
        defaultXLim = [1 G.dim(2)];
        defaultYLim = [1 G.dim(3)];
    case {'coronal',2}
        defaultXLim = [1 G.dim(1)];
        defaultYLim = [1 G.dim(3)];
    case {'axial',3}
        defaultXLim = [1 G.dim(1)];
        defaultYLim = [1 G.dim(2)];
end
defaultCAWLim = [min(nonzeros(G.A)) max(nonzeros(G.A))];

% Define validators.
isColor = @(x) isempty(x) || isequal(size(x),[1 3]);
isColormap = @(x) ismatrix(x) && (isempty(x) || size(x,2) == 3); % TODO: Allow function handles?
isMap = @(x) isvector(x);
isLim = @(x) isequal(size(x),[1,2]) && x(1) <= x(2);
isLog = @(x) islogical(x);
isXLim = @(x) isLim(x) && x(1) >= 1 && x(2) <= defaultXLim(2);
isYLim = @(x) isLim(x) && x(1) >= 1 && x(2) <= defaultYLim(2);

% Parse input.
p = inputParser;
addRequired(p,'G');
addRequired(p,'plane');
addRequired(p,'slice');
addParameter(p,'XLim',defaultXLim,isXLim);
addParameter(p,'YLim',defaultYLim,isYLim);
addParameter(p,'EdgeThreshold',0);
addParameter(p,'EdgeColor',zeros(1,3),isColor);
addParameter(p,'EdgeColorbar',false,isLog);
addParameter(p,'EdgeColormap',[],isColormap);
addParameter(p,'EdgeCLim',defaultCAWLim,isLim);
addParameter(p,'EdgeAlpha',1);
addParameter(p,'EdgeAlphamap',[],isMap);
addParameter(p,'EdgeALim',defaultCAWLim,isLim);
addParameter(p,'EdgeWidth',1);
addParameter(p,'EdgeWidthmap',[],isMap);
addParameter(p,'EdgeWLim',defaultCAWLim,isLim);
addParameter(p,'VertColor',zeros(1,3),isColor);
addParameter(p,'VertEdgeColor',[],isColor);
addParameter(p,'VertFaceColor',[],isColor);
addParameter(p,'VertAlpha',1);
addParameter(p,'VertSize',3);
addParameter(p,'VertLineWidth',1);
addParameter(p,'Underlay','');
addParameter(p,'UnderlayAlpha',1);

parse(p,G,plane,slice,varargin{:});
opts = p.Results;

opts.EdgeAlphamap = opts.EdgeAlphamap(:);
opts.EdgeWidthmap = opts.EdgeWidthmap(:);
if isempty(opts.VertEdgeColor)
    opts.VertEdgeColor = opts.VertColor;
end
if isempty(opts.VertFaceColor)
    opts.VertFaceColor = opts.VertColor;
end

% Plot underlay.
if ~isempty(opts.Underlay)
    uVol = spm_read_vols(spm_vol(hb_gunzip(opts.Underlay)));
    
    assert(isequal(size(uVol),G.dim),...
        'Underlay dimensions does not match graph');
    
    switch plane
        case {'sagittal',1}
            uSlice = squeeze(uVol(slice,:,:));
        case {'coronal',2}
            uSlice = squeeze(uVol(:,slice,:));
        case {'axial',3}
            uSlice = squeeze(uVol(:,:,slice));
    end
    xrange = opts.XLim(1):opts.XLim(2);
    yrange = opts.YLim(1):opts.YLim(2);
    uSlice = uSlice(xrange,yrange);
    imagesc(xrange,yrange,uSlice','AlphaData',opts.UnderlayAlpha);
    ax = gca;
    colormap(ax,'gray');
%     colormap(ax,[0.5 0.5 0.5; 1 1 1]);
    clims = [min(uVol(:)) max(uVol(:))];
    ax.CLim = clims;
    freezeColors_da(ax)
    hold on
end

% Find nodes in slice.
sliceVol = false(G.dim);
maskVol = false(G.dim);
maskVol(G.indices) = true;

switch plane
    case {'sagittal',1}
        maskSlice = squeeze(maskVol(slice,:,:));
        sliceVol(slice,:,:) = maskSlice;
        [x,y] = ind2sub(G.dim([2 3]),find(maskSlice));
    case {'coronal',2}
        maskSlice = squeeze(maskVol(:,slice,:));
        sliceVol(:,slice,:) = maskSlice;
        [x,y] = ind2sub(G.dim([1 3]),find(maskSlice));
    case {'axial',3}
        maskSlice = maskVol(:,:,slice);
        sliceVol(:,:,slice) = maskSlice;
        [x,y] = ind2sub(G.dim([1 2]),find(maskSlice));
end

% 2D adjacency.
sliceInd = ismember(G.indices,find(sliceVol));
sliceA = G.A(sliceInd,sliceInd);

[fromI,toI] = find(triu(sliceA));
weights = nonzeros(triu(sliceA));

% Sort edges by weight with the strong ones on top.
[weights, Isort] = sort(weights);
fromI = fromI(Isort);
toI = toI(Isort);

% Plot edges.
if (~isempty(opts.EdgeColor) || ~isempty(opts.EdgeColormap)) && opts.EdgeAlpha ~= 0
    minW = min(weights);
    maxW = max(weights);

    [cMap,cX] = prepMap(opts.EdgeColormap,opts.EdgeCLim,minW,maxW);
    [aMap,aX] = prepMap(opts.EdgeAlphamap,opts.EdgeALim,minW,maxW);
    [wMap,wX] = prepMap(opts.EdgeWidthmap,opts.EdgeWLim,minW,maxW);

    for i = 1:length(fromI)

        % Get edge start- and end-points.
        xLine = [x(fromI(i)); x(toI(i))];
        yLine = [y(fromI(i)); y(toI(i))];

        % Do not plot if edge is completely outside plotting area.
        if all((xLine > opts.XLim(2)) | (xLine < opts.XLim(1))...
                | (yLine > opts.YLim(2)) | (yLine < opts.YLim(1)))
            continue;
        end

        % Get edge weight.
        weight = weights(i);

        % Don't plot if edge weight is under threshold.
        if weight < opts.EdgeThreshold
            continue
        end

        % Set color.
        if isempty(cMap)
            color = opts.EdgeColor;
        else
            color = interp1(cX,cMap,weight);
        end

        % Set edge transparency.
        if isempty(aMap)
            color(4) = opts.EdgeAlpha;
        else
            color(4) = interp1(aX,aMap,weight);
        end

        % Set edge width.
        if isempty(wMap)
            lineWidth = opts.EdgeWidth;
        else
            lineWidth = interp1(wX,wMap,weight);
        end

        % Plot line.
        plot(xLine,yLine,'Color',color,'LineWidth',lineWidth);
        hold on
    end
end

if opts.EdgeColorbar
    colormap(ax,interp1(cX,cMap,linspace(opts.EdgeCLim(1),opts.EdgeCLim(2),200)))
    ax.CLim = opts.EdgeCLim;
    colorbar
end

% Plot vertices.
if opts.VertAlpha > 0 && ~isempty(opts.VertEdgeColor) && ~isempty(opts.VertFaceColor)
    Ivert = find(x >= opts.XLim(1) & x <= opts.XLim(2)...
        & y >= opts.YLim(1) & y <= opts.YLim(2));
    plot(x(Ivert),y(Ivert),'o',...
        'MarkerSize',opts.VertSize,...
        'MarkerEdgeColor',opts.VertEdgeColor,...
        'MarkerFaceColor',opts.VertFaceColor,...
        'LineWidth',opts.VertLineWidth);
end

% Set axis properties.
axis image
hold off
ax = gca();
ax.XTick = [];
ax.YTick = [];
ax.XLim = opts.XLim+[-0.5 0.5];
ax.YLim = opts.YLim+[-0.5 0.5];
    
% TODO: These transformations should ideally be contained within the graph
% structure G and/or be read from the underlay NIFTI header.
switch plane
    case {'sagittal',1}
        ax.XDir = 'normal';
        ax.YDir = 'normal';
    case {'coronal',2}
        ax.XDir = 'reverse';
        ax.YDir = 'normal';
    case {'axial',3}
        ax.XDir = 'reverse';
        ax.YDir = 'normal';
end
end

function [map,x] = prepMap(map,lim,minWeight,maxWeight)
x = [];
if ~isempty(map)
    x = linspace(lim(1),lim(2),size(map,1));
    if minWeight < lim(1)
        x = [minWeight x];
        map = [map(1,:); map];
    end
    if maxWeight > lim(2)
        x = [x maxWeight];
        map = [map; map(end,:)];
    end
end
end
