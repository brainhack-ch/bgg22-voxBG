function G = voxbg_adjacencymatrix_gmxx(G,opts)
% VOXBG_ADJACENCYMATRIX_GMXX Constructs gmlh/gmrh graph A matrix.
%
% Inputs:
%   G: graph structure based on HCP dataabse.   
%   
% Outputs:
%   G: updated G with fields A and indices. 
%
% Contributers: 
%    Hamid Behjat
%    Martin Larsson
% March 2017, April 2020. 

% Build A.
fprintf('\n..Determining adjacencies..')
A = voxbg_compute_adjacency(G.f.mask,26,'weight','no','sS',0);

% Prune A.
fprintf('\n..Pruning edges..');
%G = hb_prune_graph(A0,G,~opts.runPar);

f_mask = G.f.mask;

f_srfs = [
    G.f.surface.pial
    G.f.surface.white
    ];

info_srf = struct;
info_srf.name = {'pial','white'};
switch G.tissue
    case {'gmlh','gmrh'}
        info_srf.index = {1,2};
    case 'gm'
        info_srf.index = {[1 2],[3,4]};
end

outs = hb_prune_graph(...
    A,...
    f_mask,...
    f_srfs,...
    'info_srf',info_srf,...
    'parallelize',opts.parallel_pruning);%info_srf,~opts.runPar);

G.N = outs.N;
G.A = outs.A;
G.indices = outs.indices;
G.pruning = outs.pruning;
end
