function outs = hb_prune_graph(A,f_mask,f_srfs,varargin)
%HB_PRUNE_GRAPH removes edges that pass through given surface files (pial
%&/ white surfaces) and consequently removes nodes that become disconnected
%after pruning.
%
% Inputs:
%   A: Graph adjacency matrix to be pruned. 
%
%   f_mask: bw mask specifying the graph vertices in the associated
%   image/volume.
%
%   f_srfs: 1xS cell array of strings, each string specifying the absolute
%   address of a surface file (.gii).
%
%   Name-Value pair arguments:
%
%   'info_srf': (default: []) structure with fields:
%       -name: cell array of stings, specifying the surface types;
%       supported names are: 'pial' and 'white'.
%       -index: 1xS cell array 1x? integer vectors, where
%       S=length(info_srf.name) and ? is the number of surface files for
%       each surface type.
% 
%   For example, info_srf.name = {'pial','white'} info_srf.index = {[1
%   2],[3 4]}, which means that files 1 and 2 in f_srfs are pial surfaces
%   and files 3 and 4 are white surfaces. Note that in this setting,
%   pruning will first be done using the two pial surfaces, and then, a
%   second level of pruning is performed using the white surfaces.
%
%   'parallelize': (default: false) use parallel pool to conduct pruning?
%
%   'EdgeTypeToRecover': (default: 'none'). Either a single string, or a
%   cell array of strings. This option allows recovering a subset of the
%   suprious edges that have been pruned; the default is to not recover any
%   pruned edge, i.e., an edge that passes through teh given surface(s),
%   but alternative options are in the form 'x-y', where x,y<=9, e.g.:
%
%       '1-1': edges that connects two class1 voxels.
%       '2-2': edges that connects two class2 voxels.
%       '1-2': edges that connects class1 to class2 voxels.
%
%   If EdgeTypeToRecover=~'none', NodeClass should aslo be input.
%
%   'NodeClass': (optional, but required if opts.EdgeTypeToRecover is
%   anything other than 'none'. A vector of size(A,1), with integers within
%   range [1,9], specifying the class of each graph vertex (voxel).
%
% Outputs:
%   outs: struture with a buncg of fields:
%       -A:
%       -indices:
%       -pruning: a structre with a buncg of fileds that provided deatils
%       about the outcomes of pruning, indcluding a .txt file using which
%       the prunned edges can be visulaized.
%
% Examples:
%   1. outs = hb_prune_graph(A,f_m,f_s,'parallelize',1);   
%   2. info_srf = struct;
%      info_srf.name = {'pial','white'};
%      info_srf.index = {[1,2],[3,4]};
%      outs = hb_prune_graph(A,f_m,f_s,'parallelize',1,'info_srf',info_s);
%   3. NC = ones(size(A,1),1);
%      NC(1:floor(size(A,1)/3)) = 2; 
%      outs = hb_prune_graph(A,f_m,f_s,'parallelize',1,...
%          'info_srf',info_s,...
%          'EdgeTypeToRecover',{'2-2','1-2'},...
%          'NodeClass',NC);
% 
% Hamid Behjat 

p = inputParser;
addParameter(p,'parallelize',false);
addParameter(p,'info_srf',[],@(x) isstruct(x) || isempty(x));
addParameter(p,'EdgeTypeToRecover','none');
addParameter(p,'NodeClass',[]);

parse(p,varargin{:});
opts = p.Results;

if isempty(opts.info_srf)
    Nsrf = 1;
    info_srf = struct;
    info_srf.name = {'unnamed'};
    info_srf.index = {1:length(f_srfs)};
    CompatibleSurfOrder = true;
else
    % If both pial and white file given, pial files should be given first,
    % to get a correct info file; but the pruning will work correctly, even
    % if it is desired to first do pruning using white surfaces and then
    % pial surfaces; it's just the code for the info file that needs to be
    % extended.
    info_srf = opts.info_srf;
    Nsrf = length(info_srf.name);
    switch Nsrf
        case 1
            CompatibleSurfOrder = true;
        case 2
            if ~strcmp(info_srf.name{1},'pial')...
                    || ~strcmp(info_srf.name{2},'white')
                CompatibleSurfOrder = false;
                disp(' ');
                warning('Info file related to pruning will not be generated.');
            else
                CompatibleSurfOrder = true;
            end
    end
end
h_mask = spm_vol(f_mask);
mask = spm_read_vols(h_mask);
inds_i = find(mask);

A0 = A;

if isempty(info_srf)
    f = batchgifti(f_srfs,h_mask.mat);
    Ap = hb_prune_adjacency(...
        A0,...
        mask,...
        f,...
        'EdgeTypeToRecover',opts.EdgeTypeToRecover,...
        'NodeClass',opts.NodeClass,...
        'parallelize', opts.parallelize);
    d1 = length(find(A0));
    d2 = length(find(Ap));
    fprintf('\n  %d edges removed. [all surfaces]',(d1-d2)/2);
else
    AA = cell(1,Nsrf);
    for n=1:Nsrf
        if n==1
            Ainit = A0;
        else
            Ainit = AA{n-1};
        end
        f = f_srfs(info_srf.index{n});
        f = batchgifti(f,h_mask.mat);
        AA{n} = hb_prune_adjacency(...
            Ainit,...
            mask,...
            f,...
            'EdgeTypeToRecover',opts.EdgeTypeToRecover,...
            'NodeClass',opts.NodeClass,...
            'parallelize', opts.parallelize);
        d1 = length(find(Ainit));
        d2 = length(find(AA{n}));
        fprintf('\n  %d edges removed. [%s surfaces]',...
            (d1-d2)/2,info_srf.name{n});
    end
    Ap = AA{Nsrf};
end

% Clean A.
[A,indA,mask] = cleanadjacency(Ap,1,mask);
d1 = length(inds_i);
d2 = length(find(indA));
fprintf('\n  %d nodes removed.',d1-d2);
indices = inds_i(indA);
if setdiff(indices,find(mask)), error('[HB] fishy.'); end
chk = setdiff(inds_i,indices); %removed nodes
if ~isequal(chk,inds_i(~indA)),error('fishy.'); end

% Volume highlighting removed voxels from f_mask
h = h_mask;
v = spm_read_vols(h);
v(inds_i(~indA)) = 2;
[p,n,e] = fileparts(h.fname); 
h.fname = fullfile(p,[n,'.postpruning_cleaned_nodes_marked',e]); 
spm_write_vol(h,v);
gzip(h.fname);
delete(h.fname);

pruning = struct;
if CompatibleSurfOrder
    switch Nsrf
        case 1
            d = A0-AA{1};
            switch info_srf.name{1}
                case 'pial'
                    pruning.A_diff_pre_post_pial_pruning = d;
                case 'white'
                    pruning.A_diff_pre_post_white_pruning = d;
                case 'unnamed'
                    pruning.A_diff_pre_post_pruning = d;
            end
        case 2
            for n=1:Nsrf
                switch info_srf.name{n}
                    case 'pial'
                        if n==1
                            pruning.A_diff_pre_post_pial_pruning = A0-AA{1};
                        else
                            pruning.A_diff_pre_post_pial_pruning = AA{1}-AA{2};
                        end
                    case 'white'
                        if n==1
                            pruning.A_diff_pre_post_white_pruning = A0-AA{1};
                        else
                            pruning.A_diff_pre_post_white_pruning = AA{1}-AA{2};
                        end
                    otherwise
                        error('Supported surface type names: ''pial'', ''white''.');
                        % If Nsrf==2, the only supported setting i that
                        % surfe type names be 'pial' and 'white'.
                end
            end
        otherwise
            error('?')
    end
end
pruning.A_removed_post_pruning = Ap(~indA,~indA);
pruning.ind_pre_pruning_A_remained_post_pruning = indA;
pruning.indices_of_mask_cleaned_post_pruning = inds_i(~indA);
pruning.mask_with_post_pruning_cleaned_voxels_marked = h.fname;

if CompatibleSurfOrder
    
    pruning.info = fullfile(p,[n,'.pruning_info.txt']);
    
    % Info file for plotting prunned edges.
    fid = fopen(pruning.info,'wt');
    fprintf(fid,'%%==========================================================================\n');
    fprintf(fid,'%%To plot prunned edges: \n');
    fprintf(fid,'%%A0: pre-pruning A. \n');
    switch Nsrf
        case 1
            switch info_srf.name{1}
                case 'pial'
                    fprintf(fid,'%%Ap: post-pial-pruning A. \n');
                    fprintf(fid,'\n');
                    fprintf(fid,'Adiff = G.pruning.A_diff_pre_post_pial_pruning; \n');
                case 'white'
                    fprintf(fid,'%%Ap: post-white-pruning A. \n');
                    fprintf(fid,'\n');
                    fprintf(fid,'Adiff = G.pruning.A_diff_pre_post_white_pruning; \n');
                case 'unnamed'
                    fprintf(fid,'%%Ap: post-pruning A. \n');
                    fprintf(fid,'\n');
                    fprintf(fid,'Adiff = G.pruning.A_diff_pre_post_pruning; \n');
            end
            fprintf(fid,'markCleanedVoxelsPostPruning = true; %%false \n');
            fprintf(fid,'Armd = G.pruning.A_removed_post_pruning; \n');
            fprintf(fid,'indA = G.pruning.ind_pre_pruning_A_remained_post_pruning;\n');
            fprintf(fid,'I = G.pruning.indices_of_mask_cleaned_post_pruning;\n');
            fprintf(fid,'Ap = Adiff-Adiff;     %%empty matrix \n');
            fprintf(fid,'Ap(indA,indA) = G.A;    %%place cleaned A \n');
            fprintf(fid,'Ap(~indA,~indA) = Armd; %%put back removed rows&columns \n');
            fprintf(fid,'A0 = Ap+Adiff;         %%put back pial-pruned edges \n');
            fprintf(fid,'hb_gunzip(G.f.mask); \n');
            fprintf(fid,'mask = spm_read_vols(spm_vol(G.f.mask)); \n');
            fprintf(fid,'if markCleanedVoxelsPostPruning \n');
            fprintf(fid,'    d = find(mask); \n');
            fprintf(fid,'    mask(d) = 2; \n');
            fprintf(fid,'    mask(I) = 1; \n');
            fprintf(fid,'end \n');
            fprintf(fid,'A_pre = {A0}; \n');
            fprintf(fid,'A_post = {Ap}; \n');
            fprintf(fid,'d3 = {''r''}; \n');
        case 2
            switch info_srf.name{1}
                case 'pial'
                    fprintf(fid,'%%App: post-pial-pruning A. \n');
                    fprintf(fid,'%%Apw: post-white-pruning A. \n');
                    fprintf(fid,'\n');
                    fprintf(fid,'markCleanedVoxelsPostPruning = true; %%false \n');
                    fprintf(fid,'Adiffp = G.pruning.A_diff_pre_post_pial_pruning; \n');
                    fprintf(fid,'Adiffw = G.pruning.A_diff_pre_post_white_pruning; \n');
                    fprintf(fid,'Armd = G.pruning.A_removed_post_pruning; \n');
                    fprintf(fid,'indA = G.pruning.ind_pre_pruning_A_remained_post_pruning;\n');
                    fprintf(fid,'I = G.pruning.indices_of_mask_cleaned_post_pruning;\n');
                    fprintf(fid,'Apw = Adiffw-Adiffw;     %%empty matrix \n');
                    fprintf(fid,'Apw(indA,indA) = G.A;    %%place cleaned A \n');
                    fprintf(fid,'Apw(~indA,~indA) = Armd; %%put back removed rows&columns \n');
                    fprintf(fid,'App = Apw+Adiffw;        %%put back white-pruned edges \n');
                    fprintf(fid,'A0 = App+Adiffp;         %%put back pial-pruned edges \n');
                    fprintf(fid,'hb_gunzip(G.f.mask); \n');
                    fprintf(fid,'mask = spm_read_vols(spm_vol(G.f.mask)); \n');
                    fprintf(fid,'if markCleanedVoxelsPostPruning \n');
                    fprintf(fid,'    d = find(mask); \n');
                    fprintf(fid,'    mask(d) = 2; \n');
                    fprintf(fid,'    mask(I) = 1; \n');
                    fprintf(fid,'end \n');
                    fprintf(fid,'A_pre = {A0;App}; \n');
                    fprintf(fid,'A_post = {App;Apw}; \n');
                    fprintf(fid,'d3 = {''r'';''g''}; \n');
                    
                case 'white'
                    error('Extend info file to support pruning order: white > pial.');
                    % Update text info to pruning scenario where first
                    % white surfaces are used and then the pila surfaces.
                    
                    
                    
            end
    end
    fprintf(fid,'figure; \n');
    fprintf(fid,'sl = 88;       %%axial slice number \n');
    fprintf(fid,'d1 = ''grid''; %%''diff'',''none'' \n');
    fprintf(fid,'d2 = true;     %%false \n');
    fprintf(fid,'hm_plot_adjacency_diff(A_pre,A_post,mask,sl,d1,d2,d3);\n');
    fprintf(fid,'%%==========================================================================\n');
    fclose(fid);
end

outs = struct;
outs.N = size(A,1);
outs.A = A;
outs.indices = indices;
outs.pruning = pruning;
end

function surf = batchgifti(file_names,mat)
% surf = BATCHGIFTI(file_names)
%       Loads multiple GIFTI files into a single mesh. The vertices will be
%       transformed into ACPC space.
% surf = BATCHGIFTI(file_names,mat)
%       mat will be used to transform the vertices from ACPC space such
%       that newV=mat\V.
%
% Martin Larsson, March 2017.

if nargin < 2
    mat = eye(4);
end

if ischar(file_names) % HB
    file_names = {file_names};
end

n = length(file_names);

surfs = cell(n,1);
face_count = 0;
vert_count = 0;

for i=1:n
    try
        surfs{i} = gifti(file_names{i});
    catch
        % Apparently, gifti(f) doesn't work on R2022a to read a gifti file,
        % matlab throwing an error saying you cannpt read file with gifti.m
        % (I guess since it's a class). As a workaround, the following
        % batch taken from within gifti.m
        d = read_gifti_file_hb(file_names{i},giftistruct);
        surfs{i} = class(d,'gifti');
    end
    face_count = face_count + size(surfs{i}.faces,1);
    vert_count = vert_count + size(surfs{i}.vertices,1);
    surfs{i}.vertices = transformvertices(mat\surfs{i}.mat,...
        surfs{i}.vertices);
end

surf.faces = zeros(face_count,3,class(surfs{i}.faces));
surf.vertices = zeros(vert_count,3,class(surfs{i}.vertices));

vert_pos = 1;
face_pos = 1;

for i=1:n
    faces = size(surfs{i}.faces,1);
    verts = size(surfs{i}.vertices,1);
    
    surf.faces(face_pos:face_pos+faces-1,:) = surfs{i}.faces +...
        vert_pos - 1;
    surf.vertices(vert_pos:vert_pos+verts-1,:) = surfs{i}.vertices;
    
    face_pos = face_pos + faces;
    vert_pos = vert_pos + verts;
end

    function Vt = transformvertices(A,V)
        % Transforms the 3D vertices in V using the 4x4 matrix A. V can be
        % either 3xN or Nx3. Vt is the transformed vertices and has the
        % same size as V.
        %   
        % Martin Larsson, March 2017.
        
        if size(V,1) == 3
            Vt = A * [V; ones(1,size(V,2))];
            Vt = Vt(1:3,:);
        elseif size(V,2) == 3
            Vt = (A * [V'; ones(1,size(V,1))])';
            Vt = Vt(:,1:3);
        else
            error('V must have either 3 rows or 3 columns.');
        end
    end
end

function [Ac,Aind,maskc] = cleanadjacency(A,comps,mask)
% Limit the number of connected components in the adjacency matrix A by
% only keeping the comps largest ones. Ac is the cleaned adjacency matrix
% and Aind is a logical vector indicating what nodes in A is kept in Ac,
% i.e. Ac=A(Aind,Aind). mask is a tensor whose non-zero elements, indexed
% linearly, correspond to the nodes in A. maskc is a copy of mask but with
% the removed nodes set to zero, thus Ac and maskc share the same
% relationship as A and mask do.
%
% Martin Larsson, March 2017.

bins = conncomp(graph(A));
bin_count = max(bins);

% no cleaning to be done
if bin_count <= comps
    Ac = A;
    Aind = true(1,size(A,1));
    maskc = mask;
    return;
end

% sort connected componens based on size
bins_size = histcounts(bins,1:bin_count+1);
[~,bin_ind] = sort(bins_size,'descend');

% keep the comps largest connected components
Aind = ismember(bins,bin_ind(1:comps));
Ac = A(Aind,Aind);

% clean mask
indices = find(mask);
maskc = mask;
maskc(indices(~Aind)) = 0;
end
%==========================================================================
function s = giftistruct
% Taken from within gifti.m
s = struct(...
    'metadata', ...
        struct(...
            'name',       {}, ...
            'value',      {} ...
        ), ...
    'label', ...
        struct(...
            'name',       {}, ...
            'index',      {} ...
        ), ...
    'data', ...
        struct(...
            'attributes', {}, ...
            'metadata',   struct('name',{}, 'value',{}), ...
            'space',      {}, ...
            'data',       {} ...
        ) ...
    );
end