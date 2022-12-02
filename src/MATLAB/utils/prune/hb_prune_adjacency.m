function Ap = hb_prune_adjacency(A,mask,surf,varargin)
%
%--------------------------------------------------------------------------
% Extended version of ml_prune_adjacency.m, Martin Larsson, March 2017.
% ML_PRUNE_ADJACENCY Remove edges from adjacency matrix using surface.
%   Ap = ML_PRUNE_ADJACENCY(A,mask,surf)
%       Prunes edges in the adjacency matrix A that intersect the provided
%       mesh surf. mask is a 3D tensor whose non-zero elements, indexed
%       linearly, correspond to the nodes in A.
%   Ap = ML_PRUNE_ADJACENCY(A,mask,surf,...)
%       Additional arguments can be provided as key-value pairs.
%         - conn6 (true/false) specifies whether or not A only has
%           6-connectivity. If 6-connectivity is used this can be set to
%           true for increased speed. For higher connectivity this must be
%           set to false (default).
%         - parallelize (true/false) specifies whether or not to run
%           prunings using multiple surfaces in parallel. Default is false.
%         - progess (true/false) specifies whether or not to report
%           progress during pruning. Default is false. If parallelization
%           is used no progress is reported.
% -------------------------------------------------------------------------
% Inputs:
%
%   A: Graph adjacency matrix to be pruned. 
%
%   mask: bw mask specifying the graph vertices in the image/volume.
%
%   surf: structure with field 'faces' and 'vertices'; gifti file names may
%   als be given.
%
%   Name-Value pair arguments:
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
%   If EdgeTypeToRecover=~'none', NodeClass should aslo be input. Note that
%   'x-y' is equivalent to 'y-x' since voxBG are undirected graphs.    
%
%   'NodeClass': (optional, but required if opts.EdgeTypeToRecover is
%   anything other than 'none'. A vector of size(A,1), with integers within
%   range [1,9], specifying the class of each graph vertex (voxel).
%
% Examples:
%   1. Ap = hb_prune_adjacency(A,m,s,'parallelize',1);   
%   2. NC = ones(size(A,1),1);
%      NC(1:floor(size(A,1)/3)) = 2; 
%      Ap = hb_prune_adjacency(A,m,s,'parallelize',1,...
%          'EdgeTypeToRecover',{'2-2','1-2'},...
%          'NodeClass',NC);
% 
% Hamid Behjat

p = inputParser;
addParameter(p,'conn6',false);
addParameter(p,'progress',false);
addParameter(p,'parallelize',false);
addParameter(p,'EdgeTypeToRecover','none',...
    @(x) iscell(x) || ... % cell array of chars
    and(ischar(x),or(strcmp(x,'none'),and(length(x)==3,x(2)=='-')))); % 'none' or 'x-y'
addParameter(p,'NodeClass',[],...
    @(x) or(isempty(x),...
    and(all(x<10),length(x)==size(A,1)))); % Current code support upto 9 classes

parse(p,varargin{:});
opts = p.Results;

assert(exist('TriangleRayIntersection.m','file')==2,...
    sprintf('Missing helping function: %s',...
    'TriangleRayIntersection.m'));

assert(exist('polygon2voxel.m','file')==2,...
    sprintf('Missing helping functions: %s, %s, %s',...
    'polygon2voxel.m',...
    'polygon2voxel_double.m',...
    'polygon2voxel_double.c'));

if ~isstruct(surf)
    assert(or(ischar(surf),iscell(surf)));
    surf = batchgifti(surf,h_mask.mat);
end

if ~opts.parallelize
    Ap = pruneadj(A,mask,surf,opts);
else
    % Divide faces evenly between workers.
    p = gcp;
    surfs = cell(p.NumWorkers,1);
    bounds = floor((0:p.NumWorkers)*size(surf.faces,1)/p.NumWorkers)+1;
    bounds(end) = size(surf.faces,1);
    
    for i=1:p.NumWorkers
        surfs{i}.faces = surf.faces(bounds(i):bounds(i+1),:);
        surfs{i}.vertices = surf.vertices;
    end
    
    % Parallel execution.
    opts.progres = false;
    spmd
        Aps = pruneadj(A,mask,surfs{labindex},opts);
    end
    
    % Combine results from workers.
    Ap = Aps{1};
    for i=2:p.NumWorkers
        Ap = Ap & Aps{i};
    end
    Ap = Ap.*A;
end

% Recover some pruned edges?
if ~iscell(opts.EdgeTypeToRecover)
    opts.EdgeTypeToRecover = {opts.EdgeTypeToRecover};
end
ETTR = opts.EdgeTypeToRecover;
NC = opts.NodeClass;

N = size(A,1);

Ip = setdiff(find(A),find(Ap)); % all pruned edges

[x,y] = ind2sub(N,Ip);

for iR=1:length(ETTR)
    etype = ETTR{iR};
    switch etype
        case 'none'
            % No edge to recover.
        otherwise % '1-1','2-2','1-2',...
            
            assert(length(etype)==3,...
                'number of classes more than 9; extend code.');
            
            c1 = str2double(etype(1));
            c2 = str2double(etype(3));
            if isequal(c1,c2)
                I = and(NC(x)==c1,NC(y)==c1);
            else
                I = or(...
                    and(NC(x)==c1,NC(y)==c2),...
                    and(NC(x)==c2,NC(y)==c1)...
                    );
            end
            assert(nnz(Ap(Ip(I)))==0,'fishy');
            
            % recover edges
            Aadd = sparse(x(I),y(I),I(I),N,N,nnz(I));
            Ap = or(Aadd,Ap).*A;
            %Ap(Ip(I)) = A(Ip(I)); % naive approach to update Ap; takes ages.
    end
end
end

function Ap = pruneadj(A,mask,g,opts)
% Prunes edges in the adjacency matrix A that intersect the surface g. mask
% is a volume indicating where the nodes of the graph are in 3D space. The
% non-zero elements in mask, indexed linearly, correspond to the row/column
% numbers of A. g is a triangulated mesh struct with the fields:
%   - vertices (Nx3)
%   - faces (Mx3)
%
% opts.conn6 specifies whether or not A only has 6-connectivity. If
% 6-connectivity is used this can be set to true for increased speed. For
% higher connectivity this must be set to false (default).
%
% opts.progress specifies whether or not to report progress.
%
% Original author: Martin Larsson, March 2017.

Ap = A;
dim = size(mask);
face_count = size(g.faces,1);

indices = find(mask);

assert(length(indices)==size(A,1));

indices_inv = zeros(dim);
indices_inv(indices) = 1:length(indices);
[indx,indy,indz] = ind2sub(dim,indices);

% The size of to_remove might be a bit excessively large and arbitrary.
% Not preallocating does not seem to significantly affect performance,
% however it is kept in case the behavior is different for other
% versions of MATLAB.
to_remove = zeros(nnz(Ap),1);
to_remove_index = 1;

% It is faster to precalculate the triangles here than to extract the
% relevant vertices in the for loop.
all_tri = zeros(3,3,face_count);
all_tri(1,:,:) = g.vertices(g.faces(:,1),:)';
all_tri(2,:,:) = g.vertices(g.faces(:,2),:)';
all_tri(3,:,:) = g.vertices(g.faces(:,3),:)';

tri.faces = [1 2 3];

if opts.progress
    fprintf('* Pruning adjacency matrix...\nProgress: ');
    fprintf('%6d/%-6d\n',1,face_count);
end

% Handle one triangle at a time.
for i = 1:face_count
    if opts.progress && (~rem(i,500) || i == face_count)
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b%6d/%-6d\n',i,...
            face_count);
    end
    
    tri.vertices = all_tri(:,:,i);
    
    % Construct bounding box for the triangle.
    minv = floor(min(tri.vertices));
    maxv = ceil(max(tri.vertices));
    
    % Expand to accommodate dilation below.
    if ~opts.conn6
        minv = minv - 1;
        maxv = maxv + 1;
    end
    
    minv = max(minv, 1);
    maxv = min(maxv, dim);
    
    bbdim = maxv - minv + 1;
    
    % Translate vertices from the mask space to the bounding box.
    tri.vertices = tri.vertices - repmat(minv,3,1) + 1;
    
    % Voxelization of triangle.
    V = polygon2voxel(tri, bbdim, 'none', false);
    
    % Perform dilation with 6-connectivity. This is not pretty but it
    % is apparently faster than using the imdilate function. This is
    % not necessary if A only has 6-connectivity, but is required to
    % catch "diagonal" edges.
    if ~opts.conn6
        V2 = V;
        V2(1:end-1,:,:) = V2(1:end-1,:,:) | V(2:end,:,:);
        V2(:,1:end-1,:) = V2(:,1:end-1,:) | V(:,2:end,:);
        V2(:,:,1:end-1) = V2(:,:,1:end-1) | V(:,:,2:end);
        V2(2:end,:,:) = V2(2:end,:,:) | V(1:end-1,:,:);
        V2(:,2:end,:) = V2(:,2:end,:) | V(:,1:end-1,:);
        V2(:,:,2:end) = V2(:,:,2:end) | V(:,:,1:end-1);
        V = V2;
    end
    
    % Remove voxels not in the mask.
    V = V & mask(minv(1):maxv(1),minv(2):maxv(2),minv(3):maxv(3));
    
    % Translate voxels back to the mask space.
    [vx,vy,vz] = ind2sub(bbdim,find(V));
    if size(vx,1) == 0
        continue;
    end
    voxels = [vx vy vz];
    voxels = voxels + repmat(minv,size(voxels,1),1) - 1;
    
    % Convert mask space (x,y,z) to index in A.
    Aind = indices_inv(sub2ind(dim,voxels(:,1),voxels(:,2),...
        voxels(:,3)));
    
    % Define start and end point for all edges from the voxels.
    edge_index = 0;
    for j = 1:length(Aind)
        neighbors = find(Ap(:,Aind(j)));
        nn = length(neighbors);
        startv = repmat(voxels(j,:),nn,1);
        endv = [indx(neighbors) indy(neighbors) indz(neighbors)];
        
        orig(edge_index+1:edge_index+nn,:) = startv;
        dir(edge_index+1:edge_index+nn,:) = endv - startv;
        
        edges(edge_index+1:edge_index+nn) = sub2ind(size(Ap),...
            neighbors,repmat(Aind(j),nn,1));
        edge_index = edge_index + nn;
    end
    
    % There were no edges to check intersection with.
    if edge_index == 0
        continue;
    end
    
    % Check for intersections.
    intersections = TriangleRayIntersection(orig(1:edge_index,:),...
        dir(1:edge_index,:), all_tri(1,:,i), all_tri(2,:,i),...
        all_tri(3,:,i), 'lineType', 'segment', 'border', 'inclusive');
    
    edges = edges(intersections);
    
    % Append edges to list to be removed in bulk later.
    to_remove(to_remove_index+(0:length(edges)-1)) = edges;
    to_remove_index = to_remove_index + length(edges);
end

% Remove edges.
Ap(to_remove(1:to_remove_index-1)) = 0;

% Enforce symmetry.
Ap = (Ap & Ap').*Ap;
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
    surfs{i} = gifti(file_names{i});
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

