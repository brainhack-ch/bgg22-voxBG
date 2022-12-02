function [G,sts] = voxbg_adjacencymatrix_diffusion(G,opts,chk)
%HB_ADJACENCYMATRIX_DIFFUSION builds adjacency matrix using diffuion MRI data.
%
%
% Hamid Behjat

% Get data.
sts = pvt_get_diffusion(G);
if sts==false
    return;
end

% Build ODFs/FODs.
if opts.temp_skip_subjs_with_diff_problem
    try
        switch G.diffusionModel
            case 'odf'
                [G,sts] = pvt_get_odfs(G,opts,chk);
            case 'fod'
                [G,sts] = pvt_get_fods(G,opts,chk);
        end
    catch
        sts = false;
        fprintf('\n..Subject skipped: %s [problem with getting ODFs]',...
            G.subject);
    end
else
    switch G.diffusionModel
        case 'odf'
            [G,sts] = pvt_get_odfs(G,opts,chk);
        case 'fod'
            [G,sts] = pvt_get_fods(G,opts,chk);
    end
end

if sts==false
    return;
end

if opts.do_just_get_odfs
    return;
end

%-Update graph masks + ODF file.
% (if there exists some voxels in mask that had no ODF)

load(G.f.odf.fib_mat,'MaskIndsNoODF');

if not(isempty(MaskIndsNoODF))
    
    fprintf('\n..Updating masks.. [removing voxels with no ODF]');
    
    f_masks = {
        'mask'
        'mask_wm'
        'mask_gm'
        'mask_gm_ctx'
        'mask_gm_subctx'
        'mask_gm_ctx_subctx_common'
        'mask_wm_gmborder'
        'mask_gm_wmborder'
        };

    for iF=1:length(f_masks)
        d = G.f.(f_masks{iF});
        [v,h] = hb_nii_load(d);
        if strcmp(f_masks{iF},'mask')
            col_NoODF = find(ismember(find(v),MaskIndsNoODF)); % columns in ODF file to remove
        end
        v(MaskIndsNoODF) = 0;
        spm_write_vol(h,v);
    end
    
    fprintf('\n..Updating ODF file.. [removing no-ODF-voxel columns]');
    load(G.f.odf.fib_mat,'odf');
    d = odf(:,col_NoODF);
    assert(nnz(d)==0,'Fishy. These columns should be zero.');
    odf(:,col_NoODF) = []; % remove columns
    save(G.f.odf.fib_mat,'odf','-append');
    MaskIndsNoODF = []; 
    save(G.f.odf.fib_mat,'MaskIndsNoODF','-append');
end

switch G.tissue
    case 'cerebrum'
        opts = checkOptsRePostPruning(opts);
        
        % Node classes, required for:
        % - pruning
        % - removing outer-layer-neighb5 GM-GM edges
        if opts.UpdateMasksPostPruning
            [opts.NC,hmasks,vmasks] = getNC(G);
        else
            opts.NC = getNC(G);
        end
end

%-Build A----------------------------------------------------------
fprintf('\n..Determining voxel adjacencies.. [Moore neighbourhood]')
[G,sts] = pvt_get_A(G,opts);

switch G.tissue
    case 'cerebrum'
        
        %-Pruning A--------------------------------------------------------
        fprintf('\n..Pruning edges that pass through pial..');
        
        % surface info
        f_srfs = G.f.surface.pial;
        info_srf = struct;
        info_srf.name = {'pial'};
        info_srf.index = {[1 2]};
        
        % prune
        outs = hb_prune_graph(...
            G.A,...
            G.f.mask,...
            f_srfs,...
            'info_srf',info_srf,...
            'EdgeTypeToRecover',{'3-3','3-1','3-2','2-2','1-2'},... % *note
            'NodeClass',opts.NC,...
            'parallelize',~opts.runPar);
        
        % *note
        % '3-3': WM-WM edges
        % '3-1': WM-GM(ctx) edges
        % '3-2': WM-GM(sub-ctx) edges
        % '2-2': GM(sub-ctx)-GM(sub-ctx) edges
        % '1-2': GM(ctx)-GM(sub-ctx) edges
        % The above edges, if pruned, will be recovered. Thus, only edges
        % of type '1-1', i.e., GM(ctx)-GM(ctx), are allowed to be pruned.
        
        N_preprune    = G.N;
        inds_preprune = G.indices;
        
        % Update G.
        G.N = outs.N;
        G.A = outs.A;
        G.indices = outs.indices;
        G.pruning = outs.pruning;
        
        % remove voxels associated to removed nodes
        d = setdiff(...
            1:N_preprune,...
            find(G.pruning.ind_pre_pruning_A_remained_post_pruning));
        inds_rm = inds_preprune(d);
        
        hmasks.mask        = spm_vol(G.f.mask);
        hmasks.gm          = spm_vol(G.f.mask_gm);
        hmasks.gm_wmborder = spm_vol(G.f.mask_gm_wmborder);
        hmasks.wm_gmborder = spm_vol(G.f.mask_wm_gmborder);
        
        vmasks.mask        = spm_read_vols(hmasks.mask);
        vmasks.gm          = spm_read_vols(hmasks.gm);
        vmasks.gm_wmborder = spm_read_vols(hmasks.gm_wmborder);
        vmasks.wm_gmborder = spm_read_vols(hmasks.wm_gmborder);
        
        vmasks.mask(inds_rm)        = 0;
        vmasks.gm(inds_rm)          = 0;
        vmasks.gm_wmborder(inds_rm) = 0;
        vmasks.wm_gmborder(inds_rm) = 0;
        vmasks.wm(inds_rm)          = 0;
        vmasks.gm_ctx(inds_rm)      = 0;
        vmasks.gm_subctx(inds_rm)   = 0;
        vmasks.gm_cmn(inds_rm)      = 0;
        
        % within vol indices of all submasks
        G.indices_gm          = find(vmasks.gm);
        G.indices_gm_wmborder = find(vmasks.gm_wmborder);
        G.indices_wm_gmborder = find(vmasks.wm_gmborder);
        G.indices_wm          = find(vmasks.wm);
        G.indices_gm_ctx      = find(vmasks.gm_ctx);
        G.indices_gm_subctx   = find(vmasks.gm_subctx);
        G.indices_gm_common   = find(vmasks.gm_cmn);
        
        % Update masks? [no reason why not to, right?]
        if opts.UpdateMasksPostPruning
            fprintf('\n..Updating masks.. [removing voxels for removed nodes]');
            spm_write_vol(hmasks.mask,vmasks.mask);
            spm_write_vol(hmasks.gm,vmasks.gm);
            spm_write_vol(hmasks.gm_wmborder,vmasks.gm_wmborder);
            spm_write_vol(hmasks.wm_gmborder,vmasks.wm_gmborder);
            spm_write_vol(hmasks.gm_cmn,vmasks.gm_cmn);
            spm_write_vol(hmasks.gm_subctx,vmasks.gm_subctx);
            spm_write_vol(hmasks.gm_ctx,vmasks.gm_ctx);
            spm_write_vol(hmasks.wm,vmasks.wm);
        end
        
        % Update odf file [remove odf etc. associated to removed nodes]
        %---
        % Note: not necessarily a good thing to do, since when done, the
        % file can not be used for other versions of the cerebrum graph,
        % and thus, DSI_studio has to be run again.
        %---
        if opts.UpdateODFmatfilePostPruning 
         
            fprintf(['\n..Updating odf mat file.. ',...
                '[removing odfs for removed nodes]']);
            d = load(G.f.odf.fib_mat,'odf','odf_vertices');
            odf = d.odf;
            d = G.pruning.ind_pre_pruning_A_remained_post_pruning;
            odf = odf(:,d);
            save(G.f.odf.fib_mat,'odf','-append')
        end
        
        %-GM edge weights--------------------------------------------------
        % GM-GM edges
        % GM-WM edges
        switch G.gmEdgeWeightVers
            case 1
                fprintf(...
                    '\n..GM-GM/WM edges weighted just as WM-WM edges.');
            case {2,3}
                I = find(G.A);
                AI = full(G.A(I));
                [x,y] = ind2sub(G.N,I);
                NC = opts.NC;
                
                % GM-GM edges [~35% of edges]
                Igmgm = ...
                    and(NC(x)==1,NC(y)==1) |... GM(ctx)-GM(ctx)
                    and(NC(x)==2,NC(y)==2) |... GM(sub-ctx)-GM(sub-ctx)
                    and(NC(x)==1,NC(y)==2) |... GM(ctx)-GM(sub-ctx)
                    and(NC(x)==2,NC(y)==1);   % GM(sub-ctx)-GM(ctx)
                
                % GM-WM edges [~40% of edges]
                Igmwm = ...
                    and(NC(x)==1,NC(y)==3) |... GM(ctx)-WM
                    and(NC(x)==2,NC(y)==3) |... GM(sub-ctx)-WM
                    and(NC(x)==3,NC(y)==1) |... WM-GM(ctx)
                    and(NC(x)==3,NC(y)==2);   % WM-GM(sub-ctx) %*** bug fixed!! [11.04.2022]
                
                switch G.gmEdgeWeightVers
                    
                    case 2
                        fprintf(...
                            ['\n..Removing GM-GM edge weights.. ',...
                            '[weights set to 1]']);
                        
                        % makes no sense to assign all Igmwm weights to 1.
                        % GM-GM/WM edges [~75% of edges]
                        %Igm = or(Igmgm,Igmwm);
                        
                        % remove weights
                        AI(Igmgm) = 1;
                        
                    case 3
                        
                        fprintf(...
                            ['\n..Weighting GM-GM edges.. ',...
                            '[inverse Euclidean distance]']);
                        
                        % get GM node coordinates
                        gm1 = G.indices(x(Igmgm));
                        gm2 = G.indices(y(Igmgm));
                        
                        c1 = zeros(size(gm1,1),3);
                        c2 = zeros(size(gm2,1),3);
                        [c1(:,1),c1(:,2),c1(:,3)] = ind2sub(G.dim,gm1);
                        [c2(:,1),c2(:,2),c2(:,3)] = ind2sub(G.dim,gm2);
                        
                        % set inverse Euclidean distances as weight
                        AI(Igmgm) = 1./vecnorm(c1-c2,2,2);
                end
                
                % update A
                G.A = sparse(x,y,AI,G.N,G.N,nnz(AI));
        end
end
end

%==========================================================================
function opts = checkOptsRePostPruning(opts)

if not(isfield(opts,'UpdateMasksPostPruning'))
    opts.UpdateMasksPostPruning = true;
end

if not(isfield(opts,'UpdateODFmatfilePostPruning'))
    
    opts.UpdateODFmatfilePostPruning = false;
    
    % This setting is for the cerebrum graph design.
    %
    % If set to true, the ODF mat file will be updated after the pruning
    % phase to remove ODFs aociated to graph nodes(voxels) that have been
    % removed after pruning. That is, the number of the ODFs in the updated
    % file will be equal to G.N.
    %
    % Note that there is limitation in doing this update. If the update is
    % done, and then another version of the cerebrum graph is attempted to
    % be designed, the code will see that the required ODF mat file exists,
    % so it will not call DSIstudio to regenrate them. But since a number
    % of voxels have been removed, the algorithm will run into problem when
    % attempting to build the adjacency matrix. As such it is more robust
    % to not do this update, but if then the ODFs need to be inspected for
    % some post processing or whatever, it can then be very simply updated
    % just done as here in the code after loading the file as:
    %
    % odf = load(G.f.odf.fib_mat,'odf').odf;
    % odf = odf(:,G.pruning.ind_pre_pruning_A_remained_post_pruning);
end
end

%==========================================================================
function sts = pvt_get_diffusion(G)
% The diffusion data is zipped. Either DSI_studio needs to temporarily
% unzipped them during the process of computing .src file or in craetes
% some temporary files in the same directory where the data are. In either
% case, data should be transferred to a writable drive.

if strcmp(G.f.diffusion,G.f.diffusion_save)
    sts = true;
    return;
end

if exist(fullfile(G.f.diffusion_save,'data.nii.gz'),'file')
    sts = true;
    return;
end

fprintf('..Transfering diffusion data to writable local directory..\n')
fprintf('  This will take some time.. 3-10 minutes.\n');
fprintf('  strated:');
disp(datetime);

d1 = G.f.diffusion;
d2 = G.f.diffusion_save;
if ~exist(G.f.diffusion,'dir')
    warning('No diffusion data; skipping WM graph construction.');
    sts = false;
    return;
end

sts1 = copyfile(fullfile(d1,'bvals'),fullfile(d2,'bvals'));
sts2 = copyfile(fullfile(d1,'bvecs'),fullfile(d2,'bvecs'));
sts3 = copyfile(fullfile(d1,'data.nii.gz'),fullfile(d2,'data.nii.gz'));

if ~all([sts1,sts2,sts3])
    error('[HB] Error in transfering files.')
end

fprintf('  Ended:');
disp(datetime);
fprintf('..done.\n');
sts = true;
end

%==========================================================================
function [G,sts] = pvt_get_odfs(G,opts,owrt)

opts2 = struct;
opts2.N_fibers = G.N_fibers;
opts2.dir_save = G.f.diffusion_save;
opts2.prefix = [G.tissue,G.maskTypeTag,'.'];
if ~owrt
    if opts.DEBUG_on_aitchbi
        opts2.overwrite = false;
    else
        opts2.overwrite = true;
    end
else
    opts2.overwrite = true;
end

%try
    outs = hb_dsistudio_get_odfs(...
        G.f.diffusion_data,...
        G.f.mask,...
        opts.dsi_root,...
        opts2);
    sts = true;
%catch
%    sts = false;
%    return;
%end

% update G
G.f.odf.fib_mat = outs.f_fibmat;
G.f.odf.delete.srcgz = outs.f_srcgz;
if 1
    % since the naming is not adjusted for the mask used, becomes
    % problematic when doing defferent designs, with different masks. The
    % name of outs.f_fibmat is perfectly adjusted, so no problem to keep.
    % File outs.f_srcgz is also not affected by a mask, so ok to keep
    % an reuse.
    delete(outs.f_fibgz);
    delete(outs.f_fib);
else
    G.f.odf.delete.fibgz = outs.f_fibgz;
    G.f.odf.delete.fib = outs.f_fib;
end
end

%==========================================================================
function [G,sts] = pvt_get_fods(G,opts,chk)

if chk
    tag_force = ' -force ';
else
    tag_force = '';
end

switch G.maskType
    case 'onlyWM'
        fod_method = 'csd';
    case {'subcortical','wholebrain'}
        fod_method = 'msmtCsd';
end

dir_pwd = pwd;

cd(G.f.diffusion_save);

%-Estimate response function.
fprintf('\n\n..Estimating response function..');
switch fod_method
    case 'csd'
        command = [opts.mrtrix3_root, filesep, 'dwi2response tournier -fslgrad bvecs bvals -mask ',G.f.mask, tag_force, ' data.nii.gz wm_response.txt'];
        [sts,logfile] = system(command);
    case 'msmtCsd'
        command = [opts.mrtrix3_root, filesep, 'dwi2response dhollander -fslgrad bvecs bvals -mask ',G.f.mask, tag_force, ' data.nii.gz wm_response.txt gm_response.txt csf_response.txt'];
        [sts,logfile] = system(command);
end
chk_sts(sts,logfile);

%-Build fods.
fprintf('\n\n..Buildfing FODs..');
switch fod_method
    case 'csd'
        command = [opts.mrtrix3_root, filesep, 'dwi2fod csd -fslgrad bvecs bvals -mask ',G.f.mask, tag_force, ' data.nii.gz wm_response.txt fods.mif'];
        [sts,logfile] = system(command);
    case 'msmtCsd'
        command = [opts.mrtrix3_root, filesep, 'dwi2fod msmt_csd -fslgrad bvecs bvals -mask ',G.f.mask, tag_force, ' data.nii.gz wm_response.txt fods.mif gm_response.txt gm.mif csf_response.txt csf.mif'];
        [sts,logfile] = system(command);
end
chk_sts(sts,logfile);

%-Extract fixels.
d = ['fixels_',fod_method,'_max',num2str(G.N_fixels)];
switch G.maxAfd
    case 0.1
        dir_fixel = [d,'_peakvalue0point1'];
    otherwise
        error('extend..')
end
G.f.dir_fixel = fullfile(G.f.diffusion_save,dir_fixel);

fprintf('\n\n..Extracting fixels..');
command = [opts.mrtrix3_root, filesep, 'fod2fixel -fmls_peak_value ',num2str(G.maxAfd),' -maxnum ',num2str(G.N_fixels),' -mask ',G.f.mask, tag_force, ' -afd afd.mif -peak_amp peak_amplitude.mif -disp dispersion.mif fods.mif ',dir_fixel];
[sts,logfile] = system(command);
rerun = chk_sts(sts,logfile,dir_fixel);
if rerun
    [sts,logfile] = system(command);
    rerun = chk_sts(sts,logfile,dir_fixel);
    if rerun, error('fishy'); end
end

% 1D to 3D.
fprintf('\n\n..Changing 1D .mif to volume..');
[sts,logfile] = system([opts.mrtrix3_root, filesep, 'fixel2voxel ', tag_force, fullfile(dir_fixel,'afd.mif'),' none ', fullfile(dir_fixel,'afd_3d.mif')]);
%[sts,logfile] = system([opts.mrtrix3_root, filesep, 'fixel2voxel ', tag_force, 'peak_amplitude.mif none peak_amplitude_3d.mif']);
chk_sts(sts,logfile);

sts = 1;

cd(dir_pwd);

end

%==========================================================================
function rerun = chk_sts(sts,logfile,dir_fixel)
if ~exist('dir_fixel','var')
    dir_fixel = [];
end
rerun = 0;
switch sts
    case 0
        fprintf(' done.')
    case 1
        if contains(logfile,'already exists')
            fprintf(' file already exists.')
        elseif contains(logfile,'-force option cannot safely be applied on directories')
            fprintf(' erasing fixel directory contents to overwrite.. ')
            delete([dir_fixel,filesep,'*']);
            rerun = 1;
        else
            error('fishy; check logfile');
        end
    otherwise
        error('fishy; check logfile');
end
end

%==========================================================================
function [G,sts] = pvt_get_A(G,opts)

% Determine edges: find neighboring voxels & direction of edges.
[dirs,ndirs,~,npar] = ml_create_neighborhood(G.neighb,G.rempar);
Nn = size(dirs,2);
%note: Nn & npar are not the same if parallel edges not removed.

switch G.diffusionModel
    
    case 'odf'
        
        %-Tissue mask etc.
        hb_gunzip(G.f.mask);
        h_mask = spm_vol(G.f.mask);
        mask = spm_read_vols(h_mask);
        indices = find(mask);
        dim = size(mask);
        Ng = length(indices);
        
        %-ci ni lni.
        [ci,ni,lni] = get_cinilni_voxel(indices,dirs,dim,Nn);
        
        % Load ODFs.
        d = load(G.f.odf.fib_mat,'odf','odf_vertices');
        odfs = d.odf;
        vertices = d.odf_vertices;
        %[fa,odfs,index,vertices,faces] = ml_dsi_load_odfs(G.f.odf.fib_mat,true); % See Note 1.
        
        %if G.minnorm
        %   odfs = odfs-min(odfs);
        %end
        
        odfs = odfs.^(G.odfPow);
        odfs = [odfs; odfs];
        
        if G.magnitude
            d = load(G.f.odf.fib_mat,'fa0');
            fas = d.fa0;
            %fas = fa(1,fa(1,:)~=0);
        end
        
        % Determine weights.
        switch G.method
            case 'voronoi'
                cosangs = ndirs'*vertices;
                cosangs(abs(cosangs)>1) = 1;
                angs = acos(cosangs);
                imin = angs == repmat(min(angs),Nn,1);
                voronoi_set = imin./repmat(sum(imin),Nn,1);
                Pdiff = voronoi_set*odfs.*lni;
            case 'soliang'
                % see NOTE 2.
                omega    = 4*pi/npar;
                costheta = 1-omega/(2*pi);
                cone_set = ndirs'*vertices > costheta;
                cone_set = cone_set./repmat(sum(cone_set,2),1,size(cone_set,2)); %see NOTE 3.
                Pdiff = cone_set*odfs.*lni;
            case {'soliangshAverage','soliangshAverageuniform','soliangshCenter'}
                shDeg = G.shOrder; %8;
                angles = cart2sph_phys(vertices');
                [~,pT] = makePT(shDeg,angles);
                sh = pT*odfs;
                omega = 4*pi/npar;
                costheta = 1-omega/(2*pi);
                switch G.method
                    case 'soliangshAverage'
                        d = icosphere(5);
                        vertices_icosph = d.Vertices';
                        angles_icosphere = cart2sph_phys(vertices_icosph');
                        T = makePT(shDeg,angles_icosphere);
                        odfs_icosphere = T*sh;
                        cone_set = ndirs'*vertices_icosph > costheta;
                        cone_set = cone_set./repmat(sum(cone_set,2),1,size(cone_set,2));
                        Pdiff = cone_set*odfs_icosphere.*lni;
                    case 'soliangshAverageuniform'
                        [odfssum,Nt] = hb_soliangodf(odfs,vertices,ndirs,omega);
                        Pdiff = (odfssum/Nt).*lni;
                    case 'soliangshCenter'
                        angles_neighb = cart2sph_phys(ndirs');
                        T = makePT(shDeg,angles_neighb);
                        odfs_neighb = T*sh;
                        Pdiff = odfs_neighb.*lni;
                    case 'soliangIntegrate'
                        %zdirs = ndirs(:,or(all(ndirs==zdir),all(ndirs==-zdir)));
                        %d1 = any(zdirs'*vertices_icosph>costheta);
                end
        end
        
        % Prune edges. [hybrid neighb design]
        if all([...
                isequal(G.tissue,'cerebrum'),...
                G.neighb==5,...
                G.neighbHybrid...
                ])
            
            N_lni = nnz(lni);
            
            % get node coordinates in volume
            cs = zeros(3,length(indices));
            [cs(1,:),cs(2,:),cs(3,:)]= ind2sub(G.dim,indices);
            
            % determine GM nodes
            n_gm = or(opts.NC==1,opts.NC==2);
            
            % Removel phase 1/2--------------------------------------------  
            fprintf('\n...Removing long-distance GM-GM edges..')
            % only keep GM-GM edges of type neighb3, i.e., those of
            % length less than sqrt(3), i.e., edges connecting a
            % voxel to its neighb5-outer-layer voxels.
            
            % get Euclidean distance of edges
            xx = ci(:);
            yy = ni(:);
            yy(yy==0) = find(n_gm==0,1); % dirty hack; assign zeros to a non-GM node, to prevent error when calling yy as index, as in e.g. cs(*,yy) or node_gm(yy)
            d = zeros(3,length(xx));
            d(1,:) = cs(1,xx)-cs(1,yy);
            d(2,:) = cs(2,xx)-cs(2,yy);
            d(3,:) = cs(3,xx)-cs(3,yy);
            ed = vecnorm(d);
            ed = ed(:);
            
            % determine GM-GM edges
            e_gmgm = and(n_gm(xx),n_gm(yy));
            
            % anything fishy? [related to the dirty hack above]
            assert(nnz(and(ni(:)==0,e_gmgm))==0);
            d=ed;
            d(ni(:)==0)=0;
            if G.rempar
                isequal(...
                    sort(unique(d)),...
                    sqrt([0 1 2 3 5 6 9]')... % possible edge distances for 98-conn neighb5
                    );
            else
                error('extend to 124-conn neighb5 edge distances.');
            end
            
            % GM-GM edges to remove
            e_rm = and(e_gmgm,ed>sqrt(3));
            e_rm = reshape(e_rm,size(lni));
            
            % report
            fprintf(...
                ['\n....# of edges removed: %d ',...
                '[%.1f%% of total graph edges]'],...
                nnz(e_rm),100*(nnz(e_rm)/N_lni));
            
            % update lni
            lni = and(lni,not(e_rm));
            %--------------------------------------------------------------
            
            if 0
            % Removel phase 2/2--------------------------------------------
            fprintf(...
                ['\n...Removing long-distance GM-WM edges.. ',...
                '[if GM voxel is deep, not touching WM in its 26-conn]']);
            % only keep GM-WM long edges (i.e., with length larger than
            % sqrt(3)) if the GM voxel at the GM/WM boundary, since
            % otherwise it is not intuitive to connect a WM voxel to a GM
            % voxel that does not have a boundary with WM; in this
            % scenario, we remove this edge, but still, if the weight is
            % high, meaning that it really is pointing toards that
            % direction, then the weight of another edge from that WM voxel
            % pointing towards a surroung GM voxel (at the GM/WM) is also
            % notably high, and that will be kept. This procedure results
            % in a more intuitive and realistic diffusion behavior at the
            % GM/WM boundary, when implementing GSR or GSD.
            
            n_wm = opts.NC==3;
            
            xx = ci(:);
            yy = ni(:);
            
            % remove non-edge terms
            ind_ni = find(yy);
            xx = xx(ind_ni);
            yy = yy(ind_ni);
            
            % Keep edges in e_gmwm_p1 for which their GM node is at the
            % boundary to WM; the remaining edges (which will be
            % removed) connect a WM voxel to a GM voxel that is deep
            % (i.e., does not touch WM), which is not realistic to
            % keep.
            v = spm_read_vols(spm_vol(G.f.mask_gm_wmborder)); % WM border voxels in GM
            ind_gmborder = find(v);
            e_gmwm = and(n_gm(xx),n_wm(yy));            % GM-WM edges to inspect
            e_wmgm = and(n_wm(xx),n_gm(yy));            % WM-GM edges to inspect
            
            ind_gm1 = indices(xx(e_gmwm));
            ind_gm2 = indices(yy(e_wmgm));

            %ind_gm1 = indices(xx(n_gm(xx(e_gmwm))));
            %ind_gm2 = indices(yy(n_gm(yy(e_wmgm))));
            
            e_gmwm(ismember(ind_gm1,ind_gmborder)) = 0; % GM-WM edges connecting to a GM voxel that touches WM are kept
            e_wmgm(ismember(ind_gm2,ind_gmborder)) = 0; % WM-GM edges connecting to a GM voxel that touches WM are kept
            e_rm = or(e_gmwm,e_wmgm);
            
            % reshape
            d = zeros(size(ni(:)));
            d(ind_ni) = e_rm;
            e_rm = d;
            e_rm = reshape(e_rm,size(lni));
            
            % report
            fprintf(...
                ['\n....# of edges removed: %d ',...
                '[%.1f%% of total graph edges]'],...
                nnz(e_rm),100*(nnz(e_rm)/N_lni));
            
            % update lni
            lni = and(lni,not(e_rm));
            %--------------------------------------------------------------
            end
            
            % Update Pdiff
            Pdiff = Pdiff.*lni;
        end
        
        Pdiff = 0.5*Pdiff./repmat(max(Pdiff),Nn,1);
        
        %if G.thresh>0
        %    th = G.thresh/100/2;
        %    Pdiff(Pdiff<th) = 0;
        %end
        
        % Build A.
        A = sparse(ci(lni),ni(lni),Pdiff(lni),Ng,Ng);
        A = A+A';
        
        if G.thresh>0
            I = full(find(A));
            edge_weights = full(A(I));
            mu = @(x,a,b) (x*(1-a)).^b ./ ((x*(1-a)).^b + (a*(1-x)).^b);  % Tunable sigmoid function
            a = G.thresh/100; %Ath
            b = G.beta; %50;
            weights = mu(edge_weights, a, b);
            
            edge_weights2 = edge_weights .* weights;
            
            if 0 %debug check
                debug_chk1(mu,a,b,edge_weights,edge_weights2); %#ok<UNRCH>
            end
            
            % Re-build adjacency matrix
            [i,j] = ind2sub([Ng, Ng], I);
            A = sparse(i,j,edge_weights2);
            
            clear edge_weights edge_weights2
        end
        
        if G.magnitude
            fas = fas(:)/max(fas);
            fas98 = repelem(fas',Nn,1); % 98xNg
            fas98 = fas98.*lni;
            M = sparse(ci(lni),ni(lni),fas98(lni),Ng,Ng);
            M = M.*M';
            AM = M.*A;
            
            if 1 % debug check
                debug_chk2(G.f.wm_graphspace,mask,indices,A,M,AM,fas);
            end
            
            A = AM;
            clear M AM
        end
        
        % Update G.
        G.N = Ng;
        G.A = A;
        G.indices = indices;
        G.private.directions_neighb = ndirs;
        G.private.directions_odf = vertices;
        G.private.costheta = costheta;
        %G.private.cone_set = cone_set;
        %G.dim = wHeader.dim;
        %[G.xyz(:,1),G.xyz(:,2),G.xyz(:,3)] = ind2sub(G.dim,indices);
        %G.odf.odf = odf;
        %G.odf.odf_vertices = vertices';
        %G.odf.odf_faces = faces';
        %G.odf.index = index(:,indices);
        
        % Clean up.
        if isfield(opts,'DEBUG_on_aitchbi')
            if ~opts.DEBUG_on_aitchbi
                d = fieldnames(G.f.odf.delete);
                for iD=1:length(d)
                    delete(G.f.odf.delete.(d{iD}));
                end
                %gzip(G.f.odf.fib_mat); % almost no save in space by compressing
                %delete(G.f.odf.fib_mat);
            end
        end
        sts = true;
        
    case 'fod'
        
        %-Tissue mask.
        hb_gunzip(G.f.mask);
        h_mask = spm_vol(G.f.mask);
        
        % Load FODs.
        [~,afds,index,directions] = load_fod_etc(G);
        Nf = directions.dim(1); % total number of fixels
        fixel_dirs = directions.data;
        
        %-adjust volume orientation to match loaded mask.
        %[so that G.indices becomes consistent with 'indices_fod']
        d = index.data;
        if contains(index.layout(1:9),'-')
            if strcmp(index.layout(1),'-') %x
                if h_mask.mat(1,1)<0
                    d = flip(d,1);
                end
            end
            if strcmp(index.layout(4),'-') %y
                if h_mask.mat(2,2)<0
                    d = flip(d,2);
                end
            end
            if strcmp(index.layout(7),'-') %z
                if h_mask.mat(3,3)<0
                    d = flip(d,3);
                end
            end
        end
        indices_fod = find(d(:,:,:,1)); % voxels with at least 1 fixel
        
        num_fixels = d(:,:,:,1);
        num_fixels = num_fixels(indices_fod);
        
        inds_start = d(:,:,:,2)+1; % +1 is needed since in mrtrix3 functions (c++) indices start from 0. I checked multiple voxels to verify this is the case. That is, a 1 should be added to the values in index.data(:,:,:,2)
        inds_start = inds_start(indices_fod);
        
        %-original mask.
        mask = spm_read_vols(h_mask);
        indices = find(mask);
        dim = size(mask);
        
        %-exclude mask voxels with no fixels.
        [~,iNo] = setdiff(indices,indices_fod);
        indices_noFixel = indices(iNo);
        indices(iNo) = [];
        assert(isequal(indices,indices_fod));
        clear indices_fod
        Nv = length(indices);
        
        % update mask: remove voxels with no fixel.
        mask(indices_noFixel) = 0;
        assert(isequal(find(mask),indices));
        %[p,n,e] = fileparts(h_mask.fname);
        h_mask.fname = G.f.mask_afdPruned;
        spm_write_vol(h_mask,mask);
        
        %-VF: fixels in each voxel ========================================
        % Nv x G.N_fixels matrix, in which each row gives the fixel indices
        % associated to a voxel ordered as in 'indices'; rows associated to
        % voxels with less than G.N_fixels fixels have trailing zeros.
        
        VF = zeros(Nv,G.N_fixels);
        for iV=1:Nv
            d1 = inds_start(iV);
            d2 = num_fixels(iV);
            VF(iV,1:d2) = d1:d1+d2-1;
        end
        
        %-determine edges =================================================
        edges = get_edges(ndirs,fixel_dirs,Nf,G.N_fixels);
        [~,~,lni] = get_cinilni_fixel(indices,dirs,dim,VF,num_fixels);
        edges = edges.*lni; % prune out invalid fixel neighbs at mask borders
        
        opts = struct;
        opts.G = G;
        opts.VF = VF;
        opts.mask = mask;
        opts.edges = edges;
        opts.h_mask = h_mask;
        opts.num_fixels = num_fixels;
        opts.fixel_dirs = fixel_dirs;
        opts = update_all(opts);
        
        while ~isempty(opts.fixels_noNeighb)
            Nf = opts.Nf;
            VF = opts.VF;
            indices = opts.indices;
            num_fixels = opts.num_fixels;
            fixel_dirs = opts.fixel_dirs;
            edges = get_edges(ndirs,fixel_dirs,Nf,G.N_fixels);
            [ci,ni,lni] = get_cinilni_fixel(indices,dirs,dim,VF,num_fixels);
            edges = edges.*lni; % prune out invalid fixel neighbs associated to removed voxels
            opts.edges = edges;
            opts = update_all(opts);
        end
        
        %-determine edge weights ==========================================
        switch G.weightingType
            case 'none'
                weights = ones(Nn*G.N_fixels,Nf); %.*logical(edges);
            case 'afd'
                weights = repelem(afds.data',Nn*G.N_fixels,1); %.*logical(edges);
            otherwise
                error('extend..');
        end
        weights = weights.*logical(edges);
        
        assert(all(sum(logical(edges))==sum(logical(weights))));
        
        %-Build adjacency matrix ==========================================
        cilni = ci(lni);
        nilni = ni(lni);
        E = sparse(cilni,nilni,edges(lni),Nf,Nf);
        E = E+E';
        W = sparse(cilni,nilni,weights(lni),Nf,Nf);
        W = W.*W';
        A = W.*E;
        clear A W
        
        %-Update G.
        G.N = Nf;
        G.N_voxels = length(indices);
        G.indices = indices;
        G.A = A;
        G.private.directions_neighb = ndirs;
        
        d2 = inds_start(indices_fod);
        
        %... = flip(fai_inv,2);
        
        
        
        %Pdiff: should be 98xN matrix where there are 1/2 ones in each
        %column and the rest of the elements are zeros. The ones indicate
        %the orientation of the fixel.
        %
        % Note that there will also need to be an adaptation since there
        % are typically more than 1 fixel in a voxel. So the dimension N in
        % Pdiff will specify the total number of fixels in the mask rather
        % than the number of voxels in the mask.
        %
        % Note0: each fixel will be assigned to maximum 2 neighboring voxels, where the two voxels are symmetrical with respect to the focal voxels; Fixel in focal voxels that are at the borders of teh mask will have potentially only one neighboring voxel assigned to them.
        % Note1: in short, note 1 means that a fixel will be projected to 1/2 voxels (al fixels within those voxels).
        % Note2: despite Note1, a voxel (and fixels within it) can 'receive' more connections than the maximum 2 specified in Note1.
        % Note3: the maximum fixel-fixel connection weight should be 2: 1+1. But a fixel can have 'degree' which can be much larger, since many other fixels can be projecting towards it.
end
end

%==========================================================================
function [dirs,ndirs,norms,npar] = ml_create_neighborhood(dim,rempar)
% ML_CREATE_NEIGHBORHOOD Creates direction vectors for a 3D neighborhood.
%   [dirs,ndirs,norms] = ML_CREATE_NEIGHBORHOOD(dim)
%       dim specifies the size of the neighborhood cube, eg. dim=3 would
%       define a 26 neighborhood.
%   [dirs,ndirs,norms] = ML_CREATE_NEIGHBORHOOD(dim,rempar)
%       rempar specifies whether or not to remove parallel vectors, eg.
%       ML_CREATE_NEIGHBORHOOD(5,true) would return 98 and not 124 vectors.
%
%   Author:
%       Martin Larsson
%       June 2017

r = (dim-1)/2;
[y,x,z] = meshgrid(-r:r,-r:r,-r:r);
not_origin = true(1,numel(x));
not_origin((numel(x)+1)/2) = false;
dirs = [x(not_origin); y(not_origin); z(not_origin)];

if nargin > 1 && rempar || nargout > 1
    % Normalized vectors.
    norms = sqrt(sum(dirs.^2));
    ndirs = dirs./repmat(norms,size(dirs,1),1);
end

if nargout > 3 || nargin > 1 && rempar
    % Find parallel vectors.
    same = ndirs'*ndirs > 1-1e-6;
    include = false(size(dirs,2),1);
    for i=1:size(dirs,2)
        include(i) = norms(i) <= min(norms(same(i,:)));
    end
    npar = sum(include);
end

if nargin > 1 && rempar
    % Remove parallel vectors.
    dirs = dirs(:,include);
    ndirs = ndirs(:,include);
    norms = norms(include);
end
end

%==========================================================================
function [NC,h,v] = getNC(G)

h = struct;
v = struct;

% voxel (graph node) class types
h.gm_cmn    = spm_vol(G.f.mask_gm_ctx_subctx_common);
h.gm_subctx = spm_vol(G.f.mask_gm_subctx);
h.gm_ctx    = spm_vol(G.f.mask_gm_ctx);
h.wm        = spm_vol(G.f.mask_wm);

v.gm_cmn    = spm_read_vols(h.gm_cmn);
v.gm_subctx = spm_read_vols(h.gm_subctx);
v.gm_ctx    = spm_read_vols(h.gm_ctx);
v.wm        = spm_read_vols(h.wm);

indices = find(spm_read_vols(spm_vol(G.f.mask)));

assert(isequal(...
    length(indices),...
    nnz(v.wm)+nnz(v.gm_subctx)+nnz(v.gm_ctx)-nnz(v.gm_cmn)...
    ),'inconsistent tissue masks');
d = spm_vol(G.f.mask);
NC = zeros(d.dim);
NC(logical(v.gm_ctx))    = 1; % GM(ctx)
NC(logical(v.gm_subctx)) = 2; % GM(sub-ctx)
NC(logical(v.gm_cmn))    = 2; % label common voxels in ctx & sub-ctx as sub-ctx
NC(logical(v.wm))        = 3; % WM
NC = NC(indices);
end

%==========================================================================
function opts = update_all(opts)

% We remove fixels pointing towards a non-valid voxel. A non-valid voxel is
% either a voxel that is not within the mask, or it is a voxel within the
% mask that has no fixels in it (from begining or due to removing some
% invalid fixels). Note that other fixels might be pointing towards these
% fixels, but since they don't point towards any valid direction, we remove
% them. We then also remove them as destination points to other fixels,
% since we are sure they are not pointing towards any of those fixels, or
% better said, towards the fixels associated to those fixels. Based on
% these removals, we accordingly adjust the files: VF, edges, mask,
% num_fixels.

G = opts.G;
VF = opts.VF;
mask = opts.mask;
edges = opts.edges;
h_mask = opts.h_mask;
num_fixels = opts.num_fixels;
fixel_dirs = opts.fixel_dirs;

Nf = sum(num_fixels);
indices = find(mask);
[Nv,N_fixels] = size(VF);

d = not(logical(sum(edges)));
fixels_noNeighb = find(d);
Nremove = length(fixels_noNeighb);
Nf = Nf-Nremove;

%-update VF ===============================================================
% step1. mark removed fixels as '-1'
Nminus1 = nnz(VF==-1);
d = VF(:);
[d1,d2] = sort(d);
d3 = find(d1>0,1); % 1st non 0 element
d4 = d(d2(d3:end));
assert(isequal(d4(fixels_noNeighb),fixels_noNeighb'));
d4(fixels_noNeighb) = -1; % mark removed fixels -1
d(d2(d3:end)) = d4;

% step2. update fixel indices so max index = Nf not Nf_initial
[d1,d2] = sort(d);
d3 = find(d1>0,1); % 1st element not -1 or 0
d(d2(d3:end)) = 1:Nf;
indices_fixels_mrtrix_order = d1(d3:end);
VF = reshape(d,size(VF));
assert(nnz(VF==-1)==(Nminus1+Nremove));
assert(nnz(VF==0)==(Nv*N_fixels-Nf-Nremove-Nminus1));
assert(isequal(sort(unique(VF(:)')),[-1 0 1:Nf]));
clear d1 d2 d3 d4

%-remove voxels with no valid* fixels left ================================
% *valid fixel: fixel orientd towards a valid neighb voxel
% updtae: num_fixels, indices, VF, mask

num_fixels = sum(VF>0,2); % not 0 or -1
iNo = find(num_fixels==0); % voxels with no valid fixels
indices_noValidFixel = indices(iNo);
num_fixels(iNo) = [];
indices(iNo) = [];
VF(iNo,:) = [];
mask(indices_noValidFixel) = 0;
assert(isequal(find(mask),indices));
[p,n,e] = fileparts(h_mask.fname);
h_mask.fname = fullfile(p,[n,'.fixelOrientationPruned',e]);
%G.f.mask_fixelOrientationPruned = h_mask.fname; %specified in hb_get_G.m
spm_write_vol(h_mask,mask);

assert(Nf==sum(num_fixels));

% update fixel dirs
fixel_dirs(fixels_noNeighb,:) = [];
assert(size(fixel_dirs,1)==Nf);

opts.G = G;
opts.VF = VF;
opts.mask = mask;
opts.edges = edges;
opts.h_mask = h_mask;
opts.indices = indices;
opts.num_fixels = num_fixels;
opts.fixel_dirs = fixel_dirs;
opts.fixels_noNeighb = fixels_noNeighb;
opts.Nremove = Nremove;
opts.Nf = Nf;
end

%==========================================================================
function edges = get_edges(ndirs,directions,Nf,N_fixels)
Nn = size(ndirs,2);

% 2 neighb voxels for each fixel
angles = ndirs'*directions';
M = max(angles); %NnxNf
edges = repmat(M,Nn,1)-abs(angles)<1e-12;
assert(isequal(sum(edges,1),2*ones(1,Nf))); % 2 colinear directions for each fixel

% voxel res to fixel res
edges = repelem(edges,N_fixels,1); % 2xG.N_fixels 1s in each column
end

%==========================================================================
function [ci,ni,lni] = get_cinilni_voxel(indices,dirs,dim,Nn)
Nv = length(indices);
indices_inv = zeros(dim);
indices_inv(indices) = 1:Nv;
ci = repelem(indices,Nn,1);
cs = zeros(size(ci,1),3);
[cs(:,1),cs(:,2),cs(:,3)] = ind2sub(dim,ci);
ns = cs+repmat(dirs',Nv,1);
ni = sub2ind(dim,ns(:,1),ns(:,2),ns(:,3));
ni = reshape(ni,Nn,Nv);
ci = indices_inv(ci);
ni = indices_inv(ni);
lni = logical(ni);
end

%==========================================================================
function [ci,ni,lni] = get_cinilni_fixel(indices,dirs,dim,VF,num_fixels)
Nn = size(dirs,2);
Nv = length(indices);
Nf = sum(num_fixels);
N_fixels = size(VF,2);
assert(Nf==nnz(VF>0));

[~,ni] = get_cinilni_voxel(indices,dirs,dim,Nn);

% replicate voxels with max number of possible fixels in them.
ni = repelem(ni,N_fixels,1);

% change voxel orders to fixel orders
for iV=1:Nv
    [d1,d2] = unique(ni(:,iV));
    for k=1:length(d1)
        if d1(k)~=0
            ni(d2(k):d2(k)+N_fixels-1,iV) = VF(d1(k),:);
        end
    end
end

if any(VF(:)==-1)
    ni(ni==-1) = 0;
    minus1 = true;
else
    minus1 = false;
end

% replicate columns [replicate each voxel by # of fixels in it]
ni = repelem(ni,ones(1,Nn*N_fixels),num_fixels');

assert(size(ni,2)==sum(num_fixels));

% sort columns so that column N is associated with N'th MRtrix3 fixel order
d = VF';
assert(size(d,2)==Nv);
d = d(:);
d(d==0) = [];
if minus1
    d(d==-1) = [];
end
assert(length(d)==Nf);
assert(max(d)==Nf);
[d1,d2] = sort(d);
assert(isequal(d1,(1:Nf)'));
assert(isequal([min(d2) max(d2)],[1 Nf]));
ni = ni(:,d2);

lni = logical(ni);

ci = repelem(1:Nf,Nn*N_fixels,1);

% eventual sparse A matrix rows
assert(min(ci(:))==1);
assert(max(ci(:))==Nf);

% eventual sparse A matrix columns
assert(min(ni(:))==0);
assert(min(ni(ni>0))==1);
assert(max(ni(:))==Nf);

end

%==========================================================================
function [fods,afds,index,directions] = load_fod_etc(G)
dir_fixels = G.f.dir_fixel;
f_fod = fullfile(G.f.diffusion_save,'fods.mif');
f_afd = fullfile(dir_fixels,'afd.mif');
%f_afd3d = fullfile(dir_fixels,'afd_3d.mif');
f_index = fullfile(dir_fixels,'index.mif');
f_directions = fullfile(dir_fixels,'directions.mif');

fods = read_mrtrix(f_fod);
afds = read_mrtrix(f_afd);
%afds3d = read_mrtrix(f_afd3d);
index = read_mrtrix(f_index);
directions = read_mrtrix(f_directions);
end

%==========================================================================
function debug_chk1(mu,a,b,ew,ew2)
x = 0:0.001:1;
figure
subplot(311)
plot(x, mu(x,a,b))
title('Soft threshold function')

subplot(312)
hold on
bin_edges = 0:0.005:1;
bin_edges2 = 0.85:0.001:1;
histogram(ew(ew>0.01),bin_edges)
title('histogram(edgeWeights) before-after comparison (peak around 0 removed for clarity)')

subplot(313)
hold on
histogram(ew(ew>0.85),bin_edges2)
title(sprintf('histogram(edgeWeights(edgeWeights>%.2f)) before-after comparison',a));

subplot(312)
histogram(ew2(ew2>0.01),bin_edges)

subplot(313)
histogram(ew2(ew2>0.85),bin_edges2)
end

%==========================================================================
function debug_chk2(f,mask,indices,A,M,AM,fas)
h = spm_vol(f);
v = spm_read_vols(h);

d1 = mask-v;
d2 = v-mask;
d3 = v+mask;

v0 = zeros(size(v));
v0(d1==1) = 1; % only in mask (WM)
v0(d2==1) = 2; % only in aparaseg_wm
v0(d3==2) = 3; % common between the two

hh = h;
hh.fname = fullfile(fileparts(h.fname),'check___aparcaseg_wm_vs_graphNodes.nii');
spm_write_vol(hh,v0);
clear v0 d1 d2 d3
indwm = find(v);

indwm_inA = ismember(indices,indwm);

% only inter-WM connections
A_wm = A(indwm_inA,indwm_inA);
M_wm = M(indwm_inA,indwm_inA);
AM_wm = AM(indwm_inA,indwm_inA);

% only WM-GM/CSF connections
A_wmElse = A(indwm_inA,~indwm_inA);
M_wmElse = M(indwm_inA,~indwm_inA);
AM_wmElse = AM(indwm_inA,~indwm_inA);

% only inter-GM/CSF connections
A_elseElse = A(~indwm_inA,~indwm_inA);
M_elseElse = M(~indwm_inA,~indwm_inA);
AM_elseElse = AM(~indwm_inA,~indwm_inA);

%%
QAmax = max(fas);
bw = 0.01;

figure(10);
clf(10);
h1 = histcounts(fas(indwm_inA),'BinWidth',bw,'BinLimits',[0 QAmax]);
h2 = histcounts(fas(~indwm_inA),'BinWidth',bw,'BinLimits',[0 QAmax]);
d = 0+bw/2:bw:QAmax;
bar(d,[h1(:)/sum(h1) h2(:)/sum(h2)])
set(gca,'xlim',[0 1]);
legend('WM voxels','GM/CSF voxels')
title('QA')
ylabel('fraction of total # of voxels')
set(gcf,'position',[10 600 750 200]);


binMAX = full(max(max([A(:);M(:);AM(:)]),1));

figure(11);
clf(11);

t = tiledlayout(3,1);
t.TileSpacing = 'compact';

nexttile
aux_debug(A,A_wm,A_elseElse,A_wmElse,binMAX);
title('edges defined only by ODF values (David''s work)');
nexttile
aux_debug(M,M_wm,M_elseElse,M_wmElse,binMAX);
title('the magnitude term in Anjali''s approach')
nexttile
aux_debug(AM,AM_wm,AM_elseElse,AM_wmElse,binMAX);
title('edges defined by QA values & ODF values (anjali''s work, though ODF power replaced by A matrix threshoolding)');

set(gcf,'Name', 'Distribution of edge weights','position',[800 50 1000 1200]);%get(0, 'screensize'));

end

%==========================================================================
function aux_debug(P1,P2,P3,P4,binMAX)

bw = 0.025;
tt = bw/2:bw:binMAX;
cl = lines(4);

d1 = P1(P1(:)>bw); % skip showing first bin
d2 = P2(P2(:)>bw);
d3 = P3(P3(:)>bw);
d4 = P4(P4(:)>bw);

h1 = histcounts(d1,'BinWidth',bw,'BinLimits',[0 binMAX]);
h2 = histcounts(d2,'BinWidth',bw,'BinLimits',[0 binMAX]);
h3 = histcounts(d3,'BinWidth',bw,'BinLimits',[0 binMAX]);
h4 = histcounts(d4,'BinWidth',bw,'BinLimits',[0 binMAX]);

hbar = bar(tt,[h1(:)/sum(h1) h2(:)/sum(h2) h3(:)/sum(h3) h4(:)/sum(h4)]);
hbar(1).FaceColor = cl(1,:);
hbar(2).FaceColor = cl(2,:);
hbar(3).FaceColor = cl(3,:);
hbar(4).FaceColor = cl(4,:);

hbar(1).DisplayName = 'all edges';
hbar(2).DisplayName = 'WM-to-WM edges';
hbar(3).DisplayName = 'else-to-else edges';
hbar(4).DisplayName = 'WM-to-else edges';

ylabel(gca,'distribution of edge weights');

set(gca,'ylim',[0 1]);

yyaxis right
hold on
plot(tt,cumsum(h1/length(d1)),'-o','Color',cl(1,:),'MarkerFaceColor',cl(1,:),'MarkerSize',3,'HandleVisibility','off');
plot(tt,cumsum(h2/length(d2)),'-o','Color',cl(2,:),'MarkerFaceColor',cl(2,:),'MarkerSize',3,'HandleVisibility','off');
plot(tt,cumsum(h3/length(d3)),'-o','Color',cl(3,:),'MarkerFaceColor',cl(3,:),'MarkerSize',3,'HandleVisibility','off');
plot(tt,cumsum(h4/length(d4)),'-o','Color',cl(4,:),'MarkerFaceColor',cl(4,:),'MarkerSize',3,'HandleVisibility','off');
ylabel(gca,'cumulative sum of distributions');
set(gca,'xlim',[0 1])
set(gca,'ylim',[0 1])
set(gca,'xgrid','on');

legend('location','n');

end

% NOTE 1.
% Note the 2nd input. For DSI_studio output, flip is needed. For outputs
% from other software packages, it should be verified if flip across any
% direction is needed. To do: code a systematic check, some how.
%

% NOTE 2.
% omega: solid angle of cone
% theta: (apex angle of cone)/2
% Eq(1): omega = 2*pi(1-cos(theta)).
% Eq(2): cos(angle bw. vec1 & vect2) = vec1'*vec2/(norm(vec1)*norm(vec2))
% Note that the samller the angle between two vects, the larger the cosine
% of the angle.
% Thus, assuming that vec1 is one of the ndirs, i.e., a neighborhood
% directions, and vec2 is one of the vertices, i.e.  an ODF acquisition
% direction, it can be stated that if vec1'*vec2 > cos(theta), vec2 falls
% within the apex angle of the cone, where cos(theta) gives the minimum
% possible cosine of the angle between the vectors, such that the odf angle
% fales within the cone.

% NOTE 3. [April 2020]
% why divide by number of odf dirs in each neighb dir? Shouldn't it be
% instead divide by the number of neighb dirs each odf dir falls in? i.e.,
% if an odf dir falls within the solid angles of two neighb dirs, split the
% value by two.
%
% [21 May 2020]: Now that I think more, it's maybe correct to have this
% sort of normalization, to reduce bias caused by having different number
% of odf samples fall within the different solid angles. If we had instead
% integrated the ODFs, for instance using their spherical harmoic
% representation, an equal solid angle would have been integrated along
% each direction. But since we are implementing this in a discrete way,
% this normalization sort of ensures equal contribution of odf samples
% along the different directions.