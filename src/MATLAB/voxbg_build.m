function voxbg_build(ID,gtype,opts)
% VOXBG_BUILD Builds a CHC or WM voxBG; specifically for HCP data.
%
% CHC  : cerebral hemisphere cortex
% WM   : white matter
% voxBG: voxel-wise brain graph 
% HCP  : Human Connectome Project
%
% Inputs 
%   ID: HCP ubject ID
%   gtype: voxBG type
%
% Hamid Behjat

opts.parallel_pruning = true;  % run graph pruning in parallel? [for CHC]
opts.owrt             = false; % overwrite existing file?

%-Initiate G structure.
G = voxbg_get_struct(ID,gtype,opts);

%-Extract mask.
f_src = pvt_prepare_source(G);

switch G.tissue
    case {'gmlh','gmrh'}
        pvt_tissue_from_ribbonORaparcaseg(...
            f_src,...
            struct('tissue',G.tissue,'res',G.res),...
            opts.owrt,...
            G.f.mask);
    case 'wm'
        switch lower(G.maskType)
            case 'onlywm' % as in Neuroimage-2021-DSS, from aparc+aseg.
                hb_gunzip(G.f.aparcaseg);
                pvt_tissue_from_ribbonORaparcaseg(...
                    G.f.aparcaseg,...
                    struct('tissue',G.tissue,'ResliceRef',f_src),...
                    opts.owrt,...
                    G.f.mask);
            case 'subcortical'
                G = pvt_tissue_from_surface(...
                    f_src,...
                    G,...
                    opts.parallel_pruning,...
                    opts.owrt,...
                    G.f.mask);
        end
end

%-Reslice t1w, t2w & aparc+aseg to match mask.
d = '_acpc_dc_restore_brain.nii';
t1w   = hb_gunzip(fullfile(G.f.T1w,['T1w',d]),G.f.T1w_save);
t2w   = hb_gunzip(fullfile(G.f.T1w,['T2w',d]),G.f.T1w_save);
apas  = hb_gunzip(fullfile(G.f.T1w,'aparc+aseg.nii'),G.f.T1w_save);

delete(t1w);
delete(t2w);
delete(apas);

switch G.tissue
    case 'wm'
        t1wg  = hb_reslice_vol(t1w,G.f.mask,1,G.f.t1w_graphspace);
        t2wg  = hb_reslice_vol(t2w,G.f.mask,1,G.f.t2w_graphspace);
        apasg = hb_reslice_vol(apas,G.f.mask,0,G.f.aparcaseg_graphspace);
        
        gzip(t1wg);
        gzip(t2wg);
        gzip(apasg);
        
        delete(t1wg);
        delete(t2wg);
        delete(apasg);
    case {'gmlh','gmrh'}
        % t1w, t2w, and apas already in graph space.
end

d = spm_vol(hb_gunzip(G.f.mask));
G.dim = d.dim;

%-Build graph.
switch G.tissue
    case {'gmlh','gmrh'}
        G = voxbg_adjacencymatrix_gmxx(G,opts);
    case 'wm'
        G = voxbg_adjacencymatrix_diffusion(G,opts,opts.owrt);
end

%-Save graph.
save(G.fname,'G','-v7.3');

gzip(G.f.mask);
delete(G.f.mask);
end

%==========================================================================
function d = pvt_prepare_source(G)
switch G.space
    case 'T1w'
        d = hb_gunzip(G.f.source,G.f.T1w_save); % ribbon.nii
    case 'Diffusion'
        switch G.tissue
            case {'wmlh','wmrh','wm'}
                d = hb_gunzip(G.f.source,G.f.diffusion_save); % nodif_brain_mask.nii
            case 'cerebrum'
                d = cell(1,2);
                assert(endsWith(fileparts(G.f.source1),'T1w')); % ribbon.nii
                d{1} = hb_gunzip(G.f.source1,G.f.T1w_save);
                assert(endsWith(fileparts(G.f.source2),'Diffusion')); % nodif_brain_mask.nii
                d{2} = hb_gunzip(G.f.source2,G.f.diffusion_save);
        end
end
switch G.tissue
    case {'wmlh','wmrh','wm','cerebrum'}
        hb_gunzip(G.f.aparcaseg,G.f.T1w_save);
end
end

%==========================================================================
function pvt_tissue_from_ribbonORaparcaseg(f_base,opts,chk,f_o)
% Extracting tissue mask, and then resample/reslice if necessary.
% The input file is either of the following two: 
% - the T1w/ribbon.nii [for extarcting gmlh, gmrh, or subcortical mask, i.e., WM and whatever enbounded within white/gray boundary]
% - T1w/aparc+aseg.nii [for extrating subcortical WM, not all subcortical]

if ~chk && exist(f_o,'file')
    return;
end

f_tmp = fullfile(fileparts(f_o),'temp.nii');

%-Extract tissue mask.
h_base = spm_vol(f_base);
v = spm_read_vols(h_base);
switch opts.tissue
    case 'gmlh' % lh cerebral cortex
        t = 3;
    case 'gmrh' % rh cerebral cortex
        t = 42;
    case 'gm'   % cerebral cortex
        t = [3,42];
    case 'wmlh'
        t = 2; % note: I am sipping CC region; see note in line above.
    case 'wmrh'
        t = 41; % note: also skipping CC here. [14 Sep 2020]
    case 'wm'
        t = [2 41 251:255]; % 2: wmlh, 41: wmrh, [251:255]: CC.
    case 'wm_david&martin'
        t = [2 41 251:255 86:88]; % 2: wmlh, 41: wmrh, [251:255]: CC. [86:88]: I don't know what these three regions are. The division of the 5 parcels of CC is anterior-poeterior, so cannot split between Left/Right hemispheres directly.
end

v = ismember(v,t);
h0 = struct();
h0.fname = f_tmp;
h0.dim = h_base.dim;
h0.dt = [spm_type('uint8') h_base.dt(2)]; % See NOTE 2.
h0.mat = h_base.mat;
h0.pinfo = [h_base.pinfo(1:2);0]; % See NOTE 3.
h0 = spm_create_vol(h0);
spm_write_vol(h0,v);
clear v;

%-Resample/Reslice.
interp = 0;
switch opts.tissue
    case {'gm','gmlh','gmrh'}
        if isfield(opts,'ResliceToMatchRef')
            ResliceToMatchRef = opts.ResliceToMatchRef;
        else
            ResliceToMatchRef = false;
        end
        if ResliceToMatchRef % for 'cerebrum' graph design
            f_tmp = hb_reslice_vol(f_tmp,opts.ResliceRef,interp,f_tmp);
        else
            res_o = str2double(opts.res)/1e3;
            d = abs(diag(h_base.mat));
            if isequal(d(1),d(2),d(3))
                res_i = num2str(d(1));
            else
                res_i = 'n/a';
            end
            if res_o>=res_i
                apprch = 'approach1';
            else
                apprch = 'approach2';
            end
            if ~strcmp(res_i,num2str(res_o))
                f_tmp = hb_nii_resample(...
                    f_tmp,...
                    res_o,...
                    'InterpOrder',interp,...
                    'OutputFile',f_tmp,...
                    'Method',apprch);
            end
        end
    case {'wmlh','wmrh','wm'}
        f_tmp = hb_reslice_vol(f_tmp,opts.ResliceRef,interp,f_tmp);
end

%-Threshold & Binarize.
% -If f_tmp is not resampled/resliced, it's already binary.
% -If f_tmp is resampled/resliced, it's aready thresholded & binarized,
%  since we set it's header h.dt=2, i.e., 'uint8' type; See NOTE 4.
%  MOreover, we force this by ussing interp=0. 

%-Make single connected.
h = spm_vol(f_tmp);
[v,d1,d2] = hb_make_connected(h,26);
fprintf('\n..Mask made connected:');
fprintf('\n  number of voxels removed: %d (%f %%)',d1,d2*1e2);
h.fname = f_o; 
spm_write_vol(h,v);

delete(f_tmp);
end

%==========================================================================
function G = pvt_tissue_from_surface(f_mask,G,runPar,owrt,f_o)
% Extract white structure from f_mask by keeping only those voxels that are
% within the white surface.

if ~owrt
    if exist(f_o,'file')
        return;
    elseif exist([f_o,'.gz'],'file')
        hb_gunzip(f_o);
        return;
    end
end

h_mask = spm_vol(f_mask);
mask = spm_read_vols(h_mask);
mask = hb_make_connected(mask,26);
indices = find(mask);

% adjacencies across whole brain mask
A0 = hb_compute_adjacency(mask,26,'sS',0);

% remove edges passing through white surface
whites = ml_batch_gifti(G.f.surface.white,h_mask.mat);
A = ml_prune_adjacency(A0,mask,whites,'parallelize',runPar);

% get connected comps
bins = conncomp(graph(A));

% sort comps based on size
bin_count = max(bins);
bins_size = histcounts(bins,1:bin_count+1);

% remove comps other than the desired white structure (lh, rh, or both)
switch G.tissue
    case {'wmlh','wmrh'}
        Ncomp = 1; 
    case {'wm','cerebrum'}
        Ncomp = 2; 
end
check_comps(bins_size,Ncomp);
[~,d] = sort(bins_size,'descend');
switch Ncomp
    case 1
        d = ismember(bins,d(2));
    case 2
        d = ismember(bins,d(2:3));
end
inds_rm = indices(~d); % voxels to remove
mask_marked = mask;  
mask_marked(inds_rm) = 2; % volume showing removed voxels
mask(inds_rm) = 0;

% After pruning out edges that pass through the white surfaces on the
% nodif_brain_mask, we use the aparc+aseg, to adjust the resulting mask by,
%  adding: corpus collaum 
% droping: ventricles + CSF.

switch G.tissue
    case 'cerebrum'
        f_src = G.f.source2;
    otherwise
        f_src = G.f.source;
end

switch G.tissue
    
    case {'wm','cerebrum'}
        
        switch G.maskType
                
            case {'subcortical','cerebrum'}
                
                %-Get CC (to add) and ventricles+CSF (to omit).
                h_apas = spm_vol(G.f.aparcaseg);
                v = spm_read_vols(h_apas);
                t = 251:255; % CC
                v_add = ismember(v,t);
                t = [4,43,14,15,72,24]; % ventricles:4,43,14,15,72 - CSF:24 
                v_omt = ismember(v,t);
                
                %-Resample/Reslice to match mask. 
                h = struct();
                h.dim = h_apas.dim;
                h.mat = h_apas.mat;
                h.dt = [spm_type('uint8') 0]; 
                f_add = fullfile(fileparts(f_o),'temp1.nii');
                if exist(f_add,'file')
                    delete(f_add);
                end
                h.fname = f_add;
                spm_write_vol(h,v_add);
                f_add = hb_reslice_vol(f_add,f_src,1,f_add);
                f_omt = fullfile(fileparts(f_o),'temp2.nii');
                if exist(f_omt,'file')
                    delete(f_omt);
                end
                h.fname = f_omt;
                spm_write_vol(h,v_omt);
                f_omt = hb_reslice_vol(f_omt,f_src,1,f_omt);
                
                %-Update mask: add CC and omit ventricles+CSF.
                v_add = logical(spm_read_vols(spm_vol(f_add)));
                v_omt = logical(spm_read_vols(spm_vol(f_omt)));
                mask  = logical(mask);
                mask = or(mask,v_add);       % in mask or CC
                mask = and(mask,not(v_omt)); % in mask but not ventricle or CSF
                
                %-Clean-up.
                delete(f_add);
                delete(f_omt);
                
            otherwise
                
        end
        
    otherwise
        
end

%-Make mask connected.
d1 = nnz(mask);
mask = hb_make_connected(mask,26);
d2 = nnz(mask);
assert(d2>.95*d1,'fishy; too many voxels removed; check resulting mask');

%-Write mask.
h = struct();
h.fname = f_o;
h.dim = h_mask.dim;
h.dt = [spm_type('uint8') 0];
h.mat = h_mask.mat;
%h.pinfo = [h_mask.pinfo(1:2);0]; % See NOTE 3. - commented out on 19.03.2022
spm_write_vol(h,mask);

%-Write base mask with removed parts marked.
[p,n,e] = fileparts(f_o);
h.fname = fullfile(p,[n,'.base_mask_with_removed_part_marked',e]);
spm_write_vol(h,mask_marked);

gzip(h.fname);
delete(h.fname);

%-Update G.
G.pruning.A_diff_pre_post_white_pruning = A0-A;
G.pruning.mask_with_removed_part_marked = h.fname;
G.pruning.info = fullfile(p,[n,'.extracting_white_structure_info.txt']);

%-Write info file.
fid = fopen(G.pruning.info,'wt');
fprintf(fid,'=========================================================\n');
fprintf(fid,'%%To plot prunned edges resulting in detached white: \n');
fprintf(fid,'%%A0: pre-pruning A. \n');
fprintf(fid,'%%Ap : post-pruning A. \n');
fprintf(fid,'\n');
fprintf(fid,'d = hb_gunzip(G.pruning.mask_with_removed_part_marked); \n');
fprintf(fid,'mask = spm_read_vols(spm_vol(d)); \n');
fprintf(fid,'A0 = hb_compute_adjacency(logical(mask),26,''sS'',0); \n');
fprintf(fid,'Adiff = G.pruning.A_diff_pre_post_white_pruning; \n');
fprintf(fid,'Ap = A0-Adiff;\n');
fprintf(fid,'sl = 88;       %%axial slice number \n');
fprintf(fid,'d1 = ''grid''; %%''diff'',''none'' \n');
fprintf(fid,'d2 = true;     %%false \n');
fprintf(fid,'d3 = ''r'';    %%color \n');
fprintf(fid,'figure; \n');
fprintf(fid,'hm_plot_adjacency_diff(A0,Ap,mask,sl,d1,d2,d3);\n');
fprintf(fid,'%%or \n');
fprintf(fid,'figure; \n');
fprintf(fid,'hm_plot_adjacency_diff(A0,Ap,logical(mask),sl,d1,d2,d3);\n');
fprintf(fid,'=========================================================\n');
fclose(fid);
end

%==========================================================================
function check_comps(s,n)
%CHECK_COMPS assuming that the input mask has been an entire brain mask,
%inlcuding the brain stem and cerebellum, a heuristic check is done on the
%component sizes, s, to ensure that each component is a brain structure as
%follows:
% [Assume M=length(s)]
% If n==1:
% component 1   : the brain mask, excluding the desired l/r cerebral white matter.
% component 2   : the desired l/r cerebral white matter.
% components 3-M: junk stuff which sould be really small.
% If n==2
% component 1   : the  brain mask, excluding the l&r cerebral white matters.
% components 2&3: the l&r cerebral white matter.
% components 4-M: junk stuff which sould be really small.

s = sort(s,'descend');
switch n
    case 1
        sts(1) = s(1)>4*s(2);
        sts(2) = s(3)<(s(2)*1e-3);
    case 2
        sts(1) = s(1)>2*s(2);
        sts(2) = s(1)>2*s(3);
        sts(3) = s(4)<(s(3)*1e-3);
end
assert(all(sts),...
    'fishy; extracted components in brain mask arer not as expected.');
end

%==========================================================================
function extract_cerebrum_masks(f_rb,f_aa,f_rf,f_masks)
%Extracts and writess to disc 3 masks that are required for building a
%cerebrum voxBG.

%-Names of files to be saved.

f1 = f_masks.gm_ctx;
f2 = f_masks.gm_subctx;
f3 = f_masks.wm;
f_gm = f_masks.gm;
f_gm_common = f_masks.gm_ctx_subctx_common;
f_omit = strrep(f3,...
    '.wm.nii',...
    '.discarded_voxels_from_aparaseg_in_defining_mask.nii');
f_mask_all = f_masks.all;

f1_hr = strrep(f1,'.res1250','.res0700'); 
f2_hr = strrep(f2,'.res1250','.res0700'); 
f3_hr = strrep(f3,'.res1250','.res0700'); 
f_gm_hr = strrep(f_gm,'.res1250','.res0700'); 
f_gm_common_hr = strrep(f_gm_common,'.res1250','.res0700'); 
f_omit_hr = strrep(f_omit,'.res1250','.res0700'); 
f_mask_all_hr = strrep(f_mask_all,'.res1250','.res0700'); 

%-Source volumes. 
h_rb = spm_vol(f_rb);
h_aa = spm_vol(f_aa);
v_rb = spm_read_vols(h_rb);
v_aa = spm_read_vols(h_aa);

%-Extract mask1: GM cortex.
v1 = ismember(v_rb,[3,42]);                                                % lh and rh ctx

%-Extract mask2: GM sub cortex.
d1 = [9:13,17:20,26:28, 48:56, 58:60];                                     % sub-cortical GM  
d2 = [80,81,82];                                                            % non-WM hypointensities
v2 = ismember(v_aa,[d1,d2]); 

%-Extract mask3: WM.
d1 = [2,41, 85, 250, 251:255];                                             % Cerebrum WM; see label descriptions below.
d2 = [77,78,79];                                                           % WM hypointensities 
v3 = ismember(v_aa,[d1,d2]);                                     
d1 = [4,5,43,44,14,15,72, 24, 31,63, 30,62];                               % ventricles, CSF, choroid-plexus and vessels.
d2 = [1,6,40,45,21:23,25,57,75,76];                                        % seemingly crap
d3 = [7,8,46,47, 16];                                                      % cerebellum and brain stem (these should not be in subcortical mask extracted from ribbon, but we omit them just in case they are partially there)
v_omit_p1 = ismember(v_aa,[d1,d2,d3]);
v3 = and(v3,not(v_omit_p1));
v3 = and(v3,not(v1));
v3 = and(v3,not(v2));

% What-to-do-with regions. 
assert(nnz(ismember(v_aa,[29,61]))==0,'What to do with these?');           % see label descriptions below
assert(nnz(ismember(v_aa,[32,64]))==0,'What to do with these?');           % see label descriptions below
assert(nnz(ismember(v_aa,[33:39, 65:71]))==0,'What to do with these?');    % see label descriptions below 
assert(nnz(ismember(v_aa,[73,74, 83,84]))==0,'What to do with these?');    % see label descriptions below

%-Adjust mismatch bw. ribbon's and aparc+aseg's definition of cortex.
%
% These are UnAssigned (ua) voxels (a lot!) at the gray/white cerebral
% cortex boundary (CLASS1), also, at the pial/CSF boundary (CLASS2), that
% are assigned as being cortex in aparc+aseg.nii, but are not defined as
% cortex in ribbon.nii. 
% 
% One should note that the cerebral cortex as defined by ribbon is
% more accurate/smooth/neat than that defined by aparc+aseg. So nonw of
% these unassigned voxels will be labelled as cerebral cortex. They will
% either be considered as white matter (CLASS1) or will be dicarded
% (CLASS2).

% First, we assign all CLASS1 and CLASS2 voxels as being WM.
d = [1000:1035, 2000:2035];
v_ua = and(ismember(v_aa,d),not(v1)); 
v3 = or(v3,v_ua);   
if 0 % for debug chk
    figure;  %#ok<UNRCH>
    for k=1:150, imagesc(v_ua(:,:,k)); pause; end
end

% Then, we clean out CLASS2 voxels; these are isolated voxels or little
% patches of voxels, so they can be dicared via aking WM mask connected. 
[v3,~,~,v_omit_p2] = hb_make_connected(v3,26);

%-Union of ctx and subctx gray matter.
% some voxels are present in both masks, e.g., those in Hippocampus.
v_gm = or(v1,v2); 

%-Voxels present in both gm_ctx and gm_sub_ctx
v_gm_common = and(v1,v2); 

if 0 
    % NOTE: better not to remove the common voxels since I do save the common
    % voxels as a seperate mask, v_gm_common. That mask can be loaded* if
    % necessary for some particular analysis, e.g., when wanting to specify
    % classes for EdgeTypeToRecover in the pruning phase in
    % hb_adjacencymatrix_diffusion.m
    %
    % *: G.f.mask_gm_ctx_subctx_common

    %-Remove common GM voxels GM from cortex.
    % [ribbon not reliable in sub-cortical regions,
    % so better to drop for ctx than sub-ctx]
    v1(v_gm_common) = 0;
    assert(nnz(and(v1,v2))==0);
    assert(isequal(v_gm,or(v1,v2))); % no chnage in v_gm
end

% All voxels in aparc+aseg that were omitted.
v_omit = or(v_omit_p1,v_omit_p2);

% Brain graph mask: wm and gm (ctx + subctx).
v_all = or(v_gm,v3);

%-Write highres masks.
h = struct();
h.dim = h_aa.dim;
h.mat = h_aa.mat;
h.dt = [spm_type('uint8'),0]; 

h.fname = f1_hr;
spm_write_vol(h,v1);

h.fname = f2_hr;
spm_write_vol(h,v2);

h.fname = f3_hr;
spm_write_vol(h,v3);

h.fname = f_gm_hr;
spm_write_vol(h,v_gm);

h.fname = f_gm_common_hr;
spm_write_vol(h,v_gm_common);

h.fname = f_omit_hr;
spm_write_vol(h,v_omit);

h.fname = f_mask_all_hr;
spm_write_vol(h,v_all);

clear v_rb v_aa v1 v2 v3 v_omit_p1 v_omit_p2 v_omit v_gm v_gm_common v_all;

%-Reslice hr mask to diffusion space.
hb_reslice_vol(f1_hr,f_rf,0,f1,true);
hb_reslice_vol(f2_hr,f_rf,0,f2,true);
hb_reslice_vol(f3_hr,f_rf,0,f3,true);
hb_reslice_vol(f_gm_hr,f_rf,0,f_gm,true);
hb_reslice_vol(f_gm_common_hr,f_rf,0,f_gm_common,true);
hb_reslice_vol(f_omit_hr,f_rf,0,f_omit,true);
hb_reslice_vol(f_mask_all_hr,f_rf,0,f_mask_all,true);

%-Remove voxels that don't have diffusion data.
% (i.e., voxels that are not in nodif_brain_mask.nii)
h = spm_vol(f_mask_all);
v_mask = spm_read_vols(h);
v_nodif = spm_read_vols(spm_vol(f_rf));
inds_mask = find(v_mask);
inds_nodif = find(v_nodif);
d = not(ismember(inds_mask,inds_nodif));
v_mask(inds_mask(d)) = 0;
spm_write_vol(h,v_mask);

%-Make main mask connected. 
h = spm_vol(f_mask_all);
[v_mask,d1,d2] = hb_make_connected(h,26);
fprintf('\n..G.f.mask made connected:');
fprintf('\n  number of voxels removed: %d (%f %%)',d1,d2*1e2);
spm_write_vol(h,v_mask);
% remove deleted voxels also from sub-masks:
d = {
    f_masks.gm_ctx
    f_masks.gm_subctx
    f_masks.wm
    f_masks.gm
    f_masks.gm_ctx_subctx_common
    };
for k=1:length(d)
    h = spm_vol(d{k});
    v = and(spm_read_vols(h),v_mask);
    spm_write_vol(h,v);
end

%-GM/WM border masks.
h_gm = spm_vol(f_masks.gm);
h_wm = spm_vol(f_masks.wm);
v_gm = spm_read_vols(h_gm);
v_wm = spm_read_vols(h_wm);
d = ones(3,3,3);
v_gm_wmborder = and(v_gm,imdilate(v_wm,d));
v_wm_gmborder = and(v_wm,imdilate(v_gm,d));
h = h_gm;
h.fname = f_masks.gm_wmborder;
spm_write_vol(h,v_gm_wmborder);
h = h_wm;
h.fname = f_masks.wm_gmborder;
spm_write_vol(h,v_wm_gmborder);

%==========================================================================
% These are all the 89 labels from apar+aseg that are included: 
lbls = [...
    9:13,17:20,26:28, 48:56, 58:60,...
    80,81,82,...
    2,41, 85, 250, 251:255,...
    77,78,79,...
    4,5,43,44,14,15,72, 24, 31,63, 30,62,...
    1,6,40,45,21:23,25,57,75,76,...
    7,8,46,47, 16,...
    29,61,...
    32,64,...
    33:39, 65:71,...
    73, 74, 83, 84];
assert(length(unique(lbls))==89);

%---Cortical GM [in ribbon; in apar+aseg, 3 and 42 replcaed with values 1xxx and 2xxx, respectively.]:
%  3  Left-Cerebral-Cortex
% 42  Right-Cerebral-Cortex 

%---Cerebrum WM: 
%   2  Left-Cerebral-White-Matter
%  41  Right-Cerebral-White-Matter
%  85  Optic-Chiasm               [where optic nerves cross, in WM]
% 250  Fornix                     [missing in 100307]
% 251  CC_Posterior
% 252  CC_Mid_Posterior
% 253  CC_Central
% 254  CC_Mid_Anterior
% 255  CC_Anterior 

%---Sub-cortical GM:
%  9  Left-Thalamus [missing? at least in 100307]                           
% 10  Left-Thalamus-Proper*                   
% 11  Left-Caudate                            
% 12  Left-Putamen                            
% 13  Left-Pallidum                           
% 17  Left-Hippocampus                        
% 18  Left-Amygdala                           
% 19  Left-Insula           [missing in 100307]                            
% 20  Left-Operculum        [missing in 100307]        
% 26  Left-Accumbens-area                     
% 27  Left-Substancia-Nigra [missing in 100307]                   
% 28  Left-VentralDC 
%
% 48  Right-Thalamus                          
% 49  Right-Thalamus-Proper*                  
% 50  Right-Caudate                           
% 51  Right-Putamen                           
% 52  Right-Pallidum                          
% 53  Right-Hippocampus                       
% 54  Right-Amygdala                          
% 55  Right-Insula           [missing in 100307]                           
% 56  Right-Operculum        [missing in 100307]                         
% 58  Right-Accumbens-area                    
% 59  Right-Substancia-Nigra [missing in 100307]                  
% 60  Right-VentralDC
%----------------------

%---Regions to exclude [non-tissue]:
%  4  Left-Lateral-Ventricle
%  5  Left-Inf-Lat-Vent
% 43  Right-Lateral-Ventricle
% 44  Right-Inf-Lat-Vent
% 14  3rd-Ventricle                           
% 15  4th-Ventricle    
% 72  5th-Ventricle 
%
% 24  CSF
%
% 31  Left-choroid-plexus
% 63  Right-choroid-plexus
%
% 30  Left-vessel
% 62  Right-vessel 

%---Regions to exclude [tissue, but not included in 'cerebrum' voxBG]:
%  7  Left-Cerebellum-White-Matter
%  8  Left-Cerebellum-Cortex
% 46  Right-Cerebellum-White-Matter
% 47  Right-Cerebellum-Cortex
% 16  Brain-Stem     

%---Regions to exclude [seem like crap, and probably typically missing on HCP subjects' aseg]: 
%  1  Left-Cerebral-Exterior    [?; missing in 100307]
%  6  Left-Cerebellum-Exterior  [?; missing in 100307]
% 40  Right-Cerebral-Exterior   [?; missing in 100307]
% 45  Right-Cerebellum-Exterior [?; missing in 100307]
% 21  Line-1                    [?; missing in 100307]
% 22  Line-2                    [?; missing in 100307]
% 23  Line-3                    [?; missing in 100307]
% 25  Left-Lesion               [abnormality; missing in 100307]     
% 57  Right-Lesion              [abnormality; missing in 100307]  
% 75/76                         [FS description: removed. duplicates of 4/43; missing in 100307]

%---Regions of unknown nature; but since tissue, I include them in WM and subx-ctx GM:
% 77  WM-hypointensities           [in 100307]
% 78  Left-WM-hypointensities      [missing in 100307]
% 79  Right-WM-hypointensities     [missing in 100307]
% 80  non-WM-hypointensities       [in 100307]
% 81  Left-non-WM-hypointensities  [missing in 100307]
% 82  Right-non-WM-hypointensities [missing in 100307]

%---Regions that are undetermined: [WHAT TO DO with ???]
% 29  Left-undetermined  [missing in 100307]
% 61  Right-undetermined [missing in 100307]

%---Regions that are apparently cerebral cortex, which we will keep: 
% 32  Left-F3orb  [this is seemingly called the pars orbitalis as I read in a paper; missing in 100307]
% 64  Right-F3orb [this is seemingly called the pars orbitalis as I read in a paper; missing in 100307]

%---Regions I don't know what they are or what to do with: [all missing in 100307]
% 33  Left-lOg
% 34  Left-aOg
% 35  Left-mOg
% 36  Left-pOg
% 37  Left-Stellate
% 38  Left-Porg
% 39  Left-Aorg
%
% 65  Right-lOg
% 66  Right-aOg
% 67  Right-mOg
% 68  Right-pOg
% 69  Right-Stellate
% 70  Right-Porg
% 71  Right-Aorg
%
% 73  Left-Interior
% 74  Right-Interior
%
% 83  Left-F1
% 84  Right-F1
%==========================================================================
end

%==========================================================================
%-NOTES.
%--------------------------------------------------------------------------
% NOTE 2 [3,21.04.2020, 19.03.2022]
% This is the data type --- See spm_type.m
% Here, we set it to spm_type('uint8'), i.e., 2, since our file is:
% -for gmlh/rh: a logocal file, storing only 2 numbers: 0 & 1.
% -for wmlh/rh: a file storing either 2 (0,1) or only 3 (0,1,2) numbers.
% Or just a binary file (0 & 1) if we are saving a mask. 
% These few values fall within the accepatable range for 'uint8': 0-2^8-1
% Using this small data type, significantly reduces the file size:
% 0.7 mm^3 res: 84 MB ('float32') > 0.7 mm^3 res: 21 MB ('unit8').
% For instance, the 1 mm^3 res file will become 7 MB ('uint8').
%
% NOTE 3 [4,21.04.2020]
% The third element of pinfo, refers to off set in bytes; for instance, see
% the following link for some further info:
% https://nipy.org/nibabel/devel/spm_use.html
% It is imprtant to change this, since we are saving each plane as a logical
% array and thus essentially there needs not be any offset in bytes.
% If we take the value from h_rb/mask.pinfo, then that can be problematic,
% since that value is potentially associated to the other data types that
% were used for saving those volumes, probably as a float32 data type.
% spm_read_vols, etc. will run into problem when wanting to load the files
% this property is not correctly set. 
% [19.03.2022]: just don't set it, i.e., don't define that field for the
% header; spm_write_vol should properly take care of it when writing the
% volume.

%
% NOTE 4 [05.04.2020]
% If the input volume is 'logical', like the volumes I create here, then
% the output resampled volume is forced to be logical, so no need for
% thresholding. I suppose this is what spm_write_vol.m and spm_write_plane.m
% impose by seeing that the h.dt(1) is 2, i.e., 'uint8': unsigned integer,
% with min-max: 0-2^8-1; see NOTE 2. Thus, values in range [0,1] are
% rounded to closest integer value, i.e., 0 or 1.
