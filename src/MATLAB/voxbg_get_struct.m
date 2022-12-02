% same as hb_get_G but with difference that it does not save; just loads G
% if it exits; sincr the checks to comapre loaded G and created G in
% hb_get_G are not so robust.
% 
% [06.05.2022]
 

function  [G,sts] = hb_get_G_v2(ID,gtype,opts,verbose,justGf)
%HB_GET_G creates G structure. If G exists, it is loaded and verified, and
%if necessary, it's content is updated.
%
%
% Hamid Behjat

if ~exist('justGf','var') || isempty(justGf)
    justGf = false;
end

if ~exist('verbose','var') || isempty(verbose)
    verbose = true;
end

if ~ischar(ID)
    assert(isnumeric(ID),'ID should be either numeric or char')
    ID = num2str(ID);
end

hcp_root= opts.hcp_root;

if ~isfield(opts,'hcpsave_root') || isempty(opts.hcpsave_root) % 20 jan 2021
    hcpsave_root = hcp_root;
    if verbose
        fprintf('\n..Results will be save on same directory as original HCP data.');
    end
else
    hcpsave_root = opts.hcpsave_root; % 14 April 2020. Save MATLAB results on local. Then transfer to remote using rsync.
end

G = struct;
G.type = gtype;
G.subject = ID;
G.fname = [];
d = strfind(gtype,'.res');
G.tissue = gtype(1:d-1);

switch G.tissue
    case {'gmlh','wmlh','pslh','wslh'}
        G.hemisphere = 'left';
    case {'gmrh','wmrh','psrh','wsrh'}
        G.hemisphere = 'right';
    case {'gm','wm','cerebrum'}
        G.hemisphere = 'both';
end

switch G.tissue
    case {'pslh','psrh'}
        G.surftype = 'pial';
    case {'wslh','wsrh'}
        G.surftype = 'white';
    otherwise
        G.surftype = 'n/a';
end

switch G.surftype
    case {'pial','white'}
        if isfield(opts,'weighted_surface_graph')
            G.weighted = opts.weighted_surface_graph;
        else
            G.weighted = false;
        end
end

G.res = gtype(d+(4:7));

if contains(gtype,'T1w')
    G.space = 'T1w';
elseif contains(gtype,'Diffusion')
    G.space = 'Diffusion';
end

G.resTag   = ['.res',G.res];
G.spaceTag = ['.space',G.space];
G.rsTag    = [G.resTag,G.spaceTag];
G.trsTag   = [G.tissue,G.rsTag];

switch G.tissue
    case {'wm','wmlh','wmrh','cerebrum'}
        if contains(gtype,'.fod')
            G.diffusionModel = 'fod';
            G.settingsTag = gtype(strfind(gtype,'.fixels'):end);
            G.N_fixels = str2double(G.settingsTag(strfind(G.settingsTag,'.fixels')+length('.fixels'):strfind(G.settingsTag,'.maxAfd')-1));
            if contains(G.settingsTag,'.weighted')
                d1 = strfind(G.settingsTag,'.weighted');
                d2 = strfind(G.settingsTag,'.onlywm');
                G.weightingType = lower(G.settingsTag(d1+9:d2-1));
                tag = '.weighted';
            elseif contains(G.settingsTag,'.nonweighted')
                G.weightingType = 'none';
                tag = '.nonweighted';
            else
                error('extend..');
            end
            G.maxAfd = (1/100)*str2double(G.settingsTag(strfind(G.settingsTag,'.maxAfd')+length('.maxAfd'):strfind(G.settingsTag,tag)-1));
            G.neighb = 5;
            G.rempar = true;
        else
            G.diffusionModel = 'odf';
            G.settingsTag = gtype(strfind(gtype,'.neighb'):end);
            d = G.settingsTag(strfind(G.settingsTag,'.neighb')+length('.neighb'):strfind(G.settingsTag,'.pow')-1);
            if length(d)==1
                G.neighb = str2double(d);
                G.neighbHybrid = false;
            elseif isequal(d,'5+3')
                G.neighb = 5;
                G.neighbHybrid = true; % WM-WM, WM-GM: neighb5, GM-GM: neighb3
            else
                error('invalid neighborhood connectivity specifier');
            end
            G.odfPow    = str2double(G.settingsTag(strfind(G.settingsTag,'.pow')   +length('.pow')   :strfind(G.settingsTag,'.fibs')-1));
            G.N_fibers  = str2double(G.settingsTag(strfind(G.settingsTag,'.fibs')  +length('.fibs')  :strfind(G.settingsTag,'.rempar')-1));
            G.rempar    = str2double(G.settingsTag(strfind(G.settingsTag,'.rempar')+length('.rempar'):strfind(G.settingsTag,'.sh')-1));
            G.shOrder   = str2double(G.settingsTag(strfind(G.settingsTag,'.sh')    +length('.sh')    :strfind(G.settingsTag,'.th')-1));
            G.thresh    = str2double(G.settingsTag(strfind(G.settingsTag,'.th')    +length('.th')    :strfind(G.settingsTag,'.beta')-1)); % \alpha parameter in the sigmoid function
            G.beta      = str2double(G.settingsTag(strfind(G.settingsTag,'.beta')+length('.beta'):strfind(G.settingsTag,'.mag')-1)); % \beta parameter in the sigmoid function
            G.magnitude = str2double(G.settingsTag(strfind(G.settingsTag,'.mag')+length('.mag'))); %:d(end-1)-1)); % use magnitude term as in Anjali's work?
            d = strfind(G.settingsTag,'.');
            G.method = G.settingsTag(d(end-1)+1:d(end)-1);
            %switch G.tissue
            %    case {'wm','wmlh','wmrh'}
            %G.magnitude = str2double(G.settingsTag(strfind(G.settingsTag,'.mag')+length('.mag'):d(end-1)-1)); % use magnitude term as in Anjali's work?
            %G.method = G.settingsTag(d(end-1)+1:d(end)-1);
            %    case 'cerebrum'
            %G.magnitude = str2double(G.settingsTag(strfind(G.settingsTag,'.mag')+length('.mag'):d(end)-1));
            %G.method = G.settingsTag(d(end)+1:end);
            
            if strcmp(G.tissue,'cerebrum')
                
                d = str2double(G.settingsTag(strfind(G.settingsTag,'.gmEdgeWeightVers')+length('.gmEdgeWeightVers')));
                
                %-Weighting type for GM-GM/WM edges:
                % [v1] weighted based on ODFs, just as for WM-WM edges.
                % [v2] non-weighted:
                %       GM-GM edge weights: 1.
                %       GM-WM edge weights: 1.
                % [v3] hybrid weighting: 
                %       GM-GM edges weights: inverse Euclidean distance. 
                %       GM-WM edge weights: as in v1. 
                
                G.gmEdgeWeightVers = d;
                G.gmEdgeWeightTag = sprintf('.gmEdgeWeightVers%d',d);
            end
            %end
        end
        switch G.tissue
            case {'wm','wmlh','wmrh'}
                d = str2double(G.settingsTag(strfind(G.settingsTag,'.onlywm')+length('.onlywm'))); % only WM as graph nodes or include other subcortical stuff?
                switch d
                    case 1
                        G.maskType    = 'onlyWM';
                        G.maskTypeTag = '.onlywm1';
                    case 0
                        G.maskType    = 'subcortical';
                        G.maskTypeTag = '.onlywm0';
                end
            case 'cerebrum'
                G.maskType    = 'cerebrum';
                G.maskTypeTag = '';
        end
        G.maskTypeTag = sprintf('%s.neighb%d',G.maskTypeTag,G.neighb);
        
    otherwise
        
        G.settingsTag = '';
        switch G.tissue
            case {'pslh','psrh','wslh','wsrh'} % pial/white surface
                G.neighb = 'n/a';
                G.maskTypeTag = '';
            otherwise
                G.neighb = 3;
                G.maskTypeTag = sprintf('.neighb%d',G.neighb);
        end
end

G.f = [];
f = struct; 

%-Directories--------------------------------------------------------------
f.hcp_root = hcp_root;
f.hcpsave_root = hcpsave_root;
f.T1w = fullfile(hcp_root,ID,'T1w');
f.MNI = fullfile(hcp_root,ID,'MNINonLinear');
f.T1w_save         = fullfile(hcpsave_root,ID,'T1w');
f.T1w_results_save = fullfile(f.T1w_save,'Results');
f.MNI_save         = fullfile(hcpsave_root,ID,'MNINonLinear');
f.MNI_results      = fullfile(f.MNI,'Results');
f.MNI_results_save = fullfile(f.MNI_save,'Results');
f.graphmain        = fullfile(f.T1w_save,'graph'); 
f.graph            = fullfile(f.graphmain,strrep(gtype,'.','_'));
switch G.tissue
    case {'pslh','psrh','wslh','wsrh'}
        if G.weighted
            f.graph = [f.graph,'_weighted'];
        end
end
[~,~] = mkdir(f.T1w_save); % outputs set to prevent warnings
[~,~] = mkdir(f.MNI_save);
[~,~] = mkdir(f.MNI_results_save);
[~,~] = mkdir(f.graph);
switch G.tissue
    case {'wm','wmlh','wmrh','cerebrum'}
        f.diffusion        = fullfile(f.T1w,'Diffusion');
        f.diffusion_save   = fullfile(f.T1w_save,'Diffusion');
        [~,~] = mkdir(f.diffusion_save);
end
f.G = fullfile(f.graph,['G.',G.type,'.mat']);

%-Volumetric files---------------------------------------------------------
% f.source is used for mask extraction.
% f.mask is the file used for G design, which has been:
% 1. Extracted from f.source.
% 2. Resampled/resliced to desired resolution/space.
% 3. Made single connected.

%switch G.tissue
%    case {'wm','wmlh','wmrh'}
%        midtag = G.maskTypeTag;
%    otherwise
%        midtag = '';
%end
%n = [G.tissue,midtag,G.resTag]; %gtype(1:(strfind(gtype,'res')+6));
n = [G.tissue,G.maskTypeTag,G.resTag];

switch G.tissue
    case {'gm','gmlh','gmrh'}
        f.source = fullfile(f.T1w,'ribbon.nii');
        f.mask = fullfile(f.graph,[n,'.spaceT1w.nii']);
    case {'wm','wmlh','wmrh','cerebrum'}
        switch G.tissue
            case {'wm','wmlh','wmrh'}
                f.source = fullfile(f.diffusion,'nodif_brain_mask.nii');
                f.mask = fullfile(f.graph,[n,'.DiffusionSpace.nii']);
                if contains(gtype,'.fod')
                    f.mask_afdPruned = strrep(f.mask,'.nii','.afdPruned.nii');
                    f.mask_fixelOrientationPruned = strrep(f.mask,'.nii','.fixelOrientationPruned.nii');
                end
            case 'cerebrum'
                f.source1 = fullfile(f.T1w,'ribbon.nii');                 % to extract cortex
                f.source2 = fullfile(f.diffusion,'nodif_brain_mask.nii'); % to extract WM and subcortical nuclei
                f.mask           = fullfile(f.graph,[n,'.DiffusionSpace.nii']); % both gm and wm; final mask for building cerebrum graph. 
                
                f.mask_wm        = fullfile(f.graph,[n,'.DiffusionSpace.wm.nii']);
                
                f.mask_gm        = fullfile(f.graph,[n,'.DiffusionSpace.gm.nii']);
                f.mask_gm_ctx    = fullfile(f.graph,[n,'.DiffusionSpace.gm_ctx.nii']);
                f.mask_gm_subctx = fullfile(f.graph,[n,'.DiffusionSpace.gm_subctx.nii']);
                
                f.mask_gm_ctx_subctx_common = fullfile(f.graph,[n,'.DiffusionSpace.gm_ctx_subctx_common.nii']);
                
                f.mask_wm_gmborder          = fullfile(f.graph,[n,'.DiffusionSpace.wm_gmborder.nii']);
                f.mask_gm_wmborder          = fullfile(f.graph,[n,'.DiffusionSpace.gm_wmborder.nii']);
                
                f.mask  = fullfile(f.graph,[n,'.DiffusionSpace.nii']);
        end
        f.diffusion_data = fullfile(f.diffusion_save,'data.nii.gz');
        f.t1w_graphspace = fullfile(f.graph,'t1w.DiffusionSpace.nii');
        f.t2w_graphspace = fullfile(f.graph,'t2w.DiffusionSpace.nii');
        f.aparcaseg_graphspace = fullfile(f.graph,'aparc+aseg.DiffusionSpace.nii');
        f.aparcaseg = fullfile(f.T1w,'aparc+aseg.nii'); % not used in graph definition; just used to extract wm part within subcortical region used as graph mask, for debuging purposes
        f.wm_graphspace = fullfile(f.graph,'aparc+aseg_wm.DiffusionSpace.nii');
    case {'pslh','psrh','wslh','wsrh'}
        %f.source = fullfile(f.T1w,'ribbon.nii');
        %f.mask = fullfile(f.T1w,ID,'mri',[gtype(3:4),'.ribbon.nii']);
        f.ribbon_FS_1mm_mgz = fullfile(f.T1w,ID,'mri',[gtype(3:4),'.ribbon.mgz']); % 1mm cubic resolution in coregister with freesurfer orig.mgz using which synthetic data are defined [29 april 2021] 
        f.ribbon_FS_1mm_nii = fullfile(f.T1w,ID,'mri',[gtype(3:4),'.ribbon_1mm.nii']); % nii converted version of f.ribbon_FS_1mm_mgz
        % mask used as ref for reslicing preprocessed fMRI from MNI to G
        % space.
end

%-Surface files------------------------------------------------------------
% In gmxx design, pial files used to prune edges connnecting sulci banks.
% In wmxx design, white files used to define inter-cerebral mask.
%
d_surf = fullfile(hcpsave_root,ID,'T1w',ID,'surf'); %see NOTE 1.
[~,~] = mkdir(d_surf);

switch G.tissue
    case {'gmlh','gmrh','pslh','psrh','wslh','wsrh'}
        d = G.tissue(3:4);
        f.surface.pial = {fullfile(d_surf,[d,'.pial.surf.gii'])};
        f.surface.white = {fullfile(d_surf,[d,'.white.surf.gii'])};
        f.surfmesh.pial = fullfile(d_surf,[d,'.pial']);
        f.surfmesh.white = fullfile(d_surf,[d,'.white']);
        f.surfmesh.pial_ascii = fullfile(d_surf,[d,'.pial.asc']);
        f.surfmesh.white_ascii = fullfile(d_surf,[d,'.white.asc']);
    case {'wmlh','wmrh'}
        f.surface.white = {fullfile(d_surf,[gtype(3:4),'.white.surf.gii'])};
    case {'gm','cerebrum'}
        % for 'gm': pial and white surfaces are used for pruning. 
        % for 'cerebrum': pial and white surfaces used for pruning; white
        % surface is also used to determine edges between cortex and white.
        f.surface.pial = {...
            fullfile(d_surf,'lh.pial.surf.gii'),...
            fullfile(d_surf,'rh.pial.surf.gii')};
        
        f.surface.white = {...
            fullfile(d_surf,'lh.white.surf.gii'),...
            fullfile(d_surf,'rh.white.surf.gii')};
    case 'wm'
        f.surface.white = {...
            fullfile(d_surf,'lh.white.surf.gii'),...
            fullfile(d_surf,'rh.white.surf.gii')};
end

fn = fieldnames(f.surface);
for iFN=1:length(fn)
    d = f.surface.(fn{iFN});
    for iF=1:length(d)
        if ~exist(d{iF},'file')
            if exist(fullfile(hcp_root,ID,'preproc'),'dir') % ldcNAS format
                d1 = fullfile(hcp_root,ID,'preproc',ID,'T1w',ID,'surf');
            else
                d1 = fullfile(hcp_root,ID,'T1w',ID,'surf');
            end
            [~,d2,d3] = fileparts(d{iF});
            sts = copyfile(fullfile(d1,[d2,d3]),d{iF});
            if ~sts, error('[HB] problem copying surface file.'); end
        end
    end
end

%-Mapping files------------------------------------------------------------
f.xfms = fullfile(f.MNI,'xfms');
f.xfms_save = fullfile(f.MNI_save,'xfms');
[~,~] = mkdir(f.xfms_save);

% Displacement field for mapping from MNI to ACPC
f.disp_mni2acpc = fullfile(f.xfms,'standard2acpc_dc.nii');

% Displacement field for mapping from ACPC to MNI
f.disp_acpc2mni = fullfile(f.xfms,'acpc_dc2standard.nii');

%-Preprocessed fMRI volumes in MNI space-----------------------------------
% address to fMRI data that
% have been coregistered
% with graph.

%d = load('tasks9_LR.mat','tasks');
%taskSets{1} = d.tasks;

taskSets{1} = [
    {'tfMRI_EMOTION_LR'   }
    {'tfMRI_GAMBLING_LR'  }
    {'tfMRI_LANGUAGE_LR'  }
    {'tfMRI_MOTOR_LR'     }
    {'tfMRI_RELATIONAL_LR'}
    {'tfMRI_SOCIAL_LR'    }
    {'tfMRI_WM_LR'        }
    {'rfMRI_REST1_LR'     }
    {'rfMRI_REST2_LR'     }];

%d = load('tasks9_RL.mat','tasks');
%taskSets{2} = d.tasks;

taskSets{2} = [
    {'tfMRI_EMOTION_RL'   }
    {'tfMRI_GAMBLING_RL'  }
    {'tfMRI_LANGUAGE_RL'  }
    {'tfMRI_MOTOR_RL'     }
    {'tfMRI_RELATIONAL_RL'}
    {'tfMRI_SOCIAL_RL'    }
    {'tfMRI_WM_RL'        }
    {'rfMRI_REST1_RL'     }
    {'rfMRI_REST2_RL'     }];

for iS=1:2
    tasks = taskSets{iS}; 
    for t = 1:length(tasks)
        task = tasks{t};
        n = [task '.nii'];
        f.fmri_mni.(task)      = fullfile(f.MNI_results,task,n);               % on original volume
        f.fmri_mni_save.(task) = fullfile(f.MNI_results_save,task,n);          % on writable volume
    end
end

%-Resliced fMRI volumes that have 1-to-1 voxel macth with G.f.mask---------
% address to fMRI data that have been coregistered with graph.
for iS=1:2
    tasks = taskSets{iS};
    for t = 1:length(tasks)
        task = tasks{t};
        n = [task,G.rsTag,'.nii'];
        f.fmri_graph.(task) = fullfile(f.T1w_results_save,task,n);
    end
end

%-fMRI graph signals-------------------------------------------------------
% signals extracted from graph coregistered fMRI volumes.
for iS=1:2
    tasks = taskSets{iS};
    for t = 1:length(tasks)
        task = tasks{t};
        n = ['G.',G.type,'.signals_',task,'.mat'];
        f.signals.(task) = fullfile(f.graph,n);
    end
end

G.f = f;

G.fname = fullfile(f.graph,['G.',G.type,'.mat']);

if justGf
    sts = 'just G.f returned';
   return; 
end

if exist(G.fname,'file')
    d = load(G.fname);
    if isfield(d.G,'A') && isfield(d.G,'indices')
        G_ld = d.G;
        sts = 'loaded';
    else
        sts = 'new';
    end
else
    sts = 'new';
end

if strcmp(sts,'loaded')
    %-Update paths?
    % current machine may differ from that used when saving G.
    if isequal(G_ld.f,G.f)
        G = G_ld;
    else
        if isfield(G_ld.f,'odf')
            Gfodf = G_ld.f.odf; % G.f.odf is not defined in hb_get_G.m, but rather in hb_adjacencymatrix_diffusion.m. Therefore, it will result in a mistmatch if checked, so, we reomve it, check the rest of G.f, and we then replace it.
            G_ld.f = rmfield(G_ld.f,'odf');
            PutBackGfodf = true;
        else
            PutBackGfodf = false;
        end
        if isequal(G_ld.f,G.f)
            G = G_ld;
        else
            G_ld.f = G.f;
            G = G_ld;
        end
        if PutBackGfodf
            G.f.odf = Gfodf; % put back
        end
    end
end

if verbose
    fprintf('\n..Graph: %s',sts);
end
end
