% GRAph-based Spatial Decomposition (GRASD) of fMRI data using voxBG via 
% graph signal processing. 

clc
clear 
close all

which_graph = 1;
which_sess  = 1;
which_ID    = 1;
path_spm12  = []; 

%-Settings. 
%--------------------------------------------------------------------------
ID = '100307';

gtypes = {
    'gmlh.res2000.spaceT1w' % lh CHC voxBG
    'gmrh.res2000.spaceT1w' % rh ...
    'gmlh.res1250.spaceT1w' 
    'gmrh.res1250.spaceT1w' 
    };
gtype = gtypes{which_graph};

sessions = {
    'tfMRI_EMOTION_LR'
    'tfMRI_SOCIAL_LR'
    };
sess = sessions{which_sess};

d_scr     = fileparts(mfilename('fullpath'));
d_demos   = fileparts(d_scr);
d_main    = fileparts(d_demos);
d_src     = fullfile(d_main,'src','MATLAB'); 
d_hcpdata = fullfile(d_demos,'data','HCP','dataset');
addpath(genpath(d_src));

if ~isempty(path_spm12)
    addpath(paths.spm12);
end

%-Load voxBG.
%--------------------------------------------------------------------------
d_graph = fullfile(d_hcpdata,ID,'T1w','graph',strrep(gtype,'.','_'));
n_graph = sprintf('G.%s.mat',gtype);
f_load = fullfile(d_graph,n_graph);
d = load(f_load);
G = d.G;

%-Load SOSKS.
%--------------------------------------------------------------------------

%-Load graph signals.
%--------------------------------------------------------------------------

%-Do GRASD. 
%--------------------------------------------------------------------------

disp('done.');
