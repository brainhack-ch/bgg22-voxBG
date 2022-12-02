clc
clear 
close all

which_graph = 1;
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

d_scr     = fileparts(mfilename('fullpath'));
d_demos   = fileparts(d_scr);
d_hcpdata = fullfile(d_demos,'data','HCP','dataset');
d_hbrepos = '/Users/hamid/Dropbox/github';
getfiles(d_scr,d_hbrepos);
addpath(genpath(fullfile(d_scr,'utils')));
addpath(genpath(fullfile(d_scr,'resources')));

if ~isempty(path_spm12)
    addpath(path_spm12);
end

%-Load voxBG.
%--------------------------------------------------------------------------
d_graph = fullfile(d_hcpdata,ID,'T1w','graph',strrep(gtype,'.','_'));
n_graph = sprintf('G.%s.mat',gtype);
f_load = fullfile(d_graph,n_graph);
d = load(f_load);
G = d.G;

%-Visualize voxBG.
%--------------------------------------------------------------------------

disp('done.');
