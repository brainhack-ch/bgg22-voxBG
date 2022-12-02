clc
clear 
close all

which_graph = 1;

%path_spm12  = []; 
path_spm12 = '/Users/hamid/Dropbox/my/my_toolboxes/external/spm12';


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
d_main    = fileparts(d_demos);
d_src     = fullfile(d_main,'src','MATLAB'); 
d_hcpdata = fullfile(d_demos,'data','HCP','dataset');
addpath(genpath(d_src));

if ~isempty(path_spm12)
    addpath(path_spm12);
end

opts = struct('hcp_root',d_hcpdata);

%-Build voxBG.
%--------------------------------------------------------------------------
voxbg_build(ID,gtype,opts);

fprintf('\ndone.');