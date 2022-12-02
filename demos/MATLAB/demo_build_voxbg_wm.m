clc
clear 
close all

%path_spm12  = []; 
path_spm12 = '/Users/hamid/Dropbox/my/my_toolboxes/external/spm12';

%-Settings. 
%--------------------------------------------------------------------------
ID = '100307';

gtype = 'wm.res1250.spaceDiffusion.neighb5.pow1.fibs5.rempar1.sh8.th90.beta50.mag0.soliangshAverageuniform.onlywm1';

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