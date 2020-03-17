addpath('functions/ht_pnp_function');
addpath('functions/at_netvlad_function');
addpath('functions/utils');
addpath('functions/wustl_function');
addpath('functions/relja_matlab');
addpath('functions/relja_matlab/matconvnet/');
addpath('functions/netvlad/');
addpath('functions/inpaint_nans');

env = environment();
if strcmp(env, 'ciirc')
    addpath('/home/lucivpav/NetVLAD');
elseif strcmp(env, 'cmp')
    addpath('/datagrid/personal/lucivpav/NetVLAD');
elseif strcmp(env, 'laptop')
    addpath('/Volumes/GoogleDrive/Můj disk/ARTwin/InLocCIIRC_dataset/NetVLAD');
end

addpath('functions/yael_matlab_linux64_v438');
run('functions/vlfeat/toolbox/vl_setup.m');
run('functions/matconvnet/matlab/vl_setupnn.m');
