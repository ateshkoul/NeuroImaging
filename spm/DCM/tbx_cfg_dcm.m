function dcm = tbx_cfg_dcm
% Configuration file for toolbox 'AAL'
% Written by Atesh Koul, National Brain Research Center, India

% ---------------------------------------------------------------------
% spmmat Select SPM.mat
% ---------------------------------------------------------------------
spmmat         = cfg_files;
spmmat.tag     = 'spmmat';
spmmat.name    = 'Select SPM.mat';
spmmat.help    = {'Select the SPM.mat file that contains the design specification.'};
spmmat.filter  = 'mat';
spmmat.ufilter = '^SPM\.mat$';
spmmat.num     = [1 1];
% ---------------------------------------------------------------------
% ROI
% ---------------------------------------------------------------------

name        = cfg_entry;
name.tag     = 'name';
name.name    = 'Name of DCM';
name.help    = {'Name of the DCM file'};
name.val     = {''};
name.strtype = 's';
name.num     = [1 Inf];


roi        = cfg_files;
roi.tag     = 'roi';
roi.name    = 'ROIs';
roi.help    = {'Select the Region of interest .mat file'};
roi.filter  = 'mat';
roi.ufilter = '^VOI.*\.mat$';
roi.num     = [1 1];

d         = cfg_entry;
d.tag     = 'd';
d.name    = 'Non-linear effects';
d.help    = {'Non-linear effects'};
d.val     = {[0 0 0;0 0 0;0 0 0]};
d.strtype = 'e';
d.num     = [3 3];


rois         = cfg_branch;
rois.tag     = 'rois';
rois.name    = 'ROIs';
rois.val     = {roi d};
rois.help    = {
                  'Specify various timing parameters needed to construct the design matrix. This includes the units of the design specification and the interscan interval.'
                  ''
                  'Also, with longs TRs you may want to shift the regressors so that they are aligned to a particular slice.  This is effected by changing the microtime resolution and onset. '
}';

reproi         =  cfg_repeat;
reproi.tag     = 'reproi';
reproi.name    = 'ROIs';
reproi.help    = {'.'};
reproi.values  = {rois};
reproi.num     = [1 8];

% ---------------------------------------------------------------------
% Directory for the ROI
% ---------------------------------------------------------------------

a         = cfg_entry;
a.tag     = 'a';
a.name    = 'Intrinsic connectivity';
a.help    = {'Intrinsic connectivity'};
a.val     = {''};
a.strtype = 'e';
a.num     = [Inf Inf];

b        = cfg_entry;
b.tag     = 'b';
b.name    = 'Effect of condition';
b.help    = {'Effect of condition'};
b.val     = {''};
b.strtype = 'e';
b.num     = [Inf Inf];

c         = cfg_entry;
c.tag     = 'c';
c.name    = 'Direct Effect';
c.help    = {'Direct Effect'};
c.val     = {''};
c.strtype = 'e';
c.num     = [1 Inf];


effects         = cfg_branch;
effects.tag     = 'effects';
effects.name    = 'Conditions';
effects.val     = {b c};
effects.help    = {
                  'Specify various timing parameters needed to construct the design matrix. This includes the units of the design specification and the interscan interval.'
                  ''
                  'Also, with longs TRs you may want to shift the regressors so that they are aligned to a particular slice.  This is effected by changing the microtime resolution and onset. '
}';

generic         =  cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Conditions';
generic.help    = {'.'};
generic.values  = {effects };
generic.num     = [1 Inf];

delays         = cfg_entry;
delays.tag     = 'delays';
delays.name    = 'slice timings';
delays.help    = {'Specify slice timings 1 for default value'};
delays.strtype = 'n';
delays.val     = {1};
delays.num     = [1  1];
% 

te         = cfg_entry;
te.tag     = 'te';
te.name    = 'Echo time';
te.help    = {'Specify Echo time'};
te.strtype = 'n';
te.val     = {3};
te.num     = [1  1];

% 
lin         = cfg_menu;
lin.tag     = 'lin';
lin.name    = 'Linearity';
lin.help    = {'Specify Model type'};
lin.labels  =  {'Bilinear' 'Non-linear'};
lin.values  = {0 1};
lin.val     = {0};


% 
state         = cfg_menu;
state.tag     = 'state';
state.name    = 'States per region';
state.help    = {'Specify states per region'};
state.labels  =  {'one' 'two'};
state.values  = {0 1};
state.val     = {0};

sto        = cfg_menu;
sto.tag     = 'sto';
sto.name    = 'Stocastic';
sto.help    = {'Specify stocastic'};
sto.labels  =  {'No' 'Yes'};
sto.values  = {0 1};
sto.val     = {0};

cen         = cfg_menu;
cen.tag     = 'cen';
cen.name    = 'center';
cen.help    = {'Specify center'};
cen.labels  =  {'No' 'Yes'};
cen.values  = {0 1};
cen.val     = {1};

options        = cfg_branch;
options.tag     = 'options';
options.name    = 'Model options';
options.val     = {delays te lin state sto cen};
options.help    = {
                  'Specify various timing parameters needed to construct the design matrix. This includes the units of the design specification and the interscan interval.'
                  ''
                  'Also, with longs TRs you may want to shift the regressors so that they are aligned to a particular slice.  This is effected by changing the microtime resolution and onset. '
}';

% ---------------------------------------------------------------------
% titlestr
% ---------------------------------------------------------------------
titlestr         = cfg_entry;
titlestr.tag     = 'titlestr';
titlestr.name    = 'DCM specify';
titlestr.help    = {'DCM in batch mode'};
titlestr.val     = {''};
titlestr.strtype = 's';
titlestr.num     = [0 Inf];

specify         = cfg_exbranch;
specify.tag     = 'specify';
specify.name    = 'DCM specify';
specify.val     = {spmmat name reproi a generic options};
specify.help    = {' 2002.'};
specify.prog    = @batch_dcm;
% specify.vout    = @vout_specify;
% have to make the dependency thing to work

% ---------------------------------------------------------------------
% data Session
% ---------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'Session';
data.help    = {'Select scans for this session. In the coregistration step, the sessions are first realigned to each other, by aligning the first scan from each session to the first scan of the first session.  Then the images within each session are aligned to the first image of the session. The parameter estimation is performed this way because it is assumed (rightly or not) that there may be systematic differences in the images between sessions.'};
data.filter = 'mat';
data.ufilter = '.mat';
data.num     = [1 1];

% ---------------------------------------------------------------------
% estwrite Realign: Estimate & Reslice
% ---------------------------------------------------------------------
est        = cfg_exbranch;
est.tag     = 'est';
est.name    = 'Estimate';
est.val     = {data };
est.help    = {
                    'This routine realigns a time-series of images acquired from the same subject using a least squares approach and a 6 parameter (rigid body) spatial transformation/* \cite{friston95a}*/.  The first image in the list specified by the user is used as a reference to which all subsequent scans are realigned. The reference scan does not have to be the first chronologically and it may be wise to chose a "representative scan" in this role.'
                    ''
                    'The aim is primarily to remove movement artefact in fMRI and PET time-series (or more generally longitudinal studies) /* \cite{ashburner97bir}*/. The headers are modified for each of the input images, such that. they reflect the relative orientations of the data. The details of the transformation are displayed in the results window as plots of translation and rotation. A set of realignment parameters are saved for each session, named rp_*.txt. After realignment, the images are resliced such that they match the first image selected voxel-for-voxel. The resliced images are named the same as the originals, except that they are prefixed by ''r''.'
}';
est.prog = @batch_dcm_estimate;
% the program batch_dcm_estimate has to be in the path for this to work

dcm         = cfg_choice;
dcm.tag     = 'dcm';
dcm.name    = 'DCM';
dcm.help    = {' 2002.'};
dcm.values     = {specify est};


% %======================================================================
% function dep = vout_specify(job)
%     dep(1)            = cfg_dep;
%     dep(1).sname      = sprintf('new');
%     dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
%     


