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
roi.ufilter = '.mat$';
roi.num     = [1 Inf];

% % ---------------------------------------------------------------------
% % contrasts Contrast(s)
% % ---------------------------------------------------------------------
% voi         = cfg_repeat;
% voi.tag     = 'voi';
% voi.name    = 'VOI';
% voi.help    = {
%                    'The design matrix defines the experimental design and the nature of hypothesis testing to be implemented.  The design matrix has one row for each scan and one column for each effect or explanatory variable. (e.g. regressor or stimulus function).  '
%                    ''
%                    'This allows you to build design matrices with separable session-specific partitions.  Each partition may be the same (in which case it is only necessary to specify it once) or different.  Responses can be either event- or epoch related, where the latter model involves prolonged and possibly time-varying responses to state-related changes in experimental cond.  Event-related response are modelled in terms of responses to instantaneous events.  Mathematically they are both modelled by convolving a series of delta (stick) or box-car functions, encoding the input or stimulus function. with a set of hemodynamic basis functions.'
% }';
% voi.values  = {roi };
% voi.num     = [1 Inf];

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

d         = cfg_entry;
d.tag     = 'd';
d.name    = 'Non-linear effects';
d.help    = {'Non-linear effects'};
d.val     = {[0 0 0;0 0 0;0 0 0]};
d.strtype = 'e';
d.num     = [3 3];


effects         = cfg_branch;
effects.tag     = 'effects';
effects.name    = 'Timing parameters';
effects.val     = {b c d};
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
te.val     = {30};
te.num     = [1  1];

% 
lin         = cfg_menu;
lin.tag     = 'lin';
lin.name    = 'Linearity';
lin.help    = {'Specify Model type'};
lin.labels  =  {'Bilinear' 'Non-linear'};
lin.values  = {0 1};
lin.val     = {1};


% 
state         = cfg_menu;
state.tag     = 'state';
state.name    = 'States per region';
state.help    = {'Specify states per region'};
state.labels  =  {'one' 'two'};
state.values  = {0 1};
state.val     = {1};

sto        = cfg_menu;
sto.tag     = 'sto';
sto.name    = 'Stocastic';
sto.help    = {'Specify stocastic'};
sto.labels  =  {'No' 'Yes'};
sto.values  = {0 1};
sto.val     = {1};

cen         = cfg_menu;
cen.tag     = 'cen';
cen.name    = 'center';
cen.help    = {'Specify center'};
cen.labels  =  {'No' 'Yes'};
cen.values  = {0 1};
cen.val     = {1};

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
specify.name    = 'DCM';
specify.val     = {spmmat name roi a generic delays te lin state sto cen};
specify.help    = {'help 2002.'};
specify.prog    = @spm_local_dcm;
specify.vout    = @vout_specify;

data         = cfg_files;
data.tag     = 'data';
data.name    = 'Session';
data.help    = {'Select scans for this session. In the coregistration step, the sessions are first realigned to each other, by aligning the first scan from each session to the first scan of the first session.  Then the images within each session are aligned to the first image of the session. The parameter estimation is performed this way because it is assumed (rightly or not) that there may be systematic differences in the images between sessions.'};
data.filter = 'mat';
data.ufilter = '.mat';
data.num     = [1 1];

est         = cfg_exbranch;
est.tag     = 'est';
est.name    = 'Estimate';
est.val     = {data };
est.help    = {
                    'This routine realigns a time-series of images acquired from the same subject using a least squares approach and a 6 parameter (rigid body) spatial transformation/* \cite{friston95a}*/.  The first image in the list specified by the user is used as a reference to which all subsequent scans are realigned. The reference scan does not have to be the first chronologically and it may be wise to chose a "representative scan" in this role.'
                    ''
                    'The aim is primarily to remove movement artefact in fMRI and PET time-series (or more generally longitudinal studies) /* \cite{ashburner97bir}*/. The headers are modified for each of the input images, such that. they reflect the relative orientations of the data. The details of the transformation are displayed in the results window as plots of translation and rotation. A set of realignment parameters are saved for each session, named rp_*.txt. After realignment, the images are resliced such that they match the first image selected voxel-for-voxel. The resliced images are named the same as the originals, except that they are prefixed by ''r''.'
}';
est.prog = @spm_dcm_estimate;

dcm         = cfg_choice;
dcm.tag     = 'dcm';
dcm.name    = 'DCM';
dcm.help    = {'Within-subject registration of image time series.'};
dcm.values  = {specify est };
%
%

%======================================================================
function spm_local_dcm(job)

if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','DCM')); end

batch_dcm(job);


function dep = vout_specify(job)
    cdep(1)            = cfg_dep;
    cdep(1).sname      = sprintf('Realignment Param File (Sess %d)', k);
    cdep(1).src_output = substruct('.','sess', '()',{k}, '.','rpfile');
    cdep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
    dep = cdep;
end;
