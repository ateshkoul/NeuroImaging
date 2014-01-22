function dcm = tbx_cfg_dcm
% Configuration file for DCM toolbox
% Written by Atesh Koul, National Brain Research Center, India

%%

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


name        = cfg_entry;
name.tag     = 'name';
name.name    = 'Name of DCM';
name.help    = {'Name of the DCM file to be appended. It will appear as DCM_name_'};
name.val     = {''};
name.strtype = 's';
name.num     = [1 Inf];


savedcm         = cfg_files;
savedcm.tag     = 'savedcm';
savedcm.name    = 'Directory to save models';
savedcm.help    = {'Select directory to save created models'};
savedcm.filter  = 'dir';
savedcm.ufilter ='.*';
savedcm.num     = [1 Inf];

% ---------------------------------------------------------------------
% VOI
% ---------------------------------------------------------------------
voi        = cfg_files;
voi.tag     = 'voi';
voi.name    = 'VOI';
voi.help    = {'Select the Volume of interest .mat file'};
voi.filter  = 'mat';
voi.ufilter = '^VOI.*\.mat$';
voi.num     = [1 1];

%--------------------------------------------------
% Non-linear effects
%--------------------------------------------------
%d         = cfg_entry;
%d.tag     = 'd';
%d.name    = 'Non-linear effects';
%d.help    = {'Non-linear effects'};
%d.val     = {[0 0 0;0 0 0;0 0 0]};
%d.strtype = 'e';
%d.num     = [3 3];


vois         = cfg_branch;
vois.tag     = 'vois';
vois.name    = 'VOIs';
vois.val     = {voi};
vois.help    = {'Select the VOIs'};

repvoi         =  cfg_repeat;
repvoi.tag     = 'repvoi';
repvoi.name    = 'VOIs';
repvoi.help    = {'Select the VOIs'};
repvoi.values  = {vois};
repvoi.num     = [1 8];

%--------------------------------------------------
% Matrix selection
%--------------------------------------------------

matrix         = cfg_files;
matrix.tag     = 'matrix';
matrix.name    = 'Select matrix';
matrix.help    = {'Select the matrix file that contains the design specification.'
                   'It should contain values of intrinsic, modulatory and input matrices in 1 DCM struct as DCM.a, DCM.b, DCM.c, DCM.d (a b c and d parameters)'};
matrix.filter  = 'mat';
matrix.ufilter = '.*';
matrix.num     = [1 1];

%% Model Construction
%
% Select name of destination folder
%
dest         = cfg_files;
dest.tag     = 'dest';
dest.name    = 'Directory';
dest.help    = {'Select a destination directory where models will be saved.'};
dest.filter = 'dir';
dest.ufilter = '.*';
dest.num     = [1 1];
%
% Number of ROIs
%
%nROIs         = cfg_entry;
%nROIs.tag     = 'nROIs';
%nROIs.name    = 'no. of ROIs';
%nROIs.help    = {'No. of conditions'};
%nROIs.val     = {};
%nROIs.strtype = 'n';
%nROIs.num     = [1 1];


save         = cfg_menu;
save.tag     = 'save';
save.name    = 'Save generated Models';
save.help    = {'Choose whether to save generated models or not'
                 'Only model parameters a b c will be saved'};
save.labels  = {'Yes',...        
               'No'}';
save.values  = { 1 0};
save.val     = {1};

intmatrix         = cfg_entry;
intmatrix.tag     = 'intmatrix';
intmatrix.name    = 'Select general intrinsic matrix';
intmatrix.help    = {'Select the generalised intrinsic matrix that contains the design specification. This has to be in the form nxn where n is the number of VOIs.' 
                      'The sequence of VOIs has to be in the order starting from input to final region'};
intmatrix.val     = {};
intmatrix.strtype = 'e';
intmatrix.num     = [Inf Inf];


modmatrix         = cfg_entry;
modmatrix.tag     = 'modmatrix';
modmatrix.name    = 'Select general modulatory matrix';
modmatrix.help    = {'Select the generalised modulatory matrix that contains the design specification.'
                    'The matrix has to correspond to the generalised intrinsic matrix'};
modmatrix.val     = {};
modmatrix.strtype = 'e';
modmatrix.num     = [Inf Inf];

inpmatrix         = cfg_entry;
inpmatrix.tag     = 'inpmatrix';
inpmatrix.name    = 'Select general input matrix';
inpmatrix.help    = {'Select the generalised input matrix that contains the design specification.'
                     'The matrix has to be in the form nxc where n is number of VOIs and c is number of conditions.'
                     'The conditions have to coorespond to that in the selected SPM mat file'};
inpmatrix.val     = {};
inpmatrix.strtype = 'e';
inpmatrix.num     = [Inf Inf];

modcreate        = cfg_exbranch;
modcreate.tag     = 'modcreate';
modcreate.name    = 'Create Models';
modcreate.val     = {dest save intmatrix modmatrix inpmatrix};
modcreate.help    = {'This options creats all combinations of possible models given a generalised matrix of possible connections'};
%%

model         = cfg_choice;
model.tag     = 'model';
model.name    = 'Model Specification';
model.val     = {};
model.help    = {'Select whether to create models or select a model matrix'};
model.values  = {modcreate matrix};

cond         = cfg_entry;
cond.tag     = 'cond';
cond.name    = 'No. of conditions';
cond.help    = {'No. of conditions'};
cond.val     = {2};
cond.strtype = 'n';
cond.num     = [1 1];

%% Model Options
delays         = cfg_entry;
delays.tag     = 'delays';
delays.name    = 'slice timings';
delays.help    = {'Specify slice timings 1 for default value'};
delays.strtype = 'n';
delays.val     = {1};
delays.num     = [1  1];

% Echo time

te         = cfg_entry;
te.tag     = 'te';
te.name    = 'Echo time';
te.help    = {'Specify Echo time'};
te.strtype = 'n';
te.val     = {3};
te.num     = [1  1];

% bilinear or non-linear
lin         = cfg_menu;
lin.tag     = 'lin';
lin.name    = 'Linearity';
lin.help    = {'Specify Model type'};
lin.labels  =  {'Bilinear' 'Non-linear'};
lin.values  = {0 1};
lin.val     = {0};


% States per region
state         = cfg_menu;
state.tag     = 'state';
state.name    = 'States per region';
state.help    = {'Specify states per region'};
state.labels  =  {'one' 'two'};
state.values  = {0 1};
state.val     = {0};

% Stochastic 
sto        = cfg_menu;
sto.tag     = 'sto';
sto.name    = 'Stochastic';
sto.help    = {'Specify stochastic'};
sto.labels  =  {'No' 'Yes'};
sto.values  = {0 1};
sto.val     = {0};

% Specify center
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
options.help    = {'Select DCM options'};

%%

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
specify.val     = {spmmat name savedcm repvoi model cond options};
specify.help    = {'Batch specify and create DCM models'};
specify.prog    = @batch_dcm;


%% Batch DCM Estimate
%
% Estimating DCM models in batch mode
%

% ---------------------------------------------------------------------
% Select DCM files
% ---------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'DCM files';
data.help    = {'Select DCM files to estimate'};
data.filter = 'mat';
data.ufilter = '.mat';
data.num     = [1 Inf];

% ---------------------------------------------------------------------
% Estimate
% ---------------------------------------------------------------------
est        = cfg_exbranch;
est.tag     = 'est';
est.name    = 'Estimate';
est.val     = {data };
est.help    = {'Batch Estimate DCM models'};
est.prog = @batch_dcm_estimate;
% the program batch_dcm_estimate has to be in the path for this to work

%%
dcm         = cfg_choice;
dcm.tag     = 'dcm';
dcm.name    = 'DCM';
dcm.help    = {'Generate and Estimate DCM models in batch'};
dcm.values     = {specify est};
