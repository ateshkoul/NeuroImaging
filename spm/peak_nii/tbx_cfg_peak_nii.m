function peaknii = tbx_cfg_peak_nii
% Configuration file for toolbox 'AAL'
% Written by Atesh Koul, National Brain Research Center, India

% ---------------------------------------------------------------------
% spmmat Select SPM.mat
% ---------------------------------------------------------------------
image         = cfg_files;
image.tag     = 'image';
image.name    = 'Select images';
image.help    = {'Select the images for co-ordinates.'};
image.filter  = 'image';
image.ufilter = '^*.img';
image.num     = [1 1];
% ---------------------------------------------------------------------
% titlestr 
% ---------------------------------------------------------------------
uid         = cfg_entry;
uid.tag     = 'uid';
uid.name    = 'UID';
uid.help    = {'unique ID for the analysis'};
uid.val     = {''};
uid.strtype = 's';
uid.num     = [0 Inf];


analysisinfo         = cfg_branch;
analysisinfo.tag     = 'analysisinfo';
analysisinfo.name    = 'Analysis';
analysisinfo.val     = {image uid};
analysisinfo.help    = {''};

% ---------------------------------------------------------------------
% generic Contrasts
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Analysis info';
generic.help    = {''};
generic.values  = {analysisinfo};
generic.num     = [1 Inf];

% ---------------------------------------------------------------------
% contrasts Contrast(s)
% ---------------------------------------------------------------------
out         = cfg_entry;
out.tag     = 'out';
out.name    = 'Output Prefix';
out.help    = {'Output Prefix'};
out.val     = {''};
out.strtype = 's';
out.num     = [0 Inf];
% ---------------------------------------------------------------------
% threshdesc Threshold type
% ---------------------------------------------------------------------
sign         = cfg_menu;
sign.tag     = 'sign';
sign.name    = 'Sign';
sign.help    = {''};
sign.labels  = {'pos' 'neg'};
sign.values  = {'pos' 'neg'};
sign.val     = {'pos'};
% ---------------------------------------------------------------------
% thresh Threshold
% ---------------------------------------------------------------------
thresh         = cfg_entry;
thresh.tag     = 'thresh';
thresh.name    = 'Threshold';
thresh.help    = {''};
thresh.strtype = 'e';
thresh.num     = [0 1];
thresh.val     = {0.05};
% ---------------------------------------------------------------------
% extent Extent (voxels)
% ---------------------------------------------------------------------
type         = cfg_menu;
type.tag     = 'type';
type.name    = 'type of contrast';
type.help    = {'either t or f z or none'};
type.labels  =  {'t' 'f' 'z' 'none'};
type.values  = {'t' 'f' 'z' 'none'};
type.val     = ['t' 'f' 'z' 'none'];
% ---------------------------------------------------------------------
% contrasts Contrast(s)
% ---------------------------------------------------------------------
voxlimit         = cfg_entry;
voxlimit.tag     = 'voxlimit';
voxlimit.name    = 'voxlimit';
voxlimit.help    = {'No. of Voxels to show default = 1000.'};
voxlimit.strtype = 'e';
voxlimit.num     = [1 1];
voxlimit.val     = {1000};
% ---------------------------------------------------------------------
% thresh Mask threshold
% ---------------------------------------------------------------------
separation         = cfg_entry;
separation.tag     = 'separation';
separation.name    = 'separation';
separation.help    = {'distance two peaks must be separated'};
separation.strtype = 'e';
separation.num     = [1 1];
separation.val     = {8};
% ---------------------------------------------------------------------
% mtype Nature of mask
% ---------------------------------------------------------------------
SPM         = cfg_menu;
SPM.tag     = 'SPM';
SPM.name    = 'Collapsing across peaks';
SPM.help    = {''};
SPM.labels  = {'Yes' 'No'};
SPM.values  = {0 1};

conn         = cfg_menu;
conn.tag     = 'conn';
conn.name    = 'connectivity radius';
conn.help    = {''};
conn.labels  = {'6' '12' '26'};
conn.values  = {6 12 26};


cluster         = cfg_entry;
cluster.tag     = 'cluster';
cluster.name    = 'minimum cluster size';
cluster.help    = {'minimum cluster size required to keep a cluster in the results'};
cluster.strtype = 'e';
cluster.num     = [0 1];
cluster.val     = {8};

% ---------------------------------------------------------------------
% mask Mask definition
% ---------------------------------------------------------------------
df1         = cfg_entry;
df1.tag     = 'df1';
df1.name    = 'degrees of freedom';
df1.help    = {'degrees of freedom for a T-test'};
df1.strtype = 'e';
df1.num     = [0 1];
df1.val     = {''};



df2         = cfg_entry;
df2.tag     = 'df2';
df2.name    = 'denominator degrees of freedom';
df2.help    = {'denominator degrees of freedom for an F-test'};
df2.strtype = 'e';
df2.num     = [0 1];
df2.val     = {''};

mask         = cfg_files;
mask.tag     = 'mask';
mask.name    = 'Mask file';
mask.help    = {'Select the mask image.'};
mask.filter  = 'image';
mask.ufilter = '^*.img';
mask.num     = [0 1];
mask.val     = {''};

exact         = cfg_menu;
exact.tag     = 'exact';
exact.name    = 'Collapsing across peaks';
exact.help    = {''};
exact.labels  = {'0' '1'};
exact.values  = {0 1};

sphere         = cfg_entry;
sphere.tag     = 'sphere';
sphere.name    = 'radius of sphere';
sphere.help    = {'denominator degrees of freedom for an F-test'};
sphere.strtype = 'e';
sphere.num     = [0 1];


clustersphere         = cfg_entry;
clustersphere.tag     = 'clustersphere';
clustersphere.name    = 'radius of clustersphere';
clustersphere.help    = {'denominator degrees of freedom for an F-test'};
clustersphere.strtype = 'e';
clustersphere.num     = [0 1];


param         = cfg_branch;
param.tag     = 'param';
param.name    = 'Parameters';
param.val     = {out sign thresh type voxlimit separation SPM conn  cluster df1 df2 mask exact sphere clustersphere};
param.help    = {''};

% generic Masking
% ---------------------------------------------------------------------


peaknii         = cfg_exbranch;
peaknii.tag     = 'peaknii';
peaknii.name    = 'peak_nii';
peaknii.val     = {generic param};
peaknii.help    = {'AUTOMATED ANATOMICAL LABELING toolbox:Reference Automated Anatomical Labeling of Activations in SPM Using a Macroscopic Anatomical Parcellation',... 
                'of the MNI MRI Single-Subject Brain.',...
                'N. Tzourio-Mazoyer, B. Landeau, D. Papathanassiou, F. Crivello, O. Etard, N. Delcroix, B. Mazoyer, and M. Joliot. NeuroImage 2002.'};
peaknii.prog    = @spm_local_peak_nii;

%======================================================================
function spm_local_peak_nii(job)

if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','peak_nii')); end
batch_peak_nii(job);
