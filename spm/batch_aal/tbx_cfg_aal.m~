function aal = tbx_cfg_aal
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
% titlestr 
% ---------------------------------------------------------------------
titlestr         = cfg_entry;
titlestr.tag     = 'titlestr';
titlestr.name    = 'Results Title';
titlestr.help    = {'Heading on aal -determined automatically if left empty'};
titlestr.val     = {''};
titlestr.strtype = 's';
titlestr.num     = [0 Inf];
% ---------------------------------------------------------------------
% contrasts Contrast(s)
% ---------------------------------------------------------------------
contrasts         = cfg_entry;
contrasts.tag     = 'contrasts';
contrasts.name    = 'Contrast(s)';
contrasts.help    = {
                     'Index of contrast(s). If more than one number is entered, analyse a conjunction hypothesis.'
                     ''
                     'If only one number is entered, and this number is "Inf", then results are printed for all contrasts found in the SPM.mat file.'
}';
contrasts.strtype = 'e';
contrasts.num     = [1 Inf];
% ---------------------------------------------------------------------
% threshdesc Threshold type
% ---------------------------------------------------------------------
threshdesc         = cfg_menu;
threshdesc.tag     = 'threshdesc';
threshdesc.name    = 'Threshold type';
threshdesc.help    = {''};
threshdesc.labels  = {'FWE' 'FDR' 'none'};
threshdesc.values  = {'FWE' 'FDR' 'none'};
threshdesc.val     = {'FWE'};
% ---------------------------------------------------------------------
% thresh Threshold
% ---------------------------------------------------------------------
thresh         = cfg_entry;
thresh.tag     = 'thresh';
thresh.name    = 'Threshold';
thresh.help    = {''};
thresh.strtype = 'e';
thresh.num     = [1 1];
thresh.val     = {0.05};
% ---------------------------------------------------------------------
% extent Extent (voxels)
% ---------------------------------------------------------------------
extent         = cfg_entry;
extent.tag     = 'extent';
extent.name    = 'Extent (voxels)';
extent.help    = {''};
extent.strtype = 'e';
extent.num     = [1 1];
extent.val     = {0};
% ---------------------------------------------------------------------
% contrasts Contrast(s)
% ---------------------------------------------------------------------
contrasts1         = cfg_entry;
contrasts1.tag     = 'contrasts';
contrasts1.name    = 'Contrast(s)';
contrasts1.help    = {'Index of contrast(s) for masking - leave empty for no masking.'};
contrasts1.strtype = 'e';
contrasts1.num     = [1 Inf];
% ---------------------------------------------------------------------
% thresh Mask threshold
% ---------------------------------------------------------------------
thresh1         = cfg_entry;
thresh1.tag     = 'thresh';
thresh1.name    = 'Mask threshold';
thresh1.help    = {''};
thresh1.strtype = 'e';
thresh1.num     = [1 1];
thresh1.val     = {0.05};
% ---------------------------------------------------------------------
% mtype Nature of mask
% ---------------------------------------------------------------------
mtype         = cfg_menu;
mtype.tag     = 'mtype';
mtype.name    = 'Nature of mask';
mtype.help    = {''};
mtype.labels  = {'Inclusive' 'Exclusive'};
mtype.values  = {0 1};
% ---------------------------------------------------------------------
% mask Mask definition
% ---------------------------------------------------------------------
mask         = cfg_branch;
mask.tag     = 'mask';
mask.name    = 'Mask definition';
mask.val     = {contrasts1 thresh1 mtype};
mask.help    = {''};
% ---------------------------------------------------------------------
% generic Masking
% ---------------------------------------------------------------------
generic1         = cfg_repeat;
generic1.tag     = 'generic';
generic1.name    = 'Masking';
generic1.help    = {''};
generic1.values  = {mask};
generic1.num     = [0 1];
% ---------------------------------------------------------------------
% conspec Contrast query
% ---------------------------------------------------------------------
conspec         = cfg_branch;
conspec.tag     = 'conspec';
conspec.name    = 'Contrast query';
conspec.val     = {titlestr contrasts threshdesc thresh extent generic1};
conspec.help    = {''};
% ---------------------------------------------------------------------
% generic Contrasts
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Contrasts';
generic.help    = {''};
generic.values  = {conspec};
generic.num     = [1 Inf];
% ---------------------------------------------------------------------
% units Units
% ---------------------------------------------------------------------
units         = cfg_menu;
units.tag     = 'units';
units.name    = 'Data type';
units.help    = {['Data type. This option is only meaningful for M/EEG '...
    'data. Keep the default ''Volumetric'' for any other kind of data.']};
units.labels  = {'Volumetric (2D/3D)',...
                 'Scalp-Time',...
                 'Scalp-Frequency',...
                 'Time-Frequency',...
                 'Frequency-Frequency'}';
units.values  = { 1 2 3 4 5 };
units.val     = { 1 };
% ---------------------------------------------------------------------
% print Print results
% ---------------------------------------------------------------------
pr         = cfg_menu;
pr.tag     = 'pr';
pr.name    = 'Print results';
pr.help    = {''};
pr.labels  = {'Yes',...               
               'No'}';
pr.values  = { 1 0 };
pr.val     = { 1 };

% ---------------------------------------------------------------------
% Select a labelled Atlas
% ---------------------------------------------------------------------

atl         = cfg_files;
atl.tag     = 'atl';
atl.name    = 'Atlas';
atl.help    = {'Select a labelled atlas.'};
atl.filter  = 'image';
atl.ufilter = '.*';
atl.num     = [1 1];

% ---------------------------------------------------------------------
% titlestr
% ---------------------------------------------------------------------
titlestr         = cfg_entry;
titlestr.tag     = 'titlestr';
titlestr.name    = 'AAL';
titlestr.help    = {'Heading on aal - determined automatically if left empty'};
titlestr.val     = {''};
titlestr.strtype = 's';
titlestr.num     = [0 Inf];

% ---------------------------------------------------------------------
% Select labelling option
% ---------------------------------------------------------------------
option         = cfg_menu;
option.tag     = 'option';
option.name    = 'Select option';
option.help    = {['Select Labelling.']};
option.labels  = {'Local Maxima Labeling',...
                 'Extended Local Maxima Labeling',...
                 'Cluster Labeling'}';
option.values  = { 1 2 3 };
option.val     = { 1 };

% ---------------------------------------------------------------------
% Select radius for extended local maxima labelling
% ---------------------------------------------------------------------
extn         = cfg_entry;
extn.tag     = 'extn';
extn.name    = 'Radius(Extended local maxima)';
extn.help    = {'local maximum radius (mm)for only Extended local maxima'};
extn.strtype = 'e';
extn.num     = [1 1];
extn.val     = {0};

aal         = cfg_exbranch;
aal.tag     = 'aal';
aal.name    = 'AAL';
aal.val     = {spmmat generic units pr atl option extn};
aal.help    = {'AUTOMATED ANATOMICAL LABELING toolbox:Reference Automated Anatomical Labeling of Activations in SPM Using a Macroscopic Anatomical Parcellation',... 
                'of the MNI MRI Single-Subject Brain.',...
                'N. Tzourio-Mazoyer, B. Landeau, D. Papathanassiou, F. Crivello, O. Etard, N. Delcroix, B. Mazoyer, and M. Joliot. NeuroImage 2002.'};
aal.prog    = @spm_local_aal;

%======================================================================
function spm_local_aal(job)

if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','batch_aal')); end
batch_aal(job);
