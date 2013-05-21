function outstructure=peak_nii_inputs(instructure,hdrname,outputargs)
% Checks whether inputs are valid or not.
%   
%   ppi_nii_inputs.v2 last modified by Donald G. McLaren, PhD
%   (mclaren@nmr.mgh.harvard.edu)
%   GRECC, Bedford VAMC
%   Department of Neurology, Massachusetts General Hospital and Havard
%   Medical School
%
% License:
%   Copyright (c) 2011, Donald G. McLaren and Aaron Schultz
%   All rights reserved.
%
%    Redistribution, with or without modification, is permitted provided that the following conditions are met:
%    1. Redistributions must reproduce the above copyright
%        notice, this list of conditions and the following disclaimer in the
%        documentation and/or other materials provided with the distribution.
%    2. All advertising materials mentioning features or use of this software must display the following acknowledgement:
%        This product includes software developed by the Harvard Aging Brain Project.
%    3. Neither the Harvard Aging Brain Project nor the
%        names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
%    4. You are not permitted under this Licence to use these files
%        commercially. Use for which any financial return is received shall be defined as commercial use, and includes (1) integration of all 	
%        or part of the source code or the Software into a product for sale or license by or on behalf of Licensee to third parties or (2) use 	
%        of the Software or any derivative of it for research with the final aim of developing software products for sale or license to a third 	
%        party or (3) use of the Software or any derivative of it for research with the final aim of developing non-software products for sale 
%        or license to a third party, or (4) use of the Software to provide any service to an external organisation for which payment is received.
%
%   THIS SOFTWARE IS PROVIDED BY DONALD G. MCLAREN (mclaren@nmr.mgh.harvard.edu) AND AARON SCHULTZ (aschultz@nmr.mgh.harvard.edu)
%   ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
%   FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
%   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
%   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
%   TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
%   USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%   
%   In accordance with the licences of the atlas sources as being distibuted solely
%   for non-commercial use; neither this program, also soley being distributed for non-commercial use,
%   nor the atlases containe herein should therefore not be used for commercial purposes; for such
%   purposes please contact the primary co-ordinator for the relevant
%   atlas:
%       Harvard-Oxford: steve@fmrib.ox.ac.uk
%       JHU: susumu@mri.jhu.edu
%       Juelich: S.Eickhoff@fz-juelich.de
%       Thalamus: behrens@fmrib.ox.ac.uk
%       Cerebellum: j.diedrichsen@bangor.ac.uk
%       AAL_MNI_V4: maldjian@wfubmc.edu and/or bwagner@wfubmc.edu
%
%   For the program in general, please contact mclaren@nmr.mgh.harvard.edu
%
%   Change Log:
%     4/11/2001: Allows threshold to be -Inf

%% Format input instructure
while numel(fields(instructure))==1
    F=fieldnames(instructure);
    instructure=instructure.(F{1}); %Ignore coding error flag.
end

%% Check for UID
if isfield(instructure,'UID')
    if ischar(instructure.UID)
        outstructure.UID=instructure.UID;
    elseif iscell(instructure.UID)
        outstructure.UID=instructure.UID{1};
    else
        error('Maskname must be a string or cell array of 1 cell')
    end
    if ~strncmp(outstructure.UID,'_',1) && numel(outstructure.UID~=0)
        outstructure.UID=['_' outstructure.UID];
    end
else
    outstructure.UID=datestr(clock,30);
end
%% outfile
try
    outstructure.out=instructure.out;
    if isempty(outstructure.out)
        vardoesnotexist; % triggers catch statement
    end
catch
    [path,file]=fileparts(hdrname);
    if ~isempty(path)
        outstructure.out=[path filesep file '_peaks' outstructure.UID];
    else
        outstructure.out=[file '_peaks' outstructure.UID];
    end
end

%% sign of data
try
    outstructure.sign=instructure.sign;
    if ~strcmpi(outstructure.sign,'pos') && ~strcmpi(outstructure.sign,'neg')
        vardoesnotexist; % triggers catch statement
    end
catch
    outstructure.sign='pos';
    disp('Using default option of positives for sign')
end

%% threshold
try
    outstructure.thresh=instructure.thresh;
    if ~isnumeric(outstructure.thresh) || isempty(outstructure.thresh)
        vardoesnotexist; % triggers catch statement
    end
    if outstructure.thresh<0
        if strcmpi(outstructure.sign,'neg')  
            outstructure.thresh=outstructure.thresh*-1;
        elseif outstructure.thresh==-Inf
        else
            vardoesnotexist; % triggers catch statement
        end
    end
catch
    outstructure.thresh=0;
    disp('Using default option of 0 for threshold')
end

%% statistic type (F or T)
try 
    outstructure.type=instructure.type;
    if ~strcmpi(outstructure.type,'T') && ~strcmpi(outstructure.type,'F') && ~strcmpi(outstructure.type,'none') && ~strcmpi(outstructure.type,'Z')
        vardoesnotexist; % triggers catch statement
    end
catch
    if outstructure.thresh<1 && outstructure.thresh>0
        error(['Statistic must defined using: ' instructure.type])
    else
        outstructure.type='none';
    end
end

%% voxel limit
try
    outstructure.voxlimit=instructure.voxlimit;
    if ~isnumeric(outstructure.voxlimit) || outstructure.voxlimit<0
        vardoesnotexist; % triggers catch statement
    end
catch
    outstructure.voxlimit=1000;
    disp('Using default option of 1000 for voxel limit')
end

%% separation distance for peaks
try
    outstructure.separation=instructure.separation;
    if ~isnumeric(outstructure.separation) || outstructure.separation<0
        vardoesnotexist; % triggers catch statement
    end
catch
    outstructure.separation=8;
    disp('Using default option of 8mm for separate peaks')
end

%% Output peaks or collapse peaks within a cluster (0 collapse peaks closer
% than separation distance, 1 remove peaks closer than separation distance
% to mirror SPM)
try
    outstructure.SPM=instructure.SPM;
    if ~isnumeric(outstructure.SPM) || (outstructure.SPM~=0 && outstructure.SPM~=1)
        vardoesnotexist; % triggers catch statement
    end
catch
    outstructure.SPM=0;
    disp('Using default option of collapsing peaks, this is not similar to SPM')
end
%% Connectivity radius
try
    outstructure.conn=instructure.conn;
    if ~isnumeric(outstructure.conn) || (outstructure.conn~=6 && outstructure.conn~=18 && outstructure.conn~=26)
        vardoesnotexist; % triggers catch statement
    end
catch
    outstructure.conn=18;
    disp('Using default option of face connected voxels, 18-voxel neighborhood')
end
%% Cluster extent threshold
try
    outstructure.cluster=instructure.cluster;
    if ~isnumeric(outstructure.cluster) || outstructure.cluster<0
        vardoesnotexist; % triggers catch statement
    end
catch
    outstructure.cluster=0;
    disp('Using default option of 0 voxels for cluster extent. This is BAD!!!!')
end
%% mask file
try
    outstructure.mask=instructure.mask;
    if ~isempty(outstructure.mask) && ~exist(outstructure.mask,'file')
        vardoesnotexist; % triggers catch statement
    end
catch
    outstructure.mask={};
    disp('Using default option of no mask')
end

%% degrees of freedom numerator
try
    outstructure.df1=instructure.df1;
    if ~isnumeric(outstructure.df1) || outstructure.df1<1
        vardoesnotexist; % triggers catch statement
    end
catch
    if (strcmpi(outstructure.type,'T') || strcmpi(outstructure.type,'F')) && (outstructure.thresh>0 && outstructure.thresh<1)
        error('degrees of freedom numerator must be defined using df1 field; can be gotten from SPM')
    else
    outstructure.df1=[];
    end
end

%% degrees of freedom denominator
try
    outstructure.df2=instructure.df2;
    if ~isnumeric(outstructure.df2) || outstructure.df2<1
        vardoesnotexist; % triggers catch statement
    end
catch
    if (strcmpi(outstructure.type,'F')) && (outstructure.thresh>0 && outstructure.thresh<1)
        error('degrees of freedom numerator must be defined using df2 field; can be gotten from SPM')
    else
    outstructure.df2=[];
    end
end

%% Make threshold a non-decimal
if (strcmpi(outstructure.type,'T') || strcmpi(outstructure.type,'F') || strcmpi(outstructure.type,'Z')) && (outstructure.thresh>0 && outstructure.thresh<1)
    if strcmpi(outstructure.type,'T')
        outstructure.thresh = spm_invTcdf(1-outstructure.thresh,outstructure.df1);
    elseif strcmpi(outstructure.type,'F')
        outstructure.thresh = spm_invFcdf(1-outstructure.thresh,outstructure.df1,outstructure.df2);
    else 
        outstructure.thresh=norminv(1-outstructure.thresh,0,1);
    end
end

%% Check exact field
if ~isfield(instructure,'exact') || ~isnumeric(instructure.exact) || instructure.exact~=1
    outstructure.exact=0;
else
    outstructure.exact=instructure.exact;
    if outstructure.cluster<=0
        error('Cluster must be greater than 0')
    end
    img=spm_read_vols(spm_vol(outstructure.mask));
    if outstructure.cluster>sum(img(:)>0)
        error('Cluster must be smaller than mask')
    end
end

%% Check for maskname
if isfield(instructure,'maskname')
    if ischar(instructure.maskname)
        outstructure.maskname=instructure.maskname;
    elseif iscell(instructure.maskname)
        outstructure.maskname=instructure.maskname{1};
    else
        error('Maskname must be a string or cell array of 1 cell')
    end
    if ~strncmp(outstructure.maskname,'_',1)
        outstructure.maskname=['_' outstructure.maskname];
    end
else
    outstructure.maskname='';
end

%% Save structure
parameters=outstructure;
try
	parameters.label=instructure.label;
end
try
	parameters.nearest=instructure.nearest;
end
save([outstructure.out '_structure.mat'],'parameters','-v7.3'); %After this point, outstructure can be extremely large, thus it is only saved here.
clear parameters

%% Check for Peak Labels
try 
    outstructure.label.source=instructure.label;
    outstructure.label=roilabels(outstructure.label.source);
catch
    try 
        outstructure.label.source=instructure.label.source;
        outstructure.label=roilabels(outstructure.label.source);
    catch
        outstructure.label.source={};
        disp('Using default option of no labels.')
    end
end
try
    outstructure.label.nearest=instructure.nearest;
    if ~isnumeric(outstructure.label.nearest) || (outstructure.label.nearest~=0 && outstructure.label.nearest~=1)
        vardoesnotexist; % triggers catch statement
    end
catch
    outstructure.label.nearest=0;
    disp('Using default option of not using the nearest region label.')
end

% Labels needed?
if (outputargs==2 || outputargs==3) && isempty(outstructure.label.source)
        error('Labels must be specified')
        % Program exits if it gets here.
end
return

%% ROILABELS
function label=roilabels(source)
peak_nii_dir=fileparts(which('peak_nii.m'));
switch source
    case 'aal_MNI_V4'
        regionfile=[peak_nii_dir filesep 'aal_MNI_V4.img']; % From WFU_pickatlas
        load([peak_nii_dir filesep 'aal_MNI_V4_List.mat']);
    case 'Nitschke_Lab'
        regionfile=[peak_nii_dir filesep 'ControlabilityMask.nii']; % Provided by Deb Kerr
        load([peak_nii_dir filesep 'ControlabilityMask_List.mat']);
    case 'JHU_tracts'
        regionfile=[peak_nii_dir filesep 'JHU-ICBM-tracts-maxprob-thr0-1mm.nii']; % From FSL Atlas Files
        load([peak_nii_dir filesep 'JHU_tract_labels.mat']);
    case 'JHU_whitematter'
        regionfile=[peak_nii_dir filesep 'JHU-WhiteMatter-labels-1mm.nii']; % From FSL Atlas Files
        load([peak_nii_dir filesep 'JHU_labels.mat']);
    case 'Thalamus' 
        regionfile=[peak_nii_dir filesep 'Thalamus-maxprob-thr0-1mm.nii']; % From FSL Thalamus Atlas Files
        load([peak_nii_dir filesep 'Thalamus_labels.mat']);
    case 'Talairach'
        regionfile=[peak_nii_dir filesep 'Talairach-labels-1mm.nii']; % From FSL Talairach Atlas Files
        load([peak_nii_dir filesep 'Talairach_Labels.mat']);
    case 'MNI'
        regionfile=[peak_nii_dir filesep 'MNI-maxprob-thr0-1mm.nii']; % From FSL MNI Atlas Files
        load([peak_nii_dir filesep 'MNI_labels.mat']);
    case 'HarvardOxford_cortex'
        regionfile=[peak_nii_dir filesep 'HarvardOxford-cort-maxprob-thr0-1mm.nii']; % From FSL Atlas Files
        load([peak_nii_dir filesep 'HarvardOxford_cortical_labels']);
    %case 'HarvardOxford_subcortical' % Labels do not match image
    %    regionfile=[peak_nii_dir filesep 'HarvardOxford-sub-maxprob-thr0-1mm.nii']; % From FSL Atlas Files
    %    load([peak_nii_dir filesep 'HarvardOxford_subcortical_labels']);
    case 'Juelich'
        regionfile=[peak_nii_dir filesep 'Juelich-maxprob-thr0-1mm.nii']; % From FSL Atlas Files
        load([peak_nii_dir filesep 'Juelich_labels.mat']);
    case 'Cerebellum-flirt'
        regionfile=[peak_nii_dir filesep 'Cerebellum-MNIflirt-maxprob-thr0-1mm.nii']; % From FSL Atlas Files
        load([peak_nii_dir filesep 'Cerebellum_labels.mat']);
    case 'Cerebellum-fnirt'
        regionfile=[peak_nii_dir filesep 'Cerebellum-MNIfnirt-maxprob-thr0-1mm.nii']; % From FSL Atlas Files
        load([peak_nii_dir filesep 'Cerebellum_labels.mat']);
    %case 'Custom'
    %    regionfile=['ROI image file'];
    %    load(['ROIlist mat-file']);
end

label.source=source;
label.rf.hdr=spm_vol(regionfile);
[label.rf.img,label.rf.XYZmm]=spm_read_vols(label.rf.hdr);
label.ROI=ROI;
label.ROInames={ROI.Nom_C}; %ROInames is taken from ROI, where ROI is structure
if ~isfield(ROI,'ID')
		 display('ERROR: ROI must be a structure with an ID field that is the ID of values in region file')
		 return
end
