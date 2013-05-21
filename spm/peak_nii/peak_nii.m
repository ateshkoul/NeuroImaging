function [voxels regions mapparameters UID]=peak_nii(image,mapparameters)
%%
% peak_nii will write out the maximum T (or F) of the local maxima that are
% not closer than a specified separation distance.  
% SPM=0: Those that are closer are collapsed based on the COG using number 
%   of voxels at each collapsing point. The maximum T 
%   (or F) is retained. This program should be similar to peak_4dfp in use at
%   WashU (although I haven't seen their code).
% SPM=1: Eliminates the peaks closer than a specified distance to mimic
%   result tables.
%
% INPUTS:
% image string required. This should be a nii or img file.
% mapparameters is either a .mat file or a pre-load structure with the
% following fields:
%           out: output prefix, default is to define using imagefile
%          sign: 'pos' or 'neg', default is 'pos' NOTE: only can do one
%                direction at a time
%          type: statistic type, 'T' or 'F' or 'none'
%      voxlimit: number of peak voxels in image
%    separation: distance to collapse or eliminate peaks
%           SPM: 0 or 1, see above for details
%          conn: connectivity radius, either 6,18, or 26
%       cluster: cluster extent threshold in voxels
%          mask: optional to mask your data
%           df1: numerator degrees of freedom for T/F-test (if 0<thresh<1)
%           df2: denominator degrees of freedom for F-test (if 0<thresh<1)
%       nearest: 0 or 1, 0 for leaving some clusters/peaks undefined, 1 for finding the
%                nearest label
%         label: optional to label clusters, options are 'aal_MNI_V4';
%                'Nitschke_Lab'; FSL ATLASES: 'JHU_tracts', 'JHU_whitematter',
%                'Thalamus', 'Talairach', 'MNI', 'HarvardOxford_cortex', 'Cerebellum-flirt', 'Cerebellum-fnirt', and 'Juelich'. 
%                'HarvardOxford_subcortical' is not available at this time because
%                the labels don't match the image.
%                Other atlas labels may be added in the future
%        thresh: T/F statistic or p-value to threshold the data or 0
%         exact: optional. 1 if the cluster should be shrunk/expand to cluster, 0
%                do not shrink/expand cluster.
%      maskname: if your using a mask, you have the option to specify a
%                maskname to be included in the filenames
%           UID: unique ID. If this field is not specified, the UID will be
%                a timestamp based on the computer clock. If specified, it
%                must be a string. If you do not want to use a UID, then set
%                this field to ''. NOTE: By specifying a UID, there is the
%                potential to overwrite older files. UID is used to avoid
%                listing all parameters in the filenames.
%
% OUTPUTS:
%   voxels  -- table of peaks
%       cell{1}-
%         col. 1 - Cluster size
%         col. 2 - T/F-statistic
%         col. 3 - X coordinate
%         col. 4 - Y coordinate
%         col. 5 - Z coordinate
%         col. 6 - number of peaks collapsed
%         col. 7 - sorted cluster number
%       cell{2}- region names
%   regions -- region of each peak -- optional
%   mapparameters -- see above
%   UID   -- unique ID value -- user defined in mapparameters.UID or if
%   field does not exist, then this will be the datestr.
%
% NIFTI FILES SAVED:
%   *_maskname_clusters.nii:                     
%                               contains the clusters and their numbers (column 7)
%   (image)_peaks_UID_thresh*_extent*_maskname.nii:             
%                               contains the thresholded data
%   (image)_peaks_UID_thresh*_extent*_maskname_peaknumber.nii:   
%                               contains the peaks of the data,
%                               peaks are numbered by their order
%                               in the table (voxels)
%   (image)_peaks_UID_thresh*_extent*_maskname_peakcluster.nii:  
%                               contains the peaks of the data,
%                               peaks are numbered by their cluster (column 7)
%   *(image) is the image name with the the path or extension
%   NOTE: if maskname or UID are empty, then _UID and _maskname will be
%   dropped from filenames.
%
% MAT-FILES SAVED:
%   Peak_(image)_peaks_UID.mat:        contains voxelsT variable and regions, if applicable 
%   (image)_peaks_UID_structure:       contains parameter variable with
%                                       parameters used
%   *(image) is the image name with the the path or extension
%   NOTE: if maskname or UID are empty, then _UID and _maskname will be
%   dropped from filenames.
%
%
% EXAMPLE: voxels=peak_nii('imagename',mapparameters)
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
%   peak_nii.v3 -- Last modified on 12/10/2010 by Donald G. McLaren, PhD
%   (mclaren@nmr.mgh.harvard.edu)
%   Wisconsin Alzheimer's Disease Research Center - Imaging Core, Univ. of
%   Wisconsin - Madison
%   Neuroscience Training Program and Department of Medicine, Univ. of
%   Wisconsin - Madison
%   GRECC, William S. Middleton Memorial Veteren's Hospital, Madison, WI
%   GRECC, Bedford VAMC
%   Department of Neurology, Massachusetts General Hospital and Havard
%   Medical School
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

%% Program begins here
try
    if ~strcmp(spm('Ver'),'SPM8')
        disp('PROGRAM ABORTED:')
        disp('  You must use SPM8 to process your data; however, you can use SPM.mat files')
        disp('  generated with SPM2 or SPM5. In these cases, simply specify the option SPMver')
        disp('  in single qoutes followed by a comma and the version number.')
        disp(' ')
        disp('Make sure to add SPM8 to your MATLAB path before re-running.')
        return
    else
        addpath(fileparts(which('spm')))
    end
catch
    disp('PROGRAM ABORTED:')
    disp('  You must use SPM8 to process your data; however, you can use SPM.mat files')
    disp('  generated with SPM2 or SPM5. In these cases, simply specify the option SPMver')
    disp('  in single qoutes followed by a comma and the version number.')
    disp(' ')
    disp('Make sure to add SPM8 to your MATLAB path before re-running.')
    return
end

%% Check inputs
if exist(image,'file')==2
    I1=spm_vol(image);
    infoI1=I1;
    [I1,voxelcoord]=spm_read_vols(I1);
    if nansum(nansum(nansum(abs(I1))))==0
        error(['Error: ' image ' is all zeros or all NaNs'])        
    end
else
    error(['File ' image ' does not exist'])
end
if nargin==2
    if ischar(mapparameters) && exist(mapparameters,'file')==2
        mapparameters=load(mapparameters);
    end
    if ~isstruct(mapparameters)
        error('Mapparameters is not a structure OR not a file that contains a structure');
    end
else
   error('Mapparameters must be specified as a file or a structure');
end

invar=peak_nii_inputs(mapparameters,infoI1.fname,nargout);
UID=invar.UID;

%% Read in data and mask (if available)
if strcmpi(invar.sign,'neg')
    I1=-1.*I1;
    disp(['Threshold is:-' num2str(invar.thresh)])
    invar.thresh2=invar.thresh*-1;
else
    disp(['Threshold is:' num2str(invar.thresh)])
    invar.thresh2=invar.thresh;
end
mapparameters=invar;

if ~isempty(invar.mask)
    I2=spm_vol(invar.mask);
    infoI2=I2;
    if infoI1.mat==infoI2.mat
        if infoI1.dim==infoI2.dim
            I2=spm_read_vols(I2);
        else
            Vi(1)=infoI1; Vi(2)=infoI2;
            Vo=Vi(1); Vo.pinfo=Vi(2).pinfo; Vo.fname='tmp.nii'; Vo.pinfo=[1 0 352]'; Vo.dt=[spm_type('float32') spm_platform('bigend')];
            Vo=spm_imcalc(Vi,Vo,'i2');
            I2=spm_read_vols(spm_vol('tmp.nii'));
        end
    else
        Vi(1)=infoI1; Vi(2)=infoI2;
        Vo=Vi(1); Vo.pinfo=Vi(2).pinfo; Vo.fname='tmp.nii'; Vo.pinfo=[1 0 352]'; Vo.dt=[spm_type('float32') spm_platform('bigend')];
        spm_imcalc(Vi,Vo,'i2',{[],[],0});
        I2=spm_read_vols(spm_vol('tmp.nii'));
    end
            
    I2=(I2>0)+(I2<0); % Ensures that the mask is binary
    I2(I2==0)=NaN; % Since 0 can be used in computations, convert 0s in mask to NaN.
    I=I1.*I2; % Mask Dataset
    %% Mask Region file
    if ~isempty(invar.label.source)
        try
            invar.label.rf.img=invar.label.rf.img.*(I2~=0);
        catch
            Vi(1)=infoI1; Vi(2)=invar.label.rf.hdr;
            Vo=Vi(1); Vo.fname='tmp.nii'; Vo.pinfo=[1 0 352]'; Vo.dt=[spm_type('float32') spm_platform('bigend')];
            spm_imcalc(Vi,Vo,'i2',{[],[],0});
            invar.label.rf.hdr=spm_vol('tmp.nii');
            [invar.label.rf.img, invar.label.rf.XYZmm]=spm_read_vols(invar.label.rf.hdr);
            invar.label.rf.img=invar.label.rf.img.*(I2>0); % Mask Region File
        end
    end
else
    I=I1;
end

%% Program begins here
% Find significant voxels
ind=find(I>invar.thresh);
if isempty(ind)
    voxels=[]; regions={};
    disp(['NO MAXIMA ABOVE ' num2str(invar.thresh) '.'])
    if invar.exact==1
        disp('To find the cluster in this subject, please set thresh to 0')
        return
%         ind=find(I>0);
%         if isempty(ind)
%             voxels=[]; regions={};
%             disp(['NO MAXIMA ABOVE ' 0 '.'])
%             disp('This should not happen. Hmm.')
%             return
%         end
    else
        return
    end
else
   [L(1,:),L(2,:),L(3,:)]=ind2sub(infoI1.dim,ind);
end

% Cluster signficant voxels
A=peakcluster(L,invar.conn,infoI1); % A is the cluster of each voxel
A=transpose(A);
n=hist(A,1:max(A));
for ii=1:size(A,1)
    if n(A(ii))<invar.cluster % removes clusters smaller than extent threshold
        A(ii,1:2)=NaN;
    else
        A(ii,1:2)=[n(A(ii)) A(ii,1)];
    end
end

% Combine A (cluster labels) and L (voxel indicies)
L=L';
A(:,3:5)=L(:,1:3);
% Remove voxels that are in small clusters
A(any(isnan(A),2),:) = [];

% Stop if no significant clusters
if isempty(A)
    voxels=[]; regions={};
    display(['NO CLUSTERS LARGER THAN ' num2str(invar.cluster) ' voxels.'])
    if invar.exact==1
        disp(['The largest cluster in this subject @ ' num2str(invar.thresh2) 'is ' num2str(max(n)) 'voxels.'])
        disp('To find the cluster in this subject, please change the cluster size or threshold.')
    end
    return
end


% Save clusters
T=peakcluster(transpose(A(:,3:5)),invar.conn,infoI1,[invar.out '_thresh' num2str(invar.thresh2) '_extent' num2str(invar.cluster) invar.maskname]);
A(:,2)=T(:,1); clear T
% Save significant data
Iclust=spm_read_vols(spm_vol([invar.out '_thresh' num2str(invar.thresh2) '_extent' num2str(invar.cluster) invar.maskname '_clusters.nii']));
if strcmpi(invar.sign,'neg')
    Ithresh=-1.*I.*(Iclust>0);
else
    Ithresh=I.*(Iclust>0);
end
out=infoI1;
out.fname=[invar.out '_thresh' num2str(invar.thresh2) '_extent' num2str(invar.cluster) invar.maskname '.nii'];
out.descrip=['Thresholded Map @ thresh ' num2str(invar.thresh2) ' and cluster extent ' num2str(invar.cluster) ' in ' invar.maskname];
spm_write_vol(out,Ithresh);

% Find all peaks, only look at current cluster to determine the peak
Ic=zeros(infoI1.dim(1),infoI1.dim(2),infoI1.dim(3),max(A(:,2)));
for ii=1:max(A(:,2))
    Ic(:,:,:,ii)=I.*(Iclust==ii);
end
N=0;
voxelsT=zeros(size(A,1),7);
for ii=1:size(A,1)
    if A(ii,3)==1 || A(ii,4)==1 || A(ii,5)==1 || A(ii,3)==size(Ic,1) || A(ii,4)==size(Ic,2) || A(ii,5)==size(Ic,3)
    else
        if I(A(ii,3),A(ii,4),A(ii,5))==max(max(max(Ic(A(ii,3)-1:A(ii,3)+1,A(ii,4)-1:A(ii,4)+1,A(ii,5)-1:A(ii,5)+1,A(ii,2)))))
            N=N+1;
            voxind=sub2ind(infoI1.dim,A(ii,3),A(ii,4),A(ii,5));
            voxelsT(N,1)=A(ii,1);
            voxelsT(N,2)=I(voxind);
            voxelsT(N,3)=voxelcoord(1,voxind);
            voxelsT(N,4)=voxelcoord(2,voxind);
            voxelsT(N,5)=voxelcoord(3,voxind);
            voxelsT(N,6)=1;
            voxelsT(N,7)=A(ii,2);
        end
    end
end

%Remove empty rows
voxelsT=voxelsT(any(voxelsT'),:);

%Check number of peaks
if size(voxelsT,1)>invar.voxlimit
    voxelsT=sortrows(voxelsT,-2);
    voxelsT=voxelsT(1:invar.voxlimit,:); % Limit peak voxels to invar.voxlimit
end

% Sort table by cluster w/ max T then by T value within cluster (negative
% data was inverted at beginning, so we are always looking for the max).
uniqclust=unique(voxelsT(:,7));
maxT=zeros(length(uniqclust),2);
for ii=1:length(uniqclust)
    maxT(ii,1)=uniqclust(ii);
    maxT(ii,2)=max(voxelsT(voxelsT(:,7)==uniqclust(ii),2));
end
maxT=sortrows(maxT,-2);
for ii=1:size(maxT,1)
    voxelsT(voxelsT(:,7)==maxT(ii,1),8)=ii;
end
voxelsT=sortrows(voxelsT,[8 -2]);
[cluster,uniq,ind]=unique(voxelsT(:,8)); % get rows of each cluster
if invar.exact
    voxelsT(2:end,:)=[];
    A=A(A(:,2)==voxelsT(1,7),:); % Keeps most significant cluster
    voxind=zeros(size(A,1),1);
    for ii=1:size(A,1)
        voxind(ii)=sub2ind(infoI1.dim,A(ii,3),A(ii,4),A(ii,5));
        A(ii,6)=I(voxind(ii));
    end
    B=A; %in case eroding doesn't work
    while size(A,1)>invar.cluster
       A(A(:,6)==min(A(:,6)),:)=[]; 
    end
    
    %check for a single cluster
    newclust=peakcluster(A(:,3:5)',invar.conn,infoI1);
    if max(newclust)>1
         A=B; clear B;
         indmax=find(A(:,6)==max(A(:,6)));
         cluster=zeros(invar.cluster,1);
         cluster(1)=voxind(indmax);
         %voxind is the voxel indices
         indsearch=[];
         for ii=1:(invar.cluster-1)
             a={[A(indmax,3)-1:A(indmax,3)+1],[A(indmax,4)-1:A(indmax,4)+1],[A(indmax,5)-1:A(indmax,5)+1]};
             [xx yy zz]=ndgrid(a{:});
             possind=sub2ind(infoI1.dim,xx(:),yy(:),zz(:));
             addind=voxind(ismember(voxind,possind));
             indsearch=[indsearch addind']; clear addind; %#ok<AGROW>
             indsearch=setdiff(indsearch,cluster);
             indmax=find(A(:,6)==max(I(indsearch)));
             cluster(ii+1)=voxind(indmax);
         end
         A=A(ismember(voxind,cluster),:);
    end
    % Save clusters
    peakcluster2(A,voxelsT,infoI1,[invar.out '_thresh' num2str(invar.thresh2) '_extent' num2str(invar.cluster) invar.maskname]);
    % Save significant data
    Iclust=spm_read_vols(spm_vol([invar.out '_thresh' num2str(invar.thresh2) '_extent' num2str(invar.cluster) invar.maskname '_clusters.nii']));
    if strcmpi(invar.sign,'neg')
        Ithresh=-1.*I.*(Iclust>0);
    else
        Ithresh=I.*(Iclust>0);
    end
    out=infoI1;
    out.fname=[invar.out '_thresh' num2str(invar.thresh2) '_extent' num2str(invar.cluster) invar.maskname '.nii'];
    out.descrip=['Thresholded Map @ thresh ' num2str(invar.thresh2) ' and cluster extent ' num2str(invar.cluster) ' in ' invar.maskname];
    spm_write_vol(out,Ithresh);
    if strcmpi(invar.sign,'neg')
        voxelsT(:,2)=-1*voxelsT(:,2);
    end
    voxelsT(:,7)=[];
    voxels={voxelsT};
    regions={};
    savepeaks(invar,regions,voxelsT)
    return
end
  

%Collapse or eliminate peaks closer than a specified distance
voxelsF=zeros(size(voxelsT,1),size(voxelsT,2));
nn=[1 zeros(1,length(cluster)-1)];
for numclust=1:length(cluster)
    Distance=eps;
    voxelsC=voxelsT(ind==numclust,:);
    while min(min(Distance(Distance>0)))<invar.separation
            [voxelsC,Distance]=vox_distance(voxelsC);
            minD=min(min(Distance(Distance>0)));
            if minD<invar.separation
               min_ind=find(Distance==(min(min(Distance(Distance>0)))));
               [ii,jj]=ind2sub(size(Distance),min_ind(1));
               if invar.SPM==1
                    voxelsC(ii,:)=NaN; % elimate peak
               else
                    voxelsC(jj,1)=voxelsC(jj,1);
                    voxelsC(jj,2)=voxelsC(jj,2);
                    voxelsC(jj,3)=((voxelsC(jj,3).*voxelsC(jj,6))+(voxelsC(ii,3).*voxelsC(ii,6)))/(voxelsC(jj,6)+voxelsC(ii,6)); % avg coordinate
                    voxelsC(jj,4)=((voxelsC(jj,4).*voxelsC(jj,6))+(voxelsC(ii,4).*voxelsC(ii,6)))/(voxelsC(jj,6)+voxelsC(ii,6)); % avg coordinate
                    voxelsC(jj,5)=((voxelsC(jj,5).*voxelsC(jj,6))+(voxelsC(ii,5).*voxelsC(ii,6)))/(voxelsC(jj,6)+voxelsC(ii,6)); % avg coordinate
                    voxelsC(jj,6)=voxelsC(jj,6)+voxelsC(ii,6);
                    voxelsC(jj,7)=voxelsC(jj,7);
                    voxelsC(jj,8)=voxelsC(jj,8);
                    voxelsC(ii,:)=NaN; % eliminate second peak
               end
               voxelsC(any(isnan(voxelsC),2),:) = [];
            end
    end
    try
        nn(numclust+1)=nn(numclust)+size(voxelsC,1);
    end
    voxelsF(nn(numclust):nn(numclust)+size(voxelsC,1)-1,:)=voxelsC;
end
voxelsT=voxelsF(any(voxelsF'),:);
clear voxelsF voxelsC nn

% Label Peaks
if ~isempty(invar.label.source)
    regions=cell(size(voxelsT,1),2);
    for i=1:size(voxelsT,1)
             [regions{i,1} regions{i,2}]=regionname(voxelsT(i,:),invar.label.rf,invar.label.ROI,invar.label.ROInames,invar.label.nearest);
    end
end

% Modify T-values for negative
if strcmpi(invar.sign,'neg')
    voxelsT(:,2)=-1*voxelsT(:,2);
end

% Output an image of the peak coordinates (peak number and cluster number)
Iclusthdr=spm_vol([invar.out '_thresh' num2str(invar.thresh2) '_extent' num2str(invar.cluster) invar.maskname '_clusters.nii']);
[Iclust, Ixyz]=spm_read_vols(Iclusthdr);
Ipeak=Iclust.*0;
Ipeak2=Ipeak;

%Try to determine orientation
oris = [[Iclusthdr.mat(1,1) Iclusthdr.mat(2,2) Iclusthdr.mat(3,3)];...
        [Iclusthdr.mat(1,1) Iclusthdr.mat(3,2) Iclusthdr.mat(2,3)];...
        [Iclusthdr.mat(2,1) Iclusthdr.mat(1,2) Iclusthdr.mat(3,3)];...
        [Iclusthdr.mat(2,1) Iclusthdr.mat(3,2) Iclusthdr.mat(1,3)];...
        [Iclusthdr.mat(3,1) Iclusthdr.mat(1,2) Iclusthdr.mat(2,3)];...
        [Iclusthdr.mat(3,1) Iclusthdr.mat(2,2) Iclusthdr.mat(1,3)]];
ori = find(mean(abs(oris)>.000001',2) == 1);
if numel(ori)>1 % This part has not been fully tested
    oristmp=oris;
    tmp(1)=max(abs(diff([1 1 1 1; 2 1 1 1]*Iclusthdr.mat')));
    oristmp(abs(oris)==tmp(1))=1;
    tmp(2)=max(abs(diff([1 1 1 1; 1 2 1 1]*Iclusthdr.mat')));
    oristmp(abs(oris)==tmp(2))=1;
    tmp(3)=max(abs(diff([1 1 1 1; 1 1 2 1]*Iclusthdr.mat')));
    oristmp(abs(oris)==tmp(3))=1;
    ori = find(mean(oristmp,2)==1);
end

%Determine cross hair size, if value is infite, crosshairs are 2 voxels
a(1)=ceil(2*2/abs(oris(ori,1)));
a(2)=ceil(2*2/abs(oris(ori,2)));
a(3)=ceil(2*2/abs(oris(ori,3)));
a(~isfinite(a))=2;

for ii=1:size(voxelsT,1)
    [Ipeakxyz,Ipeakind] = spm_XYZreg('NearestXYZ',voxelsT(ii,3:5),Ixyz);
    [Ix,Iy,Iz]=ind2sub(size(Ipeak),Ipeakind);
    if (Ix-a(1))<=0
        Ipeak(1:Ix+a(1),Iy,Iz)=ii;
    else
        Ipeak(Ix-a(1):Iclusthdr.dim(1),Iy,Iz)=ii;
    end
    if (Iy-a(2))<=0
        Ipeak(Ix,1:Iy+a(2),Iz)=ii;
    else
        Ipeak(Ix,Iy-a(2):Iclusthdr.dim(2),Iz)=ii;
    end
    if (Iz-a(3))<=0
        Ipeak(Ix,Iy,1:Iz+a(3))=ii;
    else
        Ipeak(Ix,Iy,Iz-a(3):Iclusthdr.dim(3))=ii;
    end
    if (Ix-a(1))<=0
        Ipeak2(1:Ix+a(1),Iy,Iz)=voxelsT(ii,8);
    else
        Ipeak2(Ix-a(1):Iclusthdr.dim(1),Iy,Iz)=voxelsT(ii,8);
    end
    if (Iy-a(2))<=0
        Ipeak2(Ix,1:Iy+a(2),Iz)=voxelsT(ii,8);
    else
        Ipeak2(Ix,Iy-a(2):Iclusthdr.dim(2),Iz)=voxelsT(ii,8);
    end
    if (Iz-a(3))<=0
        Ipeak2(Ix,Iy,1:Iz+a(3))=voxelsT(ii,8);
    else
        Ipeak2(Ix,Iy,Iz-a(3):Iclusthdr.dim(3))=voxelsT(ii,8);
    end
end
out=infoI1;
out.fname=[invar.out '_thresh' num2str(invar.thresh2) '_extent' num2str(invar.cluster) invar.maskname '_peaknumber.nii'];
out.descrip=['Peaks of Thresholded Map @ thresh ' num2str(invar.thresh2) ' and cluster extent ' num2str(invar.cluster) ' in ' invar.maskname];
out.pinfo(1)=1;
spm_write_vol(out,Ipeak);
out=infoI1;
out.pinfo(1)=1;
out.fname=[invar.out '_thresh' num2str(invar.thresh2) '_extent' num2str(invar.cluster) invar.maskname '_peakcluster.nii'];
out.descrip=['Peaks of Thresholded Map @ thresh ' num2str(invar.thresh2) ' and cluster extent ' num2str(invar.cluster) ' in ' invar.maskname];
spm_write_vol(out,Ipeak2);
peakcluster2(A,voxelsT,infoI1,[invar.out '_thresh' num2str(invar.thresh2) '_extent' num2str(invar.cluster) invar.maskname]); %outputs revised cluster numbers
voxelsT(:,7)=[];

savepeaks(invar,regions,voxelsT);

try
    voxels={voxelsT {regions{:,2}}'};
catch
    voxels={voxelsT};
end
return

%% Embedded functions

%% Save Peaks
function savepeaks(invar,regions,voxelsT)
invar.datemod=date;
[path,file,ext]=fileparts(invar.out);
if ~isempty(path)
   save([path filesep 'Peak_' file invar.UID '_thresh' num2str(invar.thresh2) '_extent' num2str(invar.cluster) invar.maskname '.mat'],'regions', 'voxelsT','invar')
else
    save(['Peak_' file invar.UID '_thresh' num2str(invar.thresh2) '_extent' num2str(invar.cluster) invar.maskname '.mat'],'regions', 'voxelsT','invar')
end
return

%% Vox_distance
function [N,Distance] = vox_distance(voxelsT)
% vox_distance compute the distance between local maxima in an image
% The input is expected to be an N-M matrix with columns 2,3,4 being X,Y,Z
% coordinates
%
% pdist is only available with Statistics Toolbox in recent versions of
% MATLAB, thus, the slower code is secondary if the toolbox is unavailable.
% Speed difference is dependent on cluster sizes, 3x at 1000 peaks.
N=sortrows(voxelsT,-1);
try
    Distance=squareform(pdist(N(:,3:5)));
catch
    Distance = zeros(size(N,1),size(N,1));
    for ii = 1:size(N,1);
        TmpD = zeros(size(N,1),3);
        for kk = 1:3;
            TmpD(:,kk) = (N(:,kk+2)-N(ii,kk+2)).^2;
        end
        TmpD = sqrt(sum(TmpD,2));
        Distance(:,ii) = TmpD;
    end
end
%Distance=zeros(length(N(:,1)),length(N(:,1)))*NaN;
%for ii=1:length(N(:,1))
%    for jj=ii+1:length(N(:,1))
%           Distance(ii,jj)=((N(jj,2)-N(ii,2)).^2)+((N(jj,3)-N(ii,3)).^2)+((N(jj,4)-N(ii,4)).^2);
%    end
%end
return

%% Peakcluster
function A=peakcluster(L,conn,infoI1,out)
dim = infoI1.dim;
vol = zeros(dim(1),dim(2),dim(3));
indx = sub2ind(dim,L(1,:)',L(2,:)',L(3,:)');
vol(indx) = 1;
[cci,num] = spm_bwlabel(vol,conn);
A = cci(indx');
if nargin==4
    infoI1.fname=[out '_clusters.nii'];
    infoI1.descrip='clusters';
    infoI1.pinfo=[1 0 0]';
    A=transpose(A);
    L=transpose(L);
    A(:,2:4)=L(:,1:3);
    vol=zeros(dim(1),dim(2),dim(3));
    for ii=1:size(A,1)
        vol(A(ii,2),A(ii,3),A(ii,4))=A(ii,1);
    end
    spm_write_vol(infoI1,vol);
end
return


    

%% Regionname
function [ROInum ROIname]=regionname(voxel,rf,ROI,ROInames,nearest)
if nearest==0
        [xyz,ii] = spm_XYZreg('NearestXYZ',voxel(3:5),rf.XYZmm); % use all voxels
        try
            ROInum=rf.img(ii);
        catch
            ROInum=0;
        end
else    
        nz_ind=find(rf.img>0);
        [xyz,j] = spm_XYZreg('NearestXYZ',voxel(3:5),rf.XYZmm(:,nz_ind)); % use only voxels with a region greater than 0.
        try
            ii=nz_ind(j);
            if isempty(ii)
                invovkecatchstatement
            end
            [junk,D]=vox_distance([0 0 xyz';voxel(1:5)]); % need to pad xyz, voxel is made to be 1x4 with 2-4 being coordianates
        catch
            D(1,2)=Inf;
        end
        if D(1,2)>8 % further than 5 mm from region
           ROInum=0;
        else
           ROInum=rf.img(ii);
        end
end
if ROInum~=0
    try
       ROIind=find([ROI.ID]==ROInum);
       ROIname=ROInames{ROIind};
    catch
        keyboard
    end
else
       ROIname='undefined';
end
return
%% Peakcluster
function peakcluster2(A,voxelsT,infoI1,out)
dim = infoI1.dim;
vol = zeros(dim(1),dim(2),dim(3));
clusters=unique(voxelsT(:,8));
A(:,6)=-1; %eliminated because not in the top peaks
for ii=1:numel(clusters)
   firstvoxel=find(voxelsT(:,8)==ii);
   A(A(:,2)==voxelsT(firstvoxel(1),7),6)=ii;
end
if nargin==4
    infoI1.fname=[out '_clusters.nii'];
    infoI1.descrip='clusters';
    infoI1.pinfo=[1 0 0]';
    for ii=1:size(A,1)
        vol(A(ii,3),A(ii,4),A(ii,5))=A(ii,6);
    end
    spm_write_vol(infoI1,vol);
end
return