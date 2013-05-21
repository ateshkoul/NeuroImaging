function batch_dcm(varargin)

% ______________________________________________________________________
% 
% Batch mode file for AAL 
% The file takes input from Batch toolbox in spm8 
% It is called by config file tbx_cfg_aal
%
% Written by Atesh Koul, National Brain Research Center, India
% ___________________________________________________________________________
% 
job = varargin{1}
% Contrast from config file
spmmat = job.spmmat
% ROI title
% directory for ROI
name = job.name
a = job.a
job.effects
j = numel(job.effects)

for k = 1:j
    b(:,:,k) = job.effects(k).b
    c(:,k) = job.effects(k).c
end
job.rois(1).roi
job.rois(2).roi
d = job.rois.d
roi = job.rois.roi
RT = job.options.delays
TE = job.options.te
TE = TE/100
lin = job.options.lin
state= job.options.state
sto = job.options.sto
cen = job.options.cen


for p = 1:numel(roi)
    voi{p}=char(job.rois(p).roi)
end

m = numel(roi);

% if exist('job.rois.d')
%    for p = 1:m
%        d(:,:,p) = job.effects(p).d;
%    end
% else
%     d = zeros(k,k,0);
% end
    
       spm_dcm_specify_batch(spmmat,name,voi,a,b,c,d,RT,TE,lin,state,sto,cen,j)            
        
               
end
