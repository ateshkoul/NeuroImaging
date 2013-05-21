function batch_peak_nii(job)
% ______________________________________________________________________
% 
% Batch mode file for AAL 
% The file takes input from Batch toolbox in spm8 
% It is called by config file tbx_cfg_aal
%
% Written by Atesh Koul, National Brain Research Center, India
% ___________________________________________________________________________
% 

images = job.analysisinfo.image
mapparameters.UID = job.analysisinfo.uid
mapparameters.out= job.param.out
mapparameters.sign = job.param.sign
mapparameters.thresh=job.param.thresh
mapparameters.type=job.param.tyoe
mapparameters.voxlimit=job.param.voxlimit
mapparameters.separation=job.param.separation
mapparameters.SPM=job.param.SPM
mapparameters.conn=job.param.conn
mapparameters.cluster=job.param.cluster
mapparameters.df1=job.param.df1
mapparameters.df2=job.param.df2
mapparameters.mask=job.param.mask
mapparameters.exact=job.param.exact
mapparameters.sphere=job.param.sphere
mapparameters.clustersphere=job.param.clustersphere


[a b c d]= peak_nii(images,mapparameters)

end

