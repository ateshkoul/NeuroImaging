function batch_dcm_estimate(varargin)

% ______________________________________________________________________
% 
% Batch mode file for AAL 
% The file takes input from Batch toolbox in spm8 
% It is called by config file tbx_cfg_aal
%
% Written by Atesh Koul, National Brain Research Center, India
% ___________________________________________________________________________
% 
% checking what is the job that has been inputted; may be removed.
job = varargin{1};
job.data;
% Number of input DCM
noDCMs = numel(job.data)
for i = 1:noDCMs
% for parsing the input
DCM = job.data{i};
% for converting into char format as the default format used is struc
DCM = char(DCM)
% passing the argument to the function
spm_dcm_estimate(DCM);        
end
end
