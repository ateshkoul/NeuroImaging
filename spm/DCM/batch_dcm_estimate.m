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
job = varargin{1}
% checking what is the job that has been inputted; may be removed.
DCM = job.data
% for parsing the input
DCM = char(DCM)
% for converting into char format as the default format used is struc

spm_dcm_estimate(DCM)            
% passing the argument to the function
        
               
end
