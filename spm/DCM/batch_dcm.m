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
job = varargin{1};

% Contrast from config file
specify.spmmat = job.spmmat;

% DCM name to be appended (DCM-name)
specify.name = job.name;

% no. of conditions
specify.cond=job.cond;

% directory for ROI
voi = job.vois;

% DCM options
specify.options = job.options;
specify.savedcm = job.savedcm;

for p = 1:numel(job.vois)
    vol{p}=char(job.vois(p).voi);
end

% no. of VOIs
m = numel(voi);
job.model.modcreate.nROIs =m;
specify.voi = vol;

% Cannot use : if exist ('job.model.matrix') directly because of exist
% output; check help exist for details
if isfield (job.model,'matrix')
    % For Select matrix option with only 1 matrix as an input  
    matrix = char(job.model.matrix);
    param = load(matrix);
    specify.matrix.a = param.DCM.a; 
    specify.matrix.b = param.DCM.b; 
    specify.matrix.c = param.DCM.c; 
    specify.matrix.d = param.DCM.d;
    spm_dcm_specify_batch(specify);
else 
    % for generating matrix 
    % send to DCM_model_create and get output from that program
    % with a loop to check how many models are generated.
    models = DCM_model_create(job.model.modcreate);
    nmodel = size(models,2);
    for i=1:nmodel
        a = models(i).a; % a is a jxk matrix r8 now
        b_cell = models(i).b; % b is a 1xn cell containing n jxk matrices r8 now
        c = models(i).c;
        d = models(i).d;
        specify.matrix.a = a;
        specify.matrix.c = c;
        specify.matrix.d = d;
        %% For all the modulatory connectivity matrices for this a, transfer each model to spm_dcm_specify_batch
        s = size(b_cell,2);
        for model_number = 1:s
            if size(b_cell{1,model_number},2)>=1 % If matrix exists
            b(:,:,1) = b_cell{1,model_number};
            b(:,:,2) = b_cell{1,model_number};
            b(:,:,3) = zeros(m);
            specify.matrix.b = b;
            specify.model= strcat(num2str(i),'-',num2str(model_number));
            spm_dcm_specify_batch(specify);
            end
        end
    
    
    end


end
