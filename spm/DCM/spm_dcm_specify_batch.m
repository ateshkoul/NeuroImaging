function DCM = spm_dcm_specify_batch(specify)
% Specify inputs of a DCM
% FORMAT [DCM] = spm_dcm_specify
%
% DCM  - the DCM structure (see spm_dcm_ui)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_specify.m 4185 2011-02-01 18:46:18Z guillaume $

% Modified for Batch input
% Atesh Koul, National Brain Research Centre, India


%% Get variables
    
    % spm mat file location
    spmmatfile = specify.spmmat;
    % DCM name to be appended
    name = specify.name;
    % VOI location
    voi = specify.voi;
    % Matrix a
    a = specify.matrix.a;
    % Matrix b
    b = specify.matrix.b;
    % Matrix c
    c = specify.matrix.c;
    % Matrix d
    d = specify.matrix.d;
    % RT
    RT = specify.options.delays;
    % TE
    TEu = specify.options.te;
    % DCM options
    nonlin = specify.options.lin;
    two_st = specify.options.state;
    sto = specify.options.sto;
    cen = specify.options.cen;
    ncond = specify.cond;
    if isfield(specify,'model')
            model_number = specify.model;
    end
       
%%


%==========================================================================
% Get design and directory
%==========================================================================
if exist('spmmatfile')
    sts = 1;
    spmmatfile = char(spmmatfile);
else
[spmmatfile, sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
end

if ~sts, DCM = []; return; end
swd = spm_str_manip(spmmatfile,'H');
try
    load(fullfile(swd,'SPM.mat'))
catch
    error(['Cannot read ' fullfile(swd,'SPM.mat')]);
end

%==========================================================================
% Name
%==========================================================================
if ~exist('name')
    name  = spm_input('name for DCM_???.mat','+1','s');
end

%==========================================================================
% Outputs
%==========================================================================

%-Get cell array of region structures
%--------------------------------------------------------------------------
if ~exist('VOI')
    P = voi;
else
P = cellstr(spm_select([1 8],'^VOI.*\.mat$',{'select VOIs'},'',swd));
end
m = numel(P);
for i = 1:m
    p = load(P{i},'xY');
    xY(i) = p.xY;
end

%==========================================================================
% Inputs
%==========================================================================

%-Get (nc) 'causes' or inputs U
%--------------------------------------------------------------------------
% spm_input('Input specification:...  ',1,'d');
Sess = SPM.Sess(xY(1).Sess);
if isempty(Sess.U)
    % spontaneous activity, i.e. no stimuli
    nc = 0;
    U = [];
else
%     % with stimuli
    U.dt   = Sess.U(1).dt;
    u      = length(Sess.U);
    U.name = {};
    U.u    = [];    
    for  i = 1:ncond
        for j = 1:length(Sess.U(i).name)
%             
                U.u = [U.u Sess.U(i).u(33:end,j)];
                U.name{end + 1} = Sess.U(i).name{j};
%             str = ['include ' Sess.U(i).name{j} '?'];
%             if spm_input(str,'+1','y/n',[1 0],1)
%                 U.u             = [U.u Sess.U(i).u(33:end,j)];
%                 U.name{end + 1} = Sess.U(i).name{j};
%             end
%             end
        end
    end
    nc     = size(U.u,2);
end

%==========================================================================
% Timings
%==========================================================================

%spm_input('Timing information:...  ',-1,'d');

%-Slice timings
%--------------------------------------------------------------------------
RT     = SPM.xY.RT;
%delays=1;
if exist('RT')
    delays = repmat(RT,1,m);
else
delays = spm_input('Slice timings [s]','+1','r', repmat(RT,1,m),m,[0 RT]);
end
delays = delays';
%-Echo time (TE) of data acquisition
%--------------------------------------------------------------------------
TEu = TEu/100;
TE    = 0.04;
TE_ok = 0;
while ~TE_ok
    if exist('TEu')
        TE = TEu;
    else
        TE = spm_input('Echo time, TE [s]', '+1', 'r', TE);
    end
    if ~TE || (TE < 0) || (TE > 0.1)
        str = { 'Extreme value for TE or TE undefined.',...
            'Please re-enter TE (in seconds!)'};
        spm_input(str,'+1','bd','OK',[1],1);
    else
        TE_ok = 1;
    end
end

%==========================================================================
% Model options
%==========================================================================
if nc                                                     % there are inputs
    options.nonlinear  = nonlin;
    options.two_state  = two_st;
    options.stochastic = sto;
    options.centre     = cen;
    options.endogenous = 0;  
else 
            options.nonlinear  = 0;
            options.two_state  = 0;
            options.stochastic = 1;
            options.centre     = 1;
            options.endogenous = 1;
end

%==========================================================================
% Graph connections
%==========================================================================
if options.endogenous
    b     = zeros(m,m,1);
    c     = zeros(m,1);
end

% Response
%==========================================================================

%-Response variables & confounds (NB: the data have been whitened)
%--------------------------------------------------------------------------
n     = length(xY);                      % number of regions
v     = length(xY(1).u);                 % number of time points
Y.dt  = SPM.xY.RT;
Y.X0  = xY(1).X0;
for i = 1:n
    Y.y(:,i)  = xY(i).u;
    Y.name{i} = xY(i).name;
end

%-Error precision components (one for each region) - i.i.d. (because of W)
%--------------------------------------------------------------------------
Y.Q        = spm_Ce(ones(1,n)*v);
% 
% 
% %==========================================================================
% % DCM structure
% %==========================================================================

% Endogenous input specification
if isempty(U)
    U.u    = zeros(v,1);
    U.name = {'null'};
end

%-Store all variables in DCM structure
%--------------------------------------------------------------------------
DCM.a       = a;
DCM.b       = b;
DCM.c       = c;
DCM.d       = d;
DCM.U       = U;
DCM.Y       = Y;
DCM.xY      = xY;
DCM.v       = v;
DCM.n       = n;
DCM.TE      = TE;
DCM.delays  = delays;
DCM.options = options;

% Not Implimented: If the user wants to not input the destination folder
% for saving models, this would default to the spmmat directory
if ~isfield(specify,'savedcm')
    savedcm = swd;
else
    savedcm = char(specify.savedcm);
end


%-Save the model
%--------------------------------------------------------------------------
% The if statement for the two conditions - Generate model and Select
% matrix conditions
if isfield(specify,'model')
    if spm_check_version('matlab','7') >= 0
        save(fullfile(savedcm,['DCM_' name '_' model_number '.mat']),'-V6','DCM')
    else
        save(fullfile(savedcm,['DCM_' name '_' model_number '.mat']),'DCM');
    end
else
    if spm_check_version('matlab','7') >= 0
        save(fullfile(savedcm,['DCM_' name '.mat']),'-V6','DCM')
    else
        save(fullfile(savedcm,['DCM_' name '.mat']),'DCM');
    end
    
end






