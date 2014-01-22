function DCM_model = DCM_model_create(job)
% Script to create all possible DCM models for a given no. of ROIs resticted by user given conditions
% NOTE - Works only when modulatory effects enter all at once [Ng,Ar,En,Hi]
% Author - Vaibhav Tyagi, Atesh koul, National Brain Research Centre from Sep 3-5, 2013

% Inputs to the function- 
% A struct with the following variables:
% nROIs - No. Of ROIs in the model as nROIs variable.
% dest - Destination for saved models.
% nConds - No. of conditions in the models.
% intmatrix - A general intrinsic matrix of connections.
% modmatrix - A general modulatory matrix.
% inpmatrix - Direct Input matrix.

% Outputs from the function: A struct DCM model with following variables:
% a - 'a' cell (1 x n), each of n columns contain one possible intrinsic connectivity model/matrix
% b - 'b' cell (m x n), each m row contains 1 or 2 or 3 or....n modulation effect matrices for each intrinsic connectivity model/matrix
% c - 'c' matrix (m x n), each m row depicts ROIs and each n column depict conditions
% d - 'd' matrix (m x m x 0), each m row depicts ROIs, currently represents linear DCM

%% Example inputs
% nROIs = 4;

% nConds = 3;

% intmatrix = [2 2 0 0; 2 2 2 2; 0 1 2 0; 0 1 0 2]; 
% Matrix should be of the form of nROIs x nROIs
% a value of 0 represents connections that can't exist
% a value of 1 represent the Condition where a connection cannot exist if
% forward connection is not present
% a value of 2 represents connections can exist in any combination.

% modmatrix = [1 1 0 0; 1 1 1 1;0 1 1 0; 0 1 0 1];
% A value of 1 means modulation of that connection
% The matrix has to correspond to the intrinsic matrix format

% inpmatix = [0 0 1;0 0 0;0 0 0; 0 0 0];
% has to be of the form nROIs x ncond 
% A value of 1 means a direct connection in the ROI (here 1st ROI and 3rd
% connection) 


%%
%job

nROIs = job.nROIs;
dest = char(job.dest);
intmatrix = job.intmatrix;
modmatrix = job.modmatrix;
inpmatrix = job.inpmatrix;
save_models = job.save;

%% Change to destination folder
try 
    cd(dest);
catch
    print('Error opening folder');
end

%% Make all possible 'intrinsic connectivity - a' matrices for n ROIs
% This program is currently set for 4 ROIs - 1Fusi, 2LIPS, 3LIFG, 4RIPS

x = npermutek([1 0], nROIs*nROIs);
n = size(x,1);
y = cell(1,n);

for i = 1:n
    y{i} =(reshape(x(i,:),nROIs,nROIs))'; %To Do - Make this code generic according to no. of ROIs
end

% Apply intrinsic connectivity conditions to the martices
% Select connections that cannot exist
[p q]= find(intmatrix==0);
% Select connections that cannot exist if no forward connection is present
[u w]= find(intmatrix==1);
m = size(y,2);
z = cell(1,m);
for i = 1:m
    matrix = y{1,i};
    s = size(matrix);
    for j = 1:s
        for k = 1:s
            % Intrinsic connections between same regions
            if j == k
                matrix(j,k) = 1;
            else
            % Intrinsic connections that do not exist
            for t = 1:size(p)
                if k==q(t) & j==p(t)
                    matrix(j,k)= 0;              
                end
            end           
            end
            % Condition where a connection cannot exist if forward
            % connection is not present
            for v = 1:size(u)
                if k==w(v) & j==u(v) & j>k 
                    if matrix(j,k) ==0
                        matrix(k,j) = 0;
                    end
                end
            end
        end
    end
    % Note - Enter any additional conditions for which the model shouldn't exist at all here
    %Example- 
    % Cond - If connection from 1 to 2 is zero, then this matrix should not exist
    % May not be true in all cases (case of first region feeding to 2 distinct regions)
     if matrix (2,1) == 0
        matrix = [];
     end
     z{1,i} = matrix;
end

clear i j k m matrix n s x y p q u v w t;

%% Filter only unique models from the 'a' matrix

% Script to find numbers and index of unique matrices in a cell of
% matrices, taken from internet

% convert to strings
mcs = cellfun(@(x)(mat2str(x)),z,'uniformoutput',false);

% run unique
[~,idxOfUnique,~] = unique(mcs);
n = size(idxOfUnique,2);
a = cell(1,n-1);
for i = 1:n
    if ~ isempty(z{1,idxOfUnique(i)})
        a{1,i} = z{1,idxOfUnique(i)};
    end
end

clear i mcs n idxOfUnique z;

%% Make all possible 'Modulatory effects - b' matrix for each of these intrinsic connectivity models

x = npermutek([1 0], nROIs*nROIs);
n = size(x,1);
y = cell(1,n);
for i = 1:n
    y{i} = (reshape(x(i,:),nROIs,nROIs))'; %To Do - Make this code generic according to no. of ROIs
end

% Apply modulatory effects conditions to the matrices
[p q]= find(modmatrix==0);
m = size(y,2);
z = cell(1,m);
for i = 1:m
    matrix = y{1,i};
    s = size(matrix);
    for j = 1:s
        for k = 1:s
            % Cond1 - All modulation on self connections absent, all diagonals zero
            if j == k
                matrix(j,k) = 0;
            else
            % Cond2 - Modulatory inputs
            for t = 1:size(p)
                if k==q(t) & j==p(t)
                    matrix(j,k)= 0;              
                end
            end
            end

        end
    end
    z{1,i} = matrix;
end

clear i j k m matrix n s x y p q;

%% Filter only unique models from the 'a for b' matrix

% Script to find numbers and index of unique matrices in a cell of
% matrices, taken from internet

% convert to strings
mcs = cellfun(@(x)(mat2str(x)),z,'uniformoutput',false);

% run unique
[~,idxOfUnique,~] = unique(mcs);
n = size(idxOfUnique,2);
afb = cell(1,n-1);
for i = 1:n
    if ~ isempty(z{1,idxOfUnique(i)})
        afb{1,i} = z{1,idxOfUnique(i)};
    end
end

clear i mcs n idxOfUnique z;

%% Make all the models or in other words 'b' matrices by scalar multiplication of 'a' and 'afb'. 
%The idea is to create all combinations of possible models for each intrinsic connectivity matrix

n = size(a,2);
s = size(afb,2);
x = cell(n,s); % Put x = a' if you want to see each row end by it's corresponding intrinsic matrix in 'b' cell
for i = 1:n
    v = a{1,i};
    for j = 1:s
        t = afb{1,j};
        x{i,j} = v.*t;
    end
end

clear afb n s i v j t;

%% Filter only unique models from the 'b' matrix

s = size(x,1);
y = cell(s,1);
for i = 1:s
    [y{i},~,~] = uniquecell(x(i,:));
end

v = size(y,1);
b = cell(v);

for j = 1:s
    t = size(y{j,1},2);
    for k = 1:t
        if sum(sum(y{j,1}{1,k}))
            b{j,k} = y{j,1}{1,k};
        end
    end
end

clear s i x y v j k t;

%% Make the direct input matrix based on user defined condition(s). Note that this matrix doesn't change across models

c = inpmatrix;
%c = zeros(nROIs,nConds);
%p = find(inpmatrix==1);
% Activate(by removing %) or deactivate (by adding %) the code below as per your conditions

% Cond1 - If first region receieves input
%c(1,nConds) = 1;

% Cond2 - If second region receieves input
%c(2,nConds) = 1;

% Condn - If nth regions recieves input
%c(nROIs,nConds) = 1;

%% Create DCM structures containing one intrinsic matrix 'a', all modulatory effects matrices 'b' for this 'a' & one 'c' and 'd' matrices.

% Store 'c' and 'd' matrices in DCM structure
%--------------------------------------------%
DCM.c = c;
DCM.d = zeros(nROIs,nROIs,0); %This program is for linear DCM, change d to implement non-linear DCM

% Store 'a' and 'b' matrices in DCM structure
%--------------------------------------------%
s = size(a,2);

for i = 1:s % For each intinsic connectivity group
    DCM.a = a{1,i}; % save one a matrix
    DCM.b = b(i,:); % save one cell of b matrices
    DCM_model(i)=DCM;
    if exist('save_models')
    if spm_check_version('matlab','7') >= 0
        save(['MODEL_' num2str(i) '_' '.mat'],'-V6','DCM');
    else
        save(['MODEL_' num2str(i) '_' '.mat'],'DCM');
    end
    end
end





