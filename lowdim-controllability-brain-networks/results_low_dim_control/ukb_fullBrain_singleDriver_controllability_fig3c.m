% script to study (full network) low-dimensional controllability to the UK-Biobank connectomes

clear all;

tic
%% Batch data and subject

%%%%%% get environment variable
jobId = getenv('SLURM_ARRAY_JOB_ID');
nCpuOnNode= str2double(getenv('SLURM_CPUS_ON_NODE'));

subId = str2double(getenv('SLURM_ARRAY_TASK_ID')) ;
fprintf('arrayID = %i \n' , subId);

bilateralNetworks = str2double(getenv('BILATERALNETS'));

% wd = '/network/lustre/iss02/ukbiobank/ukbiobank-remaining/remy.benmessaoud/';
% data_dir = 'allConnectomes_healthy';


%%%%%%%%%%%
%%%%%%% local run

% nCpuOnNode = 3;
subId = 1;
% wd = '/Users/remy.benmessaoud/ownCloud/remy.benmessaoud/MATLAB/UK-biobank/';
data_dir = 'allConnectomes_healthy';


% cd(wd)
% cd the results directory

%% %%%% add dependencies
allSubDirectories = genpath('utils');
addpath(allSubDirectories);

%% get list of subjects in batch
% batch_dir = [data_dir 'batch' batch];

list_subjects = dir(data_dir);
%%%% get index of first subject folder
idxGood = 0;
notGood = 1;
while notGood
    idxGood = idxGood+1;
    nameFolder = list_subjects(idxGood).name;
    if ~strcmp(nameFolder(1) , '.')
        notGood = 0;
    end
end


%%%%%%%%%
list_subjects = list_subjects(idxGood :end);
nSubjects = length(list_subjects) ;

parc_name = 'schaefer200';
schaeferSize = 216;

subName = list_subjects(subId).name

% go to data of subject and read table
subject_dir = fullfile(data_dir , subName );

tableName = sprintf('%s_connectome_edgelist_schaefer200-yeo17_rmap_s1dild.csv' , subName ) ;
tablePath = fullfile(subject_dir , tableName);

matrix = load_UKBiobank_connectome(tablePath , schaeferSize);


n = size(matrix, 1);


%% normalize matrix to A

lambdaMax = eigs(matrix , 1);
A =  matrix - * (lambdaMax + eps) * eye(n);
% the eps is to ensure Re(Lambda(A))<0 to be stable

warning off

%% define nodes of the system

target = 1:n;
targetSize = length(target);

%%%%% big C matrix for target Control %%%%%%%
bigC = zeros(targetSize , n);
for ktarget = 1 : targetSize
    bigC(ktarget , target(ktarget)) = 1;
end


%% define  Laplacian

targetlaplac = diag(sum(matrix)) - matrix;

[VnotSorted,lambdaMatNotSorted] = eig( targetlaplac );
lambdaNotSorted = diag(lambdaMatNotSorted);

[lambda , idxAsc] = sort(lambdaNotSorted , 'ascend'); % for Laplacian
V = VnotSorted(: , idxAsc);


%% define arrays
controllableSize = 7;
orderedModesIdx = 1 : n;


lambdaMinGramAll = zeros(n , 1);
lambdaAvgGramAll = zeros(n , 1);


lambdaMinGramLap = zeros(n , controllableSize);
lambdaAvgGramLap = zeros(n , controllableSize);


%% loop over drivers

for kDriver = 1:n

    if mod(kDriver,10) == 0
        fprintf('driver %i  \n' , kDriver);
    end

    %% build matrix B
    '''% %%%%%% Drivers %%%%%%%%%''';
    driversInd = kDriver ;
    %%%%%%%%%%%%%%%%%
    Drivers = zeros(1,n);
    Drivers(driversInd) = 1;
    B = zeros(n);
    B = B - diag(diag(B)) + diag(Drivers);

    %% get Gramian metrics

    sys = ss(A,B,eye(n) , zeros(n,n));
    W = gram(sys,'c');

    W = (W + W')/2;

    lambdaAll = eig(W);
    [~ , ismallabs] = min(abs(lambdaAll));
    lambdaMinGramAll(kDriver , 1) = lambdaAll(ismallabs);
    lambdaAvgGramAll(kDriver , 1) = mean(lambdaAll);


    %% loop over output dimensions
    for kOutDim = 1 : controllableSize
        if kOutDim <= targetSize
            %                 fprintf('Subtarget %i / %i \n' , kOutDim , targetSize);

            %% define C matrix
            Cbar = V(:,orderedModesIdx(1 : kOutDim))';

            %% Gramian of low dimension
            warning off
            Wbar = Cbar * W * Cbar';

            targetGramSpec = (Wbar + Wbar')/2;
            lambdaSpecRed = eig(targetGramSpec);
            [~ , ismallabs] = min(abs(lambdaSpecRed));

            lambdaMinGramLap(kDriver , kOutDim) = lambdaSpecRed(ismallabs);
            lambdaAvgGramLap(kDriver , kOutDim) = mean(lambdaSpecRed);

        end
    end
end
%% store results

fprintf('computation done\n');
warning on

%% save  workspace

save(fullfile('controllability_ukb' , 'ClusterWorkSpaceSave', 'fullBrainControl' , ...
    sprintf("%s_fullbrainApp_gramAndControl.mat" , ...
    subName  )) , ...
     'lambdaMinGramAll' , 'lambdaMinGramAll' ,...
    'lambdaAvgGramAll' , 'lambdaMinGramLap', 'lambdaAvgGramLap');

fprintf('WorkSpace saved ...\n');

toc


