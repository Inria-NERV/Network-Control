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

wd = '/network/lustre/iss02/ukbiobank/ukbiobank-remaining/remy.benmessaoud/';
data_dir = 'allConnectomes_healthy';


%%%%%%%%%%%
%%%%%%% local run

% nCpuOnNode = 3;
% subId = 1;
% wd = '/Users/remy.benmessaoud/ownCloud/remy.benmessaoud/MATLAB/UK-biobank/';
% data_dir = 'allConnectomes_healthy';


cd(wd)

%% %%%% add dependencies
allSubDirectories = genpath('controllability_ukb');
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

% matrix = load_UKBiobank_connectome(tableName , schaeferSize);
matrix = load_UKBiobank_connectome(tablePath , schaeferSize);


%% %------------------ Control constants ------------------------------
T = 1; %in s
STEP = 0.01 ; %in s
t = 0:STEP:T;
nSamps = length(t);
%%%%%%%%%
n = size(matrix, 1);
R = eye(n,n);
Q = zeros(n,n);

n = size(matrix,1);

sigSpred = 0.1;
x0 = zeros(n , 1) + sigSpred *(randn(n , 1));

rhoOpt = 10^(-6);
multiForPropag = 2.3;


%% normalize matrix to A

lambdaMax = eigs(matrix , 1);
A = multiForPropag * ((matrix / (1.01 * lambdaMax)) - eye(n));


warning off

%% define nodes of the system

target = 1:n;
targetSize = length(target);

yf = ones(targetSize , 1) + sigSpred *(randn(targetSize , 1));

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


yfGFT = (V')* yf;
[~ , orderedModesIdx] = sort((abs(yfGFT) ), 'descend');


%% define arrays
controllableSize = 7;
orderedModesIdx = 1 : n;

%%%%%%%%%%%%%%%%%%%%%%
relativeErrorArrayOriginal = zeros(n , 1);
energyArrayOriginal = zeros(n , 1);

lambdaMinGramAll = zeros(n , 1);
lambdaAvgGramAll = zeros(n , 1);

relativeErrorArrayLap = zeros(n , controllableSize);
energyArrayLap = zeros(n , controllableSize);

lambdaMinGramLap = zeros(n , controllableSize);
lambdaAvgGramLap = zeros(n , controllableSize);

relativeErrorArraySpectral = zeros(n , controllableSize);
energyArraySpectral = zeros(n , controllableSize);

relativeErrorArraySpectral2sys = zeros(n , controllableSize);
%     relativeErrorArraySpectral2sysCosine = zeros(n , controllableSize);

%% loop over drivers

for kDriver = 1:n

    if mod(kDriver,10) == 0
        fprintf('driver %i  \n' , kDriver);
    end
    warning off

    %% build matrix B
    '''% %%%%%% Drivers %%%%%%%%%''';
    driversInd = kDriver ;
    %%%%%%%%%%%%%%%%%
    Drivers = zeros(1,n);
    Drivers(driversInd) = 1;
    B = zeros(n);
    B = B - diag(diag(B)) + diag(Drivers);

    %% get Gramian metrics
    [W,Tr] = CtrGram(A,B);

    warning on

    W = (W + W')/2;

    lambdaAll = eig(W);
    [~ , ismallabs] = min(abs(lambdaAll));
    lambdaMinGramAll(kDriver , 1) = lambdaAll(ismallabs);
    lambdaAvgGramAll(kDriver , 1) = mean(lambdaAll);

    %% control all target
    [ySimuAll, xSimuAll, Uall , P_adj , n_err , condPinvSpec , condAtilde] = ...
        myOutputControlFunction(A, B , T, t, x0, zeros(n,1), rhoOpt, ...
         R , Q , eye(n) , yf  );

    % Trajectory
    enNumAll = controlEnergy(Uall , STEP);

    relativeErrorArrayOriginal(kDriver , 1) = norm(ySimuAll(end , :)' - yf , 2)/norm(yf , 2);

    energyArrayOriginal(kDriver , 1) = enNumAll;


    %% loop over output dimensions
    for kOutDim = 1 : controllableSize
        if kOutDim <= targetSize
            %                 fprintf('Subtarget %i / %i \n' , kOutDim , targetSize);

            %% define C matrix
            Cbar = V(:,orderedModesIdx(1 : kOutDim))';
            zf = Cbar * yf;

            %% Gramian of low dimension
            warning off
            Wbar = Cbar * W * Cbar';

            targetGramSpec = (Wbar + Wbar')/2;
            lambdaSpecRed = eig(targetGramSpec);
            [~ , ismallabs] = min(abs(lambdaSpecRed));

            lambdaMinGramLap(kDriver , kOutDim) = lambdaSpecRed(ismallabs);
            lambdaAvgGramLap(kDriver , kOutDim) = mean(lambdaSpecRed);
            warning on

            %% control low dimension

            [yBarSimu, xSimuext, Ucan , P_adj , n_err , condPinvSpec , condAtilde] = ...
                myOutputControlFunction(A, B , T, t, x0, zeros(n,1) , ...
                rhoOpt , R , Q , Cbar*bigC , zf  );

            % Trajectory
            enNumEigen = controlEnergy(Ucan , STEP);
            % %%%%%  Store Results  Eigen  %%%%%%

            relativeErrorArraySpectral(kDriver , kOutDim) = norm(yBarSimu(end , :)' - zf , 2)/norm(zf , 2);

            energyArraySpectral(kDriver , kOutDim) = enNumEigen;
            %%%%%%%%%%%%%
            relativeErrorArraySpectral2sys(kDriver , kOutDim) = norm(xSimuext(end,target)' - yf , 2)/norm(yf , 2);


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
    'relativeErrorArrayOriginal' , 'energyArrayOriginal' , 'relativeErrorArrayLap', ...
    'relativeErrorArrayLap' , 'energyArrayLap' , 'lambdaMinGramAll' , 'lambdaMinGramAll' ,...
    'lambdaAvgGramAll' , 'lambdaMinGramLap', 'lambdaAvgGramLap', 'relativeErrorArraySpectral'...
    , 'energyArraySpectral', 'relativeErrorArraySpectral2sys');

fprintf('WorkSpace saved ...\n');

toc


