% script to apply output controllability to the UK-Biobank connectomes

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
% batch = '00';
% bilateralNetworks = 0;
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


%% load ROIs Info

if bilateralNetworks
    ROI_info_Table = readtable( 'controllability_ukb/Schaefer200_allinfo_bilateral_networks.csv');
    networksAssignment = table2array(ROI_info_Table(:,5));
    nNetworks = max(networksAssignment);
    save_dir = '9nets_bilateral';
else
    ROI_info_Table = readtable( 'controllability_ukb/Schaefer200_allinfo.csv');
    networksAssignment = table2array(ROI_info_Table(:,5));
    nNetworks = max(networksAssignment);
    save_dir = '9nets';
end

nSystems = nNetworks;


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



%% define arrays

relativeErrorArraySpectralCell = cell(1 , nSystems);
relativeErrorArraySpectral2sysCell = cell(1 , nSystems);
energyArraySpectralCell = cell(1 , nSystems);
lambdaMinGramSpecCell = cell(1 , nSystems);

%%%%%%%%%%%%%%%%%%%%%%
relativeErrorArrayOriginal = zeros(n , nSystems);
energyArrayOriginal = zeros(n , nSystems);
lambdaMinGramAll = zeros(n , nSystems);
lambdaAvgGramAll = zeros(n , nSystems);

%% loop over systems

controllableSize = 7;
orderedModesIdx = 1 : n;

warning off

%%%% deal with Parrallel Obj
patchJobStorageLocation
delete(gcp('nocreate'));
poolObj = parpool('local' ,nCpuOnNode);
parfor kSys = 1 : nSystems

% for kSys = 1: nSystems

    warning off
    fprintf('system %i / %i \n' , kSys , nSystems);

    %% define nodes of the system

    target = find(networksAssignment == kSys);
    targetSize = length(target);

    yf = ones(targetSize , 1) + sigSpred *(randn(targetSize , 1));

    %%%%% big C matrix for target Control %%%%%%%
    bigC = zeros(targetSize , n);
    for ktarget = 1 : targetSize
        bigC(ktarget , target(ktarget)) = 1;
    end

    %% define arrays

    lambdaMinGramSpec = zeros(n , controllableSize);
    relativeErrorArraySpectral = zeros(n , controllableSize);
%     relativeErrorArraySpectralNotNorm = zeros(n , controllableSize);
%     relativeErrorArraySpectralCosine = zeros(n , controllableSize);
    energyArraySpectral = zeros(n , controllableSize);

    %%%%%%%%%%%%%
    relativeErrorArraySpectral2sys = zeros(n , controllableSize);
%     relativeErrorArraySpectral2sysCosine = zeros(n , controllableSize);


    %% define Slepians or Laplacian
    targetNet = matrix(target , target)  ;
    targetlaplac = diag(sum(targetNet)) - targetNet;

    [VnotSorted,lambdaMatNotSorted] = eig( targetlaplac );
    lambdaNotSorted = diag(lambdaMatNotSorted);

    [lambda , idxAsc] = sort(lambdaNotSorted , 'ascend'); % for Laplacian
    V = VnotSorted(: , idxAsc);


    yfGFT = (V')* yf;
    [~ , orderedModesIdx] = sort((abs(yfGFT) ), 'descend');

    %% loop over drivers

    for kDriver = 1:n

        if kDriver == floor(n/2)
            fprintf('system %i : 50%% \n' , kSys);
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
        targetGram = (bigC * W *(bigC'));


        lambdaAll = eig(targetGram);
        [~ , ismallabs] = min(abs(lambdaAll));
        lambdaMinGramAll(kDriver , kSys) = lambdaAll(ismallabs);
        lambdaAvgGramAll(kDriver , kSys) = mean(lambdaAll);

        %% control all target
        [ySimuAll, xSimuAll, Uall , P_adj , n_err , condPinvSpec , condAtilde] = ...
            myControlFun_betterThanStiso_output(A, B , T,STEP, t, x0, zeros(n,1), rhoOpt, ...
            zeros(n,n) , eye(n) , zeros(n,n) , bigC , yf  );

        % Trajectory
        enNumAll = controlEnergy(Uall , STEP);

        relativeErrorArrayOriginal(kDriver , kSys) = norm(ySimuAll(end , :)' - yf , 2)/norm(yf , 2);

        energyArrayOriginal(kDriver , kSys) = enNumAll;


        %% loop over output dimensions
        for kOutDim = 1 : controllableSize
            if kOutDim <= targetSize
%                 fprintf('Subtarget %i / %i \n' , kOutDim , targetSize);

                %% define C matrix
                Cbar = V(:,orderedModesIdx(1 : kOutDim))';
                zf = Cbar * yf;

                %% Gramian of low dimension
                warning off
                Wbar = Cbar * targetGram * Cbar';

                targetGramSpec = (Wbar + Wbar')/2;
                lambdaSpecRed = eig(targetGramSpec);
                [~ , ismallabs] = min(abs(lambdaSpecRed));

                lambdaMinGramSpec(kDriver , kOutDim) = lambdaSpecRed(ismallabs);
                warning on

                %% control low dimension

                [yBarSimu, xSimuext, Ucan , P_adj , n_err , condPinvSpec , condAtilde] = ...
                    myControlFun_betterThanStiso_output(A, B , T,STEP, t, x0, zeros(n,1) , ...
                    rhoOpt,  zeros(n,n) , eye(n) , zeros(n,n) , Cbar*bigC , zf  );

                % Trajectory
                enNumEigen = controlEnergy(Ucan , STEP);
                % %%%%%  Store Results  Eigen  %%%%%%

                relativeErrorArraySpectral(kDriver , kOutDim) = norm(yBarSimu(end , :)' - zf , 2)/norm(zf , 2);
%                 relativeErrorArraySpectralCosine(kDriver , kTargetDim) = cosineDist(yBarSimu(end,:)' , zf);

                %         localityArraySpectral(kTask , kScale) = localRateEigen;
                energyArraySpectral(kDriver , kOutDim) = enNumEigen;
                %%%%%%%%%%%%%
                relativeErrorArraySpectral2sys(kDriver , kOutDim) = norm(xSimuext(end,target)' - yf , 2)/norm(yf , 2);
%                 relativeErrorArraySpectral2sysCosine(kDriver , kTargetDim) = cosineDist(xSimuext(end,target)' , yf);


            end
        end
        % delete(gcp('nocreate'));


        % delete(gcp('nocreate'));

    end
    %% store results

    %
    %         lambdaMinGramAvgCell{1,kSys} = lambdaMinGramAvg;
    lambdaMinGramSpecCell{1,kSys} = lambdaMinGramSpec;
    relativeErrorArraySpectralCell{1,kSys} = relativeErrorArraySpectral;
    relativeErrorArraySpectral2sysCell{1,kSys} =relativeErrorArraySpectral2sys;
    energyArraySpectralCell{1,kSys} = energyArraySpectral;

end
delete(gcp('nocreate'));
fprintf('computation done\n');
warning on

%% save  workspace

save(fullfile('controllability_ukb' , 'ClusterWorkSpaceSave', save_dir , ...
    sprintf("%s_braiApp_gramAndControl.mat" , ...
     subName  )) , ...
     'lambdaMinGramSpecCell' , 'relativeErrorArraySpectralCell' , 'relativeErrorArraySpectral2sysCell', ...
     'energyArraySpectralCell' , 'relativeErrorArrayOriginal' , 'energyArrayOriginal' , 'lambdaMinGramAll' , 'lambdaAvgGramAll');

fprintf('WorkSpace saved ...\n');

toc
