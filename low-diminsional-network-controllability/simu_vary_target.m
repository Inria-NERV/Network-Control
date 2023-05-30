% script to run the third simulation of the paper : varying target size
% with one single driver
% 
clear all;
delete(gcp('nocreate'));

fprintf('Launch Script for Simulation 1\n');
allSubDirectories = genpath('utils');
addpath(allSubDirectories);

path = pwd


%% get environment variable
jobId = getenv('SLURM_ARRAY_JOB_ID');
runId = getenv('SLURM_ARRAY_TASK_ID');
fprintf('arrayID = %s \n' , runId);
nCpuOnNode= str2double(getenv('SLURM_CPUS_ON_NODE'));

rhoExp = str2double( getenv('RHO_EXP')); % regularization parameter
gramControlVar = str2double( getenv('GRAM_CONTROL')); % control with Gramian inversion set to 0


% runId  = '1';
% nCpuOnNode = 2;
% gramControlVar = 0;
% rhoExp =  4;
% plotTraj = 0;

%%%% set the Random Numbers seed %%%%%%%%
rng(str2double(runId))
nCpus = nCpuOnNode;


%% Options Spec

normOpt = 1; dampFact = 1 ; c = 0.1 ; %10^(-3);


%% Network Size %%%%%%%%%%%%%%
newNormalization = 0;
multiForPropag = 2.3;
% multiForPropag = 20;

ordLapXf = 1;
gramControl  = gramControlVar;
directedNet = 0;
rhoOpt = 10^(-rhoExp);

%%%%%% test hierarchical Big
n = 256;
mx_lvl = 8;

%%%%%% NdArray = [1 12 30];
NdArray = [1];
targetSizeArray =  [32 64 96 128 160 192 224 256];

outputDimArray = [1 2 3 4 5 6];


%% %------------------ Time constants ------------------------------
T = 1; %in s
STEP = 0.01 ; %in s
t = 0:STEP:T;
nSamps = length(t);



%% fix initial and final states

sigSpred = 0.1;
x0 = zeros(n , 1) + sigSpred *(randn(n , 1));
xf = ones(n,1) + sigSpred *(randn(n , 1));
% xf = zeros(n,1) + sigSpred *(rand(n , 1)-0.5);


%% define arrays


nTasks =length(outputDimArray);
nOutputs = length(targetSizeArray);

relativeErrorArraySpectral = zeros(nTasks,nOutputs);
energyArraySpectral = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%
relativeErrorArraySpectral2Original = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%
relativeErrorArrayAvg = zeros(nTasks,nOutputs);
energyArrayAvg = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%
relativeErrorArrayAvg2Original = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%
errorAllTargetsOrigin = zeros(nTasks,nOutputs);
energyAllTargetsOrigin = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
relativeErrorArrayAvgReduc = zeros(nTasks,nOutputs);
energyArrayAvgReduc = zeros(nTasks,nOutputs);
relativeErrorArrayAvgReduc2orig = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%
relativeErrorArrayLapReduc = zeros(nTasks,nOutputs);
energyArrayLapReduc = zeros(nTasks,nOutputs);
relativeErrorArrayLapReduc2orig = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%
%%%%%  not normalised errors.%%%%%%
relativeErrorArraySpectralNotNorm = zeros(nTasks,nOutputs);
energyArraySpectralNotNorm = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%
relativeErrorArraySpectral2OriginalNotNorm = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%
relativeErrorArrayAvgNotNorm = zeros(nTasks,nOutputs);
energyArrayAvgNotNorm = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%
relativeErrorArrayAvg2OriginalNotNorm = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%
errorAllTargetsOriginNotNorm = zeros(nTasks,nOutputs);
energyAllTargetsOriginNotNorm = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
relativeErrorArrayAvgReducNotNorm = zeros(nTasks,nOutputs);
energyArrayAvgReducNotNorm = zeros(nTasks,nOutputs);
relativeErrorArrayAvgReduc2origNotNorm = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%
relativeErrorArrayLapReducNotNorm = zeros(nTasks,nOutputs);
energyArrayLapReducNotNorm = zeros(nTasks,nOutputs);
relativeErrorArrayLapReduc2origNotNorm = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%

%%%%% Cosine distances %%%%%%
% relativeErrorArraySpectral = zeros(nTasks,nOutputs);
energyArraySpectralCosine = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%
relativeErrorArraySpectral2OriginalCosine = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%
relativeErrorArrayAvgCosine = zeros(nTasks,nOutputs);
energyArrayAvgCosine = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%
relativeErrorArrayAvg2OriginalCosine = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%
errorAllTargetsOriginCosine = zeros(nTasks,nOutputs);
energyAllTargetsOriginCosine = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
relativeErrorArrayAvgReducCosine = zeros(nTasks,nOutputs);
energyArrayAvgReducCosine = zeros(nTasks,nOutputs);
relativeErrorArrayAvgReduc2origCosine = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%
relativeErrorArrayLapReducCosine = zeros(nTasks,nOutputs);
energyArrayLapReducCosine = zeros(nTasks,nOutputs);
relativeErrorArrayLapReduc2origCosine = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%

errorAllTargetsOriginNotCosine = zeros(nTasks,nOutputs);


graphs = cell(nTasks,nOutputs);
communities = cell(1,nOutputs);
eigenmaps = cell(1,nOutputs);
bestDrivers = cell(1,nOutputs);


modularityArray = zeros(1,nOutputs);
densityArray = zeros(1,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%

spectralObsRankArray = zeros(1,nOutputs);
averageObsRankArray = zeros(1,nOutputs);

spectralObsCondArray = zeros(1,nOutputs);
averageObsCondArray = zeros(1,nOutputs);

diffReconstructedOriginalPinvAvg = zeros(1,nOutputs);
diffReconstructedOriginalPinvSpec = zeros(1,nOutputs);

diffReconstructedOriginalPinvAvgCosine = zeros(1,nOutputs);
diffReconstructedOriginalPinvSpecCosine = zeros(1,nOutputs);


lambdaMinGramAvg = zeros(nTasks,nOutputs);
lambdaMinGramSpec = zeros(nTasks,nOutputs);
lambdaMeanGramAvg = zeros(nTasks,nOutputs);
lambdaMeanGramSpec = zeros(nTasks,nOutputs);

lambdaMinGramAll = zeros(nTasks,nOutputs);

lambdaMinGramAvgRed = zeros(nTasks,nOutputs);
lambdaMinGramSpecRed = zeros(nTasks,nOutputs);

lambdaMeanGramAvgRed = zeros(nTasks,nOutputs);
lambdaMeanGramSpecRed = zeros(nTasks,nOutputs);


condHamiltonAvg = zeros(nTasks,nOutputs);
condHamiltonSpec = zeros(nTasks,nOutputs);
condHamiltonAvgRed = zeros(nTasks,nOutputs);
condHamiltonSpecRed = zeros(nTasks,nOutputs);


%% define modular graphs parameters


Diag = 0;
%%%%
rho = 0.05; % rho = 0.2;

L = floor(rho*(n-1)*n/2); % for undirected
% L = floor(rho*(nOrig-1)*nOrig); % for directed

Karray = floor(targetSizeArray);
N1array = floor(1./Karray * n);


R = eye(n,n);
Q = zeros(n,n);

%% generate modular graph
kScale =1;

E = 2.5;

sz_cl = 2;


isConnected = 0;
% make sure the graph is connected
kInf = 0;
while ~isConnected
    %             [matrix ,V0] = GGGirvanNewman(N1,K,zi,ze,Diag); % old way
    kInf = kInf + 1;
    if mod(kInf , 100) == 0
        fprintf('100 attempts of generating modular graph in scale %i \n', kScale);
    end

    [matrix,K] = makefractalCIJ(mx_lvl,E,sz_cl);
    matrix = triu(matrix );
    matrix = matrix + matrix';

    matrix = matrix - diag(diag(matrix));


    G = graph(matrix) ;
    connBins = conncomp(G);
    differentComp = find(connBins ~= 1 , 1);
    isConnected = isempty(differentComp);
end

%%%%%%%%%%%%%
densityArray(kScale) = sum(sum(matrix))/(n-1)/n;

graphs{ 1, kScale} = sparse(matrix);




%% normalize matrix A and relax


lambdaMax = eigs(matrix , 1);
A = multiForPropag * ((matrix / (1.01 * lambdaMax)) - eye(n));


%% loop over Community outputs
% orderedModesIdx = 1 : n; %  low to high
dimNormLow = 2;
dimNormHigh = 2;
dimNorm = 2;
nDriversMax = max(NdArray);

warning off

%%%% deal with Parrallel Obj
% delete(gcp('nocreate'));
% patchJobStorageLocation
% poolObj = parpool('local' ,nCpus);
% % poolObj = parpool('SlurmProfile1' ,nCpus);

% parfor kScale = 1 : nOutputs
for kScale = 1 : nOutputs

    warning off
    %% new target
    targetSize = targetSizeArray(kScale);
    target = 1:targetSize;
    yf = (xf(target));

    %%%%%%%%%%%%%%%%%%%%%%
    %%%%% C matrix for target Control %%%%%%%
    C = zeros(targetSize , n);
    for ktarget = 1 : targetSize
        C(ktarget , target(ktarget)) = 1;
    end

    %%  define eigenModes of the target
    %     V = gsp_laplacian_eigenmaps(Ggsp, nOrig); %% not normalised laplacian

    Ggsp = gsp_graph(matrix(target, target));

    Ggsp = gsp_create_laplacian(Ggsp, 'combinatorial');

    undLap = full(Ggsp.L);
    targetlaplac = undLap;


    [VnotSorted,lambdaMatNotSorted] = eig( targetlaplac );
    lambdaNotSorted = diag(lambdaMatNotSorted);

    [lambda , idxAsc] = sort(lambdaNotSorted , 'ascend'); % for Laplacian
    V = VnotSorted(: , idxAsc);


    %  get GFT of xf;
    if ordLapXf
        xfGFT = (V')*(yf );
        x0GFT = (V')*(C*x0);
        %     [lambdaMag , idxGFTxf] = sort(abs(abs(xfGFT) - abs(x0GFT))./(abs(xfGFT) + abs(x0GFT)) , 'descend');
        [lambdaMag , idxGFTxf] = sort((abs(xfGFT) ), 'descend');
        orderedModesIdx =idxGFTxf; % descend
    else
        orderedModesIdx =1:n;
    end


    VnotSorted = nan; targetlaplac = nan; Ggsp = nan;

    eigenmaps{1,kScale} = V;

    %% new driver place and matrix B

    %% rank drivers acording to betweenness centrality
    %     BC =betweenness_wei(1./matrix);
    BC =betweenness_bin(matrix(target , target )) ;
    [driversBC , idxBestDrivers] = sort(BC , 'descend');
    %     [driversBC , idxBestDrivers] = sort(BC , 'descend'  );

    bestDrivers{1,kScale} = idxBestDrivers;

    % idxBestDrivers = (1:n);


    %% build matrix B

    nDrivers = 1;
    '''% %%%%%% Drivers %%%%%%%%%''';
    driversInd = idxBestDrivers(1: nDrivers);

    %%%%%%%%%%%%%%%%%
    Drivers = zeros(1,n);
    Drivers(driversInd) = 1;
    % B = zeros(n);
    if ~relaxB
        B = zeros(n);
    else
        B = sigma * randn(n) + mu * ones(n);
    end
    B = B - diag(diag(B)) + diag(Drivers);

    %% global Gramian
    [W,Tr] = CtrGram(A,B);

    %% Get Gramian indices
    targetGram = C*W*C';

    lambdaAll = eig(targetGram);
    lambdaMinGramAll(1 , kScale) = min(lambdaAll);

    %% control all targets
    if gramControl
        [ySimu, xSimu, U ] = myOutputControlFunction_GramInv(A, B, T,STEP, t, x0, yf , C);
    else
        [ySimu, xSimu, U , P_adj ,n_err , condPinv , condAtilde] = ...
            myOutputControlFunction(A, B, T, t, x0, xf, rhoOpt , R , Q , C , yf  );
    end

    % Trajectory
    %             enNumSpan = controlEnergy(U , STEP);
    errorAllTargetsOrigin(1 , kScale) = norm(ySimu(end,:)' - yf , dimNorm)/norm(yf, dimNorm);
    errorAllTargetsOriginNotNorm(1 , kScale) = norm(ySimu(end,:)' - yf , dimNorm);
    errorAllTargetsOriginNotCosine(1 , kScale) = cosineDist(ySimu(end,:)', yf);
    energyAllTargetsOrigin(1 , kScale) = controlEnergy(U , STEP);


    %% loop over Tasks
    for kTask = 1 : nTasks
        fprintf('\n scale: %i : task %i/%i ;\n' , kScale , kTask , nTasks);

        warning off
        outputDim = outputDimArray(kTask);
        K = outputDim;
        N1 = targetSize/K;

        %% deal with communities

        myCommunities = zeros(1 , targetSize);
        for kCom = 1 : K
            if kCom < K
                myCommunities(N1*(kCom - 1)+1 : N1*(kCom - 1)+ N1) = kCom;
            else
                myCommunities(N1*(kCom - 1)+1 : end) = kCom;
            end
        end
        %%%%%%%%%%


        %%%%% C matrix for Communities %%%%%%%
        CavgTarget = zeros(outputDim , targetSize);
        for kCom = 1 : outputDim
            nodesInCom = find(myCommunities == kCom )';
            CavgTarget(kCom , nodesInCom) = 1/sqrt(length(nodesInCom));

        end
        %     Cavg = Cavg./(sum(Cavg, 2)*ones(1,n));
        %%%%%%%%%%%%%%%%%%%%%
        communities{kScale} = myCommunities;

        Cavg = CavgTarget*C;
        yfAvg = Cavg*xf;

        %% deal with eigenmaps of the target


        modes2take = orderedModesIdx(1 : outputDim);
        Cbar = (V( : ,modes2take)')* C;
        zf = Cbar *xf;


        %% Get Gramian indices

        avgGram = (Cavg * W *(Cavg'));
        specGram = (Cbar * W *(Cbar'));

        lambdaAvg = eig(avgGram);
        lambdaSpec = eig(specGram);


        lambdaMinGramAvg(kTask , kScale) = min(lambdaAvg);
        lambdaMinGramSpec(kTask , kScale) = min(lambdaSpec);

        lambdaMeanGramAvg(kTask , kScale) = sum(lambdaAvg);
        lambdaMeanGramSpec(kTask , kScale) = sum(lambdaSpec);



        %% Control  average
        if gramControl
            [ySimuAvg, xSimuAverage, Uavg ] = myOutputControlFunction_GramInv(A, B, T, t, x0, yfAvg , Cavg);
        else
            [ySimuAvg, xSimuAverage, Uavg , P_adj ,n_err , condPinvAvg , condAtilde] = ...
                myOutputControlFunction(A, B, T, t, x0, xf, rhoOpt , R , Q , Cavg , yfAvg  );
        end

        % Trajectory
        enNumAvg = controlEnergy(Uavg , STEP);

        % Store results Average
        %%%%%%%%%%%%%%%%%%%%%%
        relativeErrorArrayAvg(kTask , kScale) = norm(ySimuAvg(end ,:)' - yfAvg , dimNormLow)/norm(yfAvg , dimNormLow);
        relativeErrorArrayAvgNotNorm(kTask , kScale) = norm(ySimuAvg(end ,:)' - yfAvg , dimNormLow);
        relativeErrorArrayAvgCosine(kTask , kScale) = cosineDist(ySimuAvg(end,:)' , yfAvg);

        %         localityArrayAvg = localRateAvg;
        energyArrayAvg(kTask , kScale) = enNumAvg;
        %%%%%%%%%%%%%
        relativeErrorArrayAvg2Original(kTask , kScale) = norm(xSimuAverage(end,target)' - yf , dimNorm)/norm(yf , dimNorm);
        relativeErrorArrayAvg2OriginalNotNorm(kTask , kScale) = norm(xSimuAverage(end,target)' - yf , dimNorm);
        relativeErrorArrayAvg2OriginalCosine(kTask , kScale) = cosineDist(xSimuAverage(end,target)' , yf);
        %         condHamiltonAvg(kTask , kScale) = condPinvAvg;

        %%%%%%%%%%%%%%%%%%%%%%
        """ -------------------------- """;
        """ move on to the Eigen Space """;
        """ -------------------------- """;

        %% Control eigen-modes
        if gramControl
            [yBarSimu, xSimuext, Ucan ] = myOutputControlFunction_GramInv(A, B, T, t, x0, zf , Cbar);
        else
            [yBarSimu, xSimuext, Ucan , P_adj , n_err , condPinvSpec , condAtilde] = ...
                myOutputControlFunction(A, B , T, t, x0, xf, rhoOpt, R , Q , Cbar , zf  );
        end
        % Trajectory
        enNumEigen = controlEnergy(Ucan , STEP);


        % %%%%%  Store Results  Eigen  %%%%%%

        relativeErrorArraySpectral(kTask , kScale) = norm(yBarSimu(end , :)' - zf , dimNormLow)/norm(zf , dimNormLow);
        relativeErrorArraySpectralNotNorm(kTask , kScale) = norm(yBarSimu(end , :)' - zf , dimNormLow);
        relativeErrorArraySpectralCosine(kTask , kScale) = cosineDist(yBarSimu(end,:)' , zf);

        %         localityArraySpectral(kTask , kScale) = localRateEigen;
        energyArraySpectral(kTask , kScale) = enNumEigen;
        %%%%%%%%%%%%%
        relativeErrorArraySpectral2Original(kTask , kScale) = norm(xSimuext(end,target)' - yf , dimNorm)/norm(yf , dimNorm);
        relativeErrorArraySpectral2OriginalNotNorm(kTask , kScale) = norm(xSimuext(end,target)' - yf , dimNorm);

        relativeErrorArraySpectral2OriginalCosine(kTask , kScale) = cosineDist(xSimuext(end,target)' , yf);

        %         condHamiltonSpec(kTask , kScale) = condPinvSpec;


    end
    warning on
end
warning on
delete(gcp('nocreate'));

fprintf('computation done ...\n')

%% save workSpace

if gramControl
    save(sprintf(...
        "ClusterWorkSpaceSave/single_driver_vary_target/varyTarget_driver_central_changing_und_hierarchical_strengthFactor_%i_ordLapXf_gram_Control_centralDrivers_Nd_%i_%i_%i_outputs_%i_n_%i_rho_%i_percent_Xf_02_runID_%s.mat" , ...
        multiForPropag , NdArray(1) , NdArray(2) , NdArray(3)  , nOutputs , n , floor(100 * rho) , runId));
else
    save(sprintf(...
        "ClusterWorkSpaceSave/single_driver_vary_target/varyTarget_driver_central_changing_und_hierarchical_strengthFactor_%iordLapXf_rhoOpt_%.2e_Control_centralDrivers_out_%i_%i_%i_outputs_%i_n_%i_rho_%i_percent_Xf_02_runID_%s.mat" , ...
        multiForPropag , rhoOpt , outputDimArray(1) , outputDimArray(2) , outputDimArray(3)  , nOutputs , n , floor(100 * rho) , runId));
end

fprintf('WorkSpace saved ...\n')
