% script to run the third simulation of the paper : varying target size
% with each node taken as single driver :

clear all;

fprintf('Launch Script for Simulation 1\n');
allSubDirectories = genpath('utils');
addpath(allSubDirectories);

didi = pwd


%% get environment variable
jobId = getenv('SLURM_ARRAY_JOB_ID');
runId = getenv('SLURM_ARRAY_TASK_ID');
fprintf('arrayID = %s \n' , runId);
nCpuOnNode= str2double(getenv('SLURM_CPUS_ON_NODE'));
rhoExp = str2double( getenv('RHO_EXP'));


runId  = '1';
% nCpuOnNode = 2;
rhoExp =  4;

%%%% set the Random Numbers seed %%%%%%%%
rng(str2double(runId))
nCpus = nCpuOnNode;


%% Network Options

ordLapXf = 1;
directedNet = 0;
rhoOpt = 10^(-rhoExp);

%%%%%% test hierarchical Big
n = 256;
mx_lvl = 8;

%%%%%% NdArray = [1 12 30];
NdArray = [1];
targetSizeArray =  [32 64 96 128 160 192 224 256];
targetSizeArray =  [4 8 16 32 64 128 256];

outputDimArray = [2 8 32];

outputDimArray = [4];

%% %------------------ Time constants ------------------------------
T = 1; %in s
STEP = 0.01 ; %in s
t = 0:STEP:T;
nSamps = length(t);



%% fix initial and final states

sigSpred = 10;
x0 = zeros(n , 1) + 0 *(randn(n , 1));
xf = ones(n,1) + sigSpred *(randn(n , 1));
% xf = zeros(n,1) + sigSpred *(rand(n , 1)-0.5);


%% define arrays

nTasks =length(outputDimArray);
nOutputs = length(targetSizeArray);

relativeErrorArraySpectral = nan(nTasks,nOutputs);
energyArraySpectral = nan(nTasks,nOutputs);
%%%%%%%%%%%%%
relativeErrorArraySpectral2Original = nan(nTasks,nOutputs);

errorAllTargetsOrigin = nan(nTasks,nOutputs);
energyAllTargetsOrigin = nan(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%

%%%%%  not normalised errors.%%%%%%
relativeErrorArraySpectralNotNorm = nan(nTasks,nOutputs);
energyArraySpectralNotNorm = nan(nTasks,nOutputs);
%%%%%%%%%%%%%
relativeErrorArraySpectral2OriginalNotNorm = nan(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%
relativeErrorArrayAvgNotNorm = nan(nTasks,nOutputs);
energyArrayAvgNotNorm = nan(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%
relativeErrorArrayAvg2OriginalNotNorm = nan(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%
errorAllTargetsOriginNotNorm = nan(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% Cosine distances %%%%%%
% relativeErrorArraySpectral = zeros(nTasks,nOutputs);
relativeErrorArraySpectralCosine = nan(nTasks,nOutputs);

relativeErrorArrayAvg2OriginalCosine = nan(nTasks,nOutputs);


errorAllTargetsOriginCosine = nan(nTasks,nOutputs);


graphs = cell(nTasks,nOutputs);
communities = cell(1,nOutputs);
eigenmaps = cell(1,nOutputs);
bestDrivers = cell(1,nOutputs);



relativeErrorArraySpectralIn =  nan(nTasks,nOutputs);
relativeErrorArraySpectralNotNormIn = nan(nTasks,nOutputs);
relativeErrorArraySpectralCosineIn = nan(nTasks,nOutputs);
energyArraySpectralIn = nan(nTasks,nOutputs);
relativeErrorArraySpectral2OriginalIn  = nan(nTasks,nOutputs);
relativeErrorArraySpectral2OriginalNotNormIn = nan(nTasks,nOutputs);
relativeErrorArraySpectral2OriginalCosineIn = nan(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%

relativeErrorArraySpectralOut =  nan(nTasks,nOutputs);
relativeErrorArraySpectralNotNormOut = nan(nTasks,nOutputs);
relativeErrorArraySpectralCosineOut = nan(nTasks,nOutputs);
energyArraySpectralOut = nan(nTasks,nOutputs);
relativeErrorArraySpectral2OriginalOut  = nan(nTasks,nOutputs);
relativeErrorArraySpectral2OriginalNotNormOut = nan(nTasks,nOutputs);
relativeErrorArraySpectral2OriginalCosineOut = nan(nTasks,nOutputs);



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

    %         [matrixOrig ,V0] = GGGirvanNewman(N1,K ,zi,ze,Diag); % new way to have similar topology from scale to other
    %         Ggsp = gsp_stochastic_block_graph(n, K, params);
    %
    %         matrix = full(Ggsp.W);


    [matrix,K] = makefractalCIJ(mx_lvl,E,sz_cl);
    matrix = triu(matrix );
    matrix = matrix + matrix';

    matrix = matrix - diag(diag(matrix));


    G = graph(matrix) ;
    connBins = conncomp(G);
    differentComp = find(connBins ~= 1 , 1);
    isConnected = isempty(differentComp);
end


densityArray(kScale) = sum(sum(matrix))/(n-1)/n;

%     matrix = [extraNet connectingBloc ; connectingBloc' matrixOrig];
%     matrix = [zeros(size(extraNet)) zeros(size(connectingBloc)) ; connectingBloc' matrixOrig];


graphs{ 1, kScale} = sparse(matrix);

%% rank drivers acording to betweenness centrality
%     BC =betweenness_wei(1./matrix);
BC =betweenness_bin(matrix) ;
[driversBC , idxBestDrivers] = sort(BC , 'descend');
%     [driversBC , idxBestDrivers] = sort(BC , 'descend'  );

bestDrivers{1,kScale} = idxBestDrivers;


%% normalize matrix A and relax


lambdaMax = eigs(matrix , 1);
A =  matrix - 1.001 * lambdaMax* eye(n);
% the coef of 1.001 is to ensure Re(Lambda(A))<0 to be stable


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
    % Ggsp = gsp_create_laplacian(Ggsp, 'normalized');
    %     Ggsp = gsp_create_laplacian(Ggsp, 'chung');

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


    %% control all targets
    % if gramControl
    %     [ySimu, xSimu, U ] = myMinEneregyControlFun_works_new(A, B, T,STEP, t, x0, yf , C);
    % else
    %     [ySimu, xSimu, U , P_adj ,n_err , condPinv , condAtilde] = ...
    %         myControlFun_betterThanStiso_output(A, B, T,STEP, t, x0, xf, rhoOpt,  zeros(n,n) , R , Q , C , yf  );
    % end
    % %%%% try signal on non normalized matrix
    % %         if apply2real
    % %             sys_xp = ss(matrix,B,C ,0);
    % %             [ySimu,tt,xSimuOriginal] = lsim(sys_xp, U,t,x0);
    % %         end
    %
    % % Trajectory
    % %             enNumSpan = controlEnergy(U , STEP);
    % errorAllTargetsOrigin(1 , kScale) = norm(ySimu(end,:)' - yf , dimNorm)/norm(yf, dimNorm);
    % errorAllTargetsOriginNotNorm(1 , kScale) = norm(ySimu(end,:)' - yf , dimNorm);
    % errorAllTargetsOriginNotCosine(1 , kScale) = cosineDist(ySimu(end,:)', yf);
    % energyAllTargetsOrigin(1 , kScale) = controlEnergy(U , STEP);


    %% loop over Tasks
    for kTask = 1 : nTasks
        fprintf('\n scale: %i : task %i/%i ;\n' , kScale , kTask , nTasks);

        warning off
        outputDim = outputDimArray(kTask);

        if outputDim <= targetSize
            %% deal with eigenmaps of the target


            modes2take = orderedModesIdx(1 : outputDim);
            Cbar = (V( : ,modes2take)')* C;
            zf = Cbar *xf;


            %% loop over drivers

            nDrivers2try = 256/16;
            drivers2try = floor(linspace(1 , n , nDrivers2try));


            relativeErrorArraySpectralAllDrivers = zeros(1,nDrivers2try);
            relativeErrorArraySpectralNotNormAllDrivers = zeros(1,nDrivers2try);
            relativeErrorArraySpectralCosineAllDrivers = zeros(1,nDrivers2try);
            energyArraySpectralAllDrivers = zeros(1,nDrivers2try);
            relativeErrorArraySpectral2OriginalAllDrivers = zeros(1,nDrivers2try);
            relativeErrorArraySpectral2OriginalNotNormAllDrivers = zeros(1,nDrivers2try);
            relativeErrorArraySpectral2OriginalCosineAllDrivers = zeros(1,nDrivers2try);

            %%%%%%%%%%%%%%%%%%
            for kDriver = 1 : nDrivers2try
                %% new driver place and matrix B

                %% rank drivers acording to betweenness centrality
                % %     BC =betweenness_wei(1./matrix);
                % BC =betweenness_bin(matrix(target , target )) ;
                % [driversBC , idxBestDrivers] = sort(BC , 'descend');
                % %     [driversBC , idxBestDrivers] = sort(BC , 'descend'  );
                %
                % bestDrivers{1,kScale} = idxBestDrivers;

                % idxBestDrivers = (1:n);


                %% build matrix B

                nDrivers = 1;
                '''% %%%%%% Drivers %%%%%%%%%''';

                %         driversInd = 1: nDrivers;
                %         driversInd = floor(linspace(1 , nExtra , nDrivers));

                % driversInd = idxBestDrivers(1: nDrivers);


                driversInd = drivers2try(kDriver);
                %%%%%%%%%%%%%%%%%
                Drivers = zeros(1,n);
                Drivers(driversInd) = 1;
                B = zeros(n);

                B = B - diag(diag(B)) + diag(Drivers);

                %% Control eigen-modes

                    [yBarSimu, xSimuext, Ucan , P_adj , n_err , condPinvSpec , condAtilde] = ...
                        myControlFun_betterThanStiso_output(A, B , T,STEP, t, x0, xf, rhoOpt,  zeros(n,n) , R , zeros(n,n) , Cbar , zf  );
                
                % Trajectory
                enNumEigen = controlEnergy(Ucan , STEP);

                relativeErrorArraySpectralAllDrivers(1,kDriver) = norm(yBarSimu(end , :)' - zf , dimNormLow)/norm(zf , dimNormLow);
                relativeErrorArraySpectralNotNormAllDrivers(1,kDriver) = norm(yBarSimu(end , :)' - zf , dimNormLow);
                relativeErrorArraySpectralCosineAllDrivers(1,kDriver)  = cosineDist(yBarSimu(end,:)' , zf);
                energyArraySpectralAllDrivers(1,kDriver)  = enNumEigen;
                %%%%%%%%%%%%%
                relativeErrorArraySpectral2OriginalAllDrivers(1,kDriver)  = norm(xSimuext(end,target)' - yf , dimNorm)/norm(yf , dimNorm);
                relativeErrorArraySpectral2OriginalNotNormAllDrivers(1,kDriver)  = norm(xSimuext(end,target)' - yf , dimNorm);
                relativeErrorArraySpectral2OriginalCosineAllDrivers(1,kDriver)  = cosineDist(xSimuext(end,target)' , yf);



            end
            %% store average of drivers
            % notTarget = setdiff(1:n , target);

            driversInIdx = find(drivers2try<=targetSize);
            driversOutIdx = setdiff(1:nDrivers2try , driversInIdx);

            relativeErrorArraySpectralIn(kTask , kScale) = mean(relativeErrorArraySpectralAllDrivers( driversInIdx));
            relativeErrorArraySpectralNotNormIn(kTask , kScale) = mean(relativeErrorArraySpectralNotNormAllDrivers( driversInIdx));
            relativeErrorArraySpectralCosineIn(kTask , kScale) = mean(relativeErrorArraySpectralCosineAllDrivers( driversInIdx ));
            energyArraySpectralIn(kTask , kScale) = mean(energyArraySpectralAllDrivers( driversInIdx));
            relativeErrorArraySpectral2OriginalIn(kTask , kScale) = mean(relativeErrorArraySpectral2OriginalAllDrivers( driversInIdx ));
            relativeErrorArraySpectral2OriginalNotNormIn(kTask , kScale) = mean(relativeErrorArraySpectral2OriginalNotNormAllDrivers( driversInIdx ));
            relativeErrorArraySpectral2OriginalCosineIn(kTask , kScale) = mean(relativeErrorArraySpectral2OriginalCosineAllDrivers( driversInIdx ));


            relativeErrorArraySpectralOut(kTask , kScale) = mean(relativeErrorArraySpectralAllDrivers( driversOutIdx ));
            relativeErrorArraySpectralNotNormOut(kTask , kScale) = mean(relativeErrorArraySpectralNotNormAllDrivers( driversOutIdx ));
            relativeErrorArraySpectralCosineOut(kTask , kScale) = mean(relativeErrorArraySpectralCosineAllDrivers( driversOutIdx ));
            energyArraySpectralOut(kTask , kScale) = mean(energyArraySpectralAllDrivers( driversOutIdx ));
            relativeErrorArraySpectral2OriginalOut(kTask , kScale) = mean(relativeErrorArraySpectral2OriginalAllDrivers( driversOutIdx ));
            relativeErrorArraySpectral2OriginalNotNormOut(kTask , kScale) = mean(relativeErrorArraySpectral2OriginalNotNormAllDrivers( driversOutIdx ));
            relativeErrorArraySpectral2OriginalCosineOut(kTask , kScale) = mean(relativeErrorArraySpectral2OriginalCosineAllDrivers( driversOutIdx ));



        end
    end
    warning on
end
warning on
delete(gcp('nocreate'));

fprintf('computation done ...\n')

%% save workSpace


save(sprintf(...
    "ClusterWorkSpaceSave/single_driver_vary_target_allDriversInOut_rho4_xf10/varyTarget_driver_central_changing_ordLapXf_rhoOpt_%.2e_Control_centralDrivers_out_r4_outputs_%i_n_%i_Xf_10_runID_%s.mat" , ...
    rhoOpt  , nOutputs , n  , runId));

fprintf('WorkSpace saved ...\n')
