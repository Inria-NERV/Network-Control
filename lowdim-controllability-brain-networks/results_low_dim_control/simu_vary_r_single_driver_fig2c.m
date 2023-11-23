% script to run the second simulation of the paper : varying output
% dimension with single drivers
%
clear all;
delete(gcp('nocreate'));

fprintf('Launch Script for Simulation 1\n');
allSubDirectories = genpath('utils');
addpath(allSubDirectories);

path = pwd


%% get environment variable
% you can call this script from a .sh with the following variables


jobId = getenv('SLURM_ARRAY_JOB_ID');
runId = getenv('SLURM_ARRAY_TASK_ID');
fprintf('arrayID = %s \n' , runId);
nCpuOnNode= str2double(getenv('SLURM_CPUS_ON_NODE'));


rhoExp = str2double( getenv('RHO_EXP')); % regularization parameter

% uncomment here and choose Id run = random seed ; default rho = 4 and gramControlVar = 0;
runId  = '1';
% nCpuOnNode = 4;
rhoExp = 4;


%%%% set the Random Numbers seed %%%%%%%%
rng(str2double(runId))
nCpus = nCpuOnNode;


%% Network Options


ordLapXf= 1;


directedNet = 0;
rhoOpt = 10^(-rhoExp);


n = 256;
outputsComponents =  [1 2 4 8 16 32 64 128 256];



%% %------------------ Time constants ------------------------------
T = 1; %in s
STEP = 0.001 ; %in s
t = 0:STEP:T;
nSamps = length(t);

%% fix initial and final states


sigSpred = 10; % sigmaf
x0 = zeros(n , 1) + sigSpred *(randn(n , 1));
xf = ones(n,1) + sigSpred *(randn(n , 1));


%%  Target Stuff
'''% %%%%%% myTarget %%%%%%%%%''';

target = 1 :n;
% target = 2:20;
targetSize = length(target);
notTarget = setdiff(1:n , target);

%%%%%%%%%%%%%%%%%%%%%%
%%%%% C matrix for target Control %%%%%%%
C = zeros(targetSize , n);
for ktarget = 1 : targetSize
    C(ktarget , target(ktarget)) = 1;
end


%% define arrays


% nTasks =length(NdArray);
nTasks =n;

nOutputs = length(outputsComponents);

relativeErrorArraySpectral = zeros(nTasks,nOutputs);
energyArraySpectral = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%
relativeErrorArraySpectral2Original = zeros(nTasks,nOutputs);

relativeErrorArrayAvg2Original = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%
errorAllTargetsOrigin = zeros(nTasks,nOutputs);
energyAllTargetsOrigin = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%

%%%%%  not normalised errors.%%%%%%
relativeErrorArraySpectralNotNorm = zeros(nTasks,nOutputs);
energyArraySpectralNotNorm = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%
relativeErrorArraySpectral2OriginalNotNorm = zeros(nTasks,nOutputs);

%%%%%%%%%%%%%%%%%%%%%%
errorAllTargetsOriginNotNorm = zeros(nTasks,nOutputs);
energyAllTargetsOriginNotNorm = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
%%%%% Cosine distances %%%%%%
% relativeErrorArraySpectral = zeros(nTasks,nOutputs);
relativeErrorArraySpectralCosine = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%
relativeErrorArraySpectral2OriginalCosine = zeros(nTasks,nOutputs);


errorAllTargetsOriginCosine = zeros(nTasks,nOutputs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


graphs = cell(nTasks,nOutputs);
communities = cell(1,nOutputs);
eigenmaps = cell(1,nOutputs);
bestDrivers = cell(1,nOutputs);


modularityArray = zeros(1,nOutputs);
densityArray = zeros(1,nOutputs);


%% define modular graphs parameters

Diag = 0;
%%%%
rho = 0.05; % rho = 0.2;

L = floor(rho*(n-1)*n/2); % for undirected
% L = floor(rho*(nOrig-1)*nOrig); % for directed

Karray = floor(outputsComponents);
N1array = floor(1./Karray * n);


R = eye(n,n);
Q = zeros(n,n);
yf = (xf(target));
% dimNorm = max([2 length(yf)]);


kScale = 1;

myCommunities = ones(1,n);

%% generate modular graph

mx_lvl = 8;

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

    G = graph(matrix) ;
    connBins = conncomp(G);
    differentComp = find(connBins ~= 1 , 1);
    isConnected = isempty(differentComp);
end

Ggsp = gsp_graph(double(matrix));

%%%%%%%%%%%%%

densityArray(kScale) = sum(sum(matrix))/(n^2);


G = nan;

graphs{ 1, kScale} = sparse(matrix);


%% normalize matrix A and relax

lambdaMax = eigs(matrix , 1);
A = matrix - 1.001 * lambdaMax* eye(n);

% the coef of 1.001 is to ensure Re(Lambda(A))<0 to be stable


%% rank drivers acording to betweenness centrality
%     BC =betweenness_wei(1./matrix);

BC =betweenness_bin(matrix) /(n-1)/(n-2);
[driversBC , idxBestDrivers] = sort(BC , 'ascend');
%     [driversBC , idxBestDrivers] = sort(BC , 'descend'  );


idxBestDrivers = 1:n;


bestDrivers{1,kScale} = idxBestDrivers(1:n);

%%  define eigenModes of the target
%     V = gsp_laplacian_eigenmaps(Ggsp, nOrig); %% not normalised laplacian

Ggsp = gsp_create_laplacian(Ggsp, 'combinatorial');
% Ggsp = gsp_create_laplacian(Ggsp, 'normalized');
%     Ggsp = gsp_create_laplacian(Ggsp, 'chung');

undLap = full(Ggsp.L);
targetlaplac = undLap;


[VnotSorted,lambdaMatNotSorted] = eig( targetlaplac );
lambdaNotSorted = diag(lambdaMatNotSorted);

[lambda , idxAsc] = sort(lambdaNotSorted , 'ascend'); % for Laplacian
V = VnotSorted(: , idxAsc);


%% order eigenvectors

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
matrix = nan; G = nan;


eigenmaps{1,kScale} = V;


%% loop over Community outputs
% orderedModesIdx = 1 : n; %  low to high
dimNormLow = 2;
dimNormHigh = 2;
dimNorm = 2;

nDriversMax = n;



warning off

% you can choose to run it with parallel processing

%%%% deal with Parrallel Obj
% delete(gcp('nocreate'));
% poolObj = parpool('local' ,nCpus);
% % poolObj = parpool('SlurmProfile1' ,nCpus);
%
% parfor kScale = 1 : nOutputs
for kScale = 1 : nOutputs

    %% adapt scale
    rng(100 * str2double(runId) + kScale )

    warning off

    N1 = N1array(kScale);
    K = Karray(kScale);

    modes2take = orderedModesIdx(1 : outputsComponents(kScale));


    Cbar = (V( : ,modes2take)')* C;
    zf = Cbar *yf; % zf = GFT of yf




    %% loop over Tasks
    for kDriver = 1 : nTasks
        fprintf('\n scale: %i : task %i/%i ;\n' , kScale , kDriver , nTasks);
        %% build matrix B
        '''% %%%%%% Drivers %%%%%%%%%''';

        driversInd = idxBestDrivers(kDriver);

        %%%%%%%%%%%%%%%%%
        Drivers = zeros(1,n);
        Drivers(driversInd) = 1;
        B = zeros(n);

        B = B - diag(diag(B)) + diag(Drivers);



        %% control all targets
        if kScale == 1

            [ySimu, xSimu, U , P_adj ,n_err , condPinv , condAtilde] = ...
                myOutputControlFunction(A, B, T, t, x0, xf, rhoOpt , R , Q , C , yf  );

            % Trajectory
            %             enNumSpan = controlEnergy(U , STEP);
            errorAllTargetsOrigin(kDriver , kScale) = norm(ySimu(end,:)' - yf , dimNorm)/norm(yf, dimNorm);
            errorAllTargetsOriginNotNorm(kDriver , kScale) = norm(ySimu(end,:)' - yf , dimNorm);
            errorAllTargetsOriginNotCosine(kDriver , kScale) = cosineDist(ySimu(end,:)', yf);
            energyAllTargetsOrigin(kDriver , kScale) = controlEnergy(U , STEP);
        end


        %% Control eigen-modes

        [yBarSimu, xSimuext, Ucan , P_adj , n_err , condPinvSpec , condAtilde] = ...
            myOutputControlFunction(A, B , T, t, x0, xf, rhoOpt , R , Q , Cbar , zf  );
        % Trajectory
        enNumEigen = controlEnergy(Ucan , STEP);

        % %%%%%  Store Results  Eigen  %%%%%%

        relativeErrorArraySpectral(kDriver , kScale) = norm(yBarSimu(end , :)' - zf , dimNormLow)/norm(zf , dimNormLow);
        relativeErrorArraySpectralNotNorm(kDriver , kScale) = norm(yBarSimu(end , :)' - zf , dimNormLow);
        relativeErrorArraySpectralCosine(kDriver , kScale) = cosineDist(yBarSimu(end,:)' , zf);

        %         localityArraySpectral(kTask , kScale) = localRateEigen;
        energyArraySpectral(kDriver , kScale) = enNumEigen;
        %%%%%%%%%%%%%
        relativeErrorArraySpectral2Original(kDriver , kScale) = norm(xSimuext(end,target)' - yf , dimNorm)/norm(yf , dimNorm);
        relativeErrorArraySpectral2OriginalNotNorm(kDriver , kScale) = norm(xSimuext(end,target)' - yf , dimNorm);
        relativeErrorArraySpectral2OriginalCosine(kDriver , kScale) = cosineDist(xSimuext(end,target)' , yf);



    end
    warning on
end
warning on
delete(gcp('nocreate'));

fprintf('computation done ...\n')

%% save workSpace


save(sprintf(...
    "ClusterWorkSpaceSave/simu2/und_Hier_ordLapXf_%i_rhoOpt_%.2e_Control_outputs_%i_n_%i_percent_Xf_sig_%.2f__runID_%s.mat" , ...
    ordLapXf , rhoOpt  , nOutputs , n  ,sigSpred,  runId));

fprintf('WorkSpace saved ...\n')
