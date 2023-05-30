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
jobId = getenv('SLURM_ARRAY_JOB_ID');
runId = getenv('SLURM_ARRAY_TASK_ID');
fprintf('arrayID = %s \n' , runId);
nCpuOnNode= str2double(getenv('SLURM_CPUS_ON_NODE'));


rhoExp = str2double( getenv('RHO_EXP')); % regularization parameter
gramControlVar = str2double( getenv('GRAM_CONTROL')); % control with Gramian inversion set to 0


% runId  = '1';
% nCpuOnNode = 4;
% gramControlVar = 1;  connectFactor = 4;
% rhoExp = 6;


%%%% set the Random Numbers seed %%%%%%%%
rng(str2double(runId))
nCpus = nCpuOnNode;


%% Options Spec

normOpt = 1; dampFact = 0 ; c = 0 ; %10^(-3);


%% Network Size %%%%%%%%%%%%%%

ordLapXf= 1;
newNormalization = 0;
multiForPropag = connectFactor;

gramControl  = gramControlVar;

directedNet = 0;
rhoOpt = 10^(-rhoExp);



n = 128;
outputsComponents =  [1 2 4 8 16 32 64 128];



%% %------------------ Time constants ------------------------------
T = 1; %in s
STEP = 0.001 ; %in s
t = 0:STEP:T;
nSamps = length(t);

%% fix initial and final states


x0 = 0.5*(zeros(n , 1) +   (rand(n , 1) ));
xf =  0.5*(ones(n , 1) +   (rand(n , 1) ));

% xf(1:32) =  xf(1:32) +  ones(32,1);
xf(1:48) =  xf(1:48) +  ones(48,1);

% xf = zeros(n,1);
% xf(floor(n/2) : n) = ones(length(floor(n/2) : n) , 1);


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

%%%%%%%%%%%%%%%%%%%%%%

%%%%% Cosine distances %%%%%%
% relativeErrorArraySpectral = zeros(nTasks,nOutputs);
relativeErrorArraySpectralCosine = zeros(nTasks,nOutputs);
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

Karray = floor(outputsComponents);
N1array = floor(1./Karray * n);


piArray = 1.5 * log(N1array)./N1array ;
piMinArray = ones(1,nOutputs)./N1array *L/n;
piArray = max([piArray ; piMinArray]);

peArray = (L *2/n - piArray.*(N1array - 1))./(n - N1array);
% peArray = max([peArray ; zeros(1,nOutputs) ]);

if peArray(2) < 0
    peArray(2) = abs(peArray(2));
end

%%%%%%%%%%%%%%%%%
% if first scale is one module
if outputsComponents(1) == 1
    peArray(1) = 0;
    piArray(1) = rho;
end

%%%%%%%%%%%%%%%%%
ziArray = piArray .*(N1array - 1);
zeArray = peArray .*N1array.*(Karray - 1);

Ltheory = n/2*(piArray .*(N1array - 1)  + peArray .*((n-1) - (N1array - 1) ));
%%%%%%%%%%%%%%%%%%%%%%%%%%


R = eye(n,n);
Q = zeros(n,n);
yf = (xf(target));
% dimNorm = max([2 length(yf)]);


kScale = 1;

myCommunities = ones(1,n);

%% generate modular graph


mx_lvl = 7;

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
A = multiForPropag * ((matrix / (1.01 * lambdaMax)) - eye(n));


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



%


%% loop over Community outputs
% orderedModesIdx = 1 : n; %  low to high
dimNormLow = 2;
dimNormHigh = 2;
dimNorm = 2;

nDriversMax = n;



warning off

%%%% deal with Parrallel Obj
delete(gcp('nocreate'));
patchJobStorageLocation
poolObj = parpool('local' ,nCpus);
% poolObj = parpool('SlurmProfile1' ,nCpus);

parfor kScale = 1 : nOutputs
% for kScale = 2:3

    %% adapt scale
    rng(100 * str2double(runId) + kScale )

    warning off

    N1 = N1array(kScale);
    K = Karray(kScale);

    modes2take = orderedModesIdx(1 : outputsComponents(kScale));

    
    Cbar = (V( : ,modes2take)')* C;
    zf = Cbar *yf; % zf = GFT of yf


    %% deal with communities
    fprintf('\n scale %i / %i \n' , kScale , nOutputs);

    scale = outputsComponents(kScale);

    myCommunities = zeros(1 , n);
    for kCom = 1 : K
        if kCom < K
            myCommunities(N1*(kCom - 1)+1 : N1*(kCom - 1)+ N1) = kCom;
        else
            myCommunities(N1*(kCom - 1)+1 : end) = kCom;
        end
    end
    %%%%%%%%%%


    %%%%% C matrix for Communities %%%%%%%
    Cavg = zeros(scale , n);
    for kCom = 1 : scale
        nodesInCom = find(myCommunities == kCom )';
        Cavg(kCom , nodesInCom) = 1/sqrt(length(nodesInCom));

    end
    %     Cavg = Cavg./(sum(Cavg, 2)*ones(1,n));
    %%%%%%%%%%%%%%%%%%%%%
    communities{kScale} = myCommunities;

    yfAvg = Cavg*yf;


    %% retrieve High Dim from low Dim
    reconstructYf = pinv(V( : ,modes2take)') * zf;
    diffReconstructedOriginalPinvSpec( 1,kScale) = norm(yf - reconstructYf , dimNorm)/ norm(yf , dimNorm) ;
    diffReconstructedOriginalPinvSpecCosine( 1,kScale) = (1 - cosineDist(yf , reconstructYf ))/2;


    reconstructYf = pinv(Cavg(:,target)) * yfAvg;
    diffReconstructedOriginalPinvAvg( 1,kScale) = norm(yf - reconstructYf , dimNorm)/ norm(yf , dimNorm) ;
    diffReconstructedOriginalPinvAvgCosine( 1,kScale) = (1 - cosineDist(yf , reconstructYf ))/2;


    %% loop over Tasks
    for kDriver = 1 : nTasks
        fprintf('\n scale: %i : task %i/%i ;\n' , kScale , kDriver , nTasks);
        %% build matrix B
        '''% %%%%%% Drivers %%%%%%%%%''';

        driversInd = idxBestDrivers(kDriver);

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

        %% Get Gramian indices
        
        [W,Tr] = CtrGram(A,B);

        avgGram = (Cavg * W *(Cavg'));
        specGram = (Cbar * W *(Cbar'));

        lambdaAvg = eig(avgGram);
        lambdaSpec = eig(specGram);
        lambdaAll = eig(W);

        lambdaMinGramAvg(kDriver , kScale) = min(lambdaAvg);
        lambdaMinGramSpec(kDriver , kScale) = min(lambdaSpec);

        lambdaMeanGramAvg(kDriver , kScale) = sum(lambdaAvg);
        lambdaMeanGramSpec(kDriver , kScale) = sum(lambdaSpec);

        lambdaMinGramAll(kDriver , kScale) = min(lambdaAll);

        %% control all targets
        if kScale == 1
            if gramControl
                [ySimu, xSimu, U ] = myOutputControlFunction_GramInv(A, B, T, t, x0, yf , C);
            else
                [ySimu, xSimu, U , P_adj ,n_err , condPinv , condAtilde] = ...
                    myOutputControlFunction(A, B, T, t, x0, xf, rhoOpt , R , Q , C , yf  );
            end

            % Trajectory
            %             enNumSpan = controlEnergy(U , STEP);
            errorAllTargetsOrigin(kDriver , kScale) = norm(ySimu(end,:)' - yf , dimNorm)/norm(yf, dimNorm);
            errorAllTargetsOriginNotNorm(kDriver , kScale) = norm(ySimu(end,:)' - yf , dimNorm);
            errorAllTargetsOriginNotCosine(kDriver , kScale) = cosineDist(ySimu(end,:)', yf);
            energyAllTargetsOrigin(kDriver , kScale) = controlEnergy(U , STEP);
        end
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
        relativeErrorArrayAvg(kDriver , kScale) = norm(ySimuAvg(end ,:)' - yfAvg , dimNormLow)/norm(yfAvg , dimNormLow);
        relativeErrorArrayAvgNotNorm(kDriver , kScale) = norm(ySimuAvg(end ,:)' - yfAvg , dimNormLow);
        relativeErrorArrayAvgCosine(kDriver , kScale) = cosineDist(ySimuAvg(end,:)' , yfAvg);

        %         localityArrayAvg = localRateAvg;
        energyArrayAvg(kDriver , kScale) = enNumAvg;
        %%%%%%%%%%%%%
        relativeErrorArrayAvg2Original(kDriver , kScale) = norm(xSimuAverage(end,target)' - yf , dimNorm)/norm(yf , dimNorm);
        relativeErrorArrayAvg2OriginalNotNorm(kDriver , kScale) = norm(xSimuAverage(end,target)' - yf , dimNorm);
        relativeErrorArrayAvg2OriginalCosine(kDriver , kScale) = cosineDist(xSimuAverage(end,target)' , yf);
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
                myOutputControlFunction(A, B , T, t, x0, xf, rhoOpt , R , Q , Cbar , zf  );
        end
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

if gramControl
    save(sprintf(...
    "ClusterWorkSpaceSave/single_driver_newNorm/und_Hier_newNorm_4_ordLapXf_%i_gram_Control_outputs_%i_n_%i_rho_%i_percent_Xf_sig_%.2f_runID_%s.mat" , ...
       ordLapXf , nOutputs , n , floor(100 * rho) , sigSpred, runId));
else
    save(sprintf(...
    "ClusterWorkSpaceSave/single_driver_newNorm/und_Hier_ordLapXf_%i_rhoOpt_%.2e_Control_outputs_%i_n_%i_rho_%i_percent_Xf_sig_%.2f__runID_%s.mat" , ...
    ordLapXf , rhoOpt  , nOutputs , n , floor(100 * rho) ,sigSpred,  runId));
end

fprintf('WorkSpace saved ...\n')
