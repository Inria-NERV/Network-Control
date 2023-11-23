clear all;

%% add path 
subdir = genpath("fig_utils");
addpath(subdir);
subdir = genpath("results_data");
addpath(subdir);

%% 
E = 2.5;

sz_cl = 2;
mx_lvl = 8;

%%%%%%%
% params = struct;
% params.directed = false;
% params.force_full = true;
% params.z  = myCommunities;
% params.p = pi;
% params.q = pe;
%%%%%%%%%%%%%%%

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

%% plot 
figure; imagesc(1-matrix); (colormap("gray"))


axis square
ax = gca;
ax.XTick  = 2.^[3 5 6 7 8];
ax.YTick  = 2.^[3 5 6 7 8];

% ax.XScale = "log";
% ax.YScale = "log";
ax.FontSize = 14;

%% colorbar

figure;
imagesc(1:7) ;colorbar ; colormap("hot")


%% plot  markers
figure; spy(matrix , 'k'); 


axis square
ax = gca;
ax.XTick  = 2.^[3 5 6 7 8];
ax.YTick  = 2.^[3 5 6 7 8];

% ax.XScale = "log";
% ax.YScale = "log";
ax.FontSize = 16;