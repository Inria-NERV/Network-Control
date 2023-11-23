clear all

%% add path 
subdir = genpath("fig_utils");
addpath(subdir);
subdir = genpath("results_data");
addpath(subdir);

%% 
%%%% rho 4 sgma 10
load('results_data/singleDrivers_new_onlyEigen_hier_rho-04_Control_ordLapXf_1_averageOver_400_Rep_n_256_Xf_mu1_sigma10_Tf_1.mat')
idx2keep  =[1 2 3 5 9 12 13 14];
%
% %%%% rho 4 sgma 1
% load('results_data/singleDrivers_new_onlyEigen_hier_rho-04_Control_ordLapXf_1_averageOver_100_Rep_n_256_Xf_mu1_sigma1_Tf_1.mat')
% %%%% rho 4 sgma 1.5
% % load('singleDrivers_new_onlyEigen_hier_rho-04_Control_ordLapXf_1_averageOver_100_Rep_n_256_Xf_mu1_sigma1point5_Tf_1.mat')
% % comment this if sigma10
% idx2keep = 1:length(outputsComponents);


%%%% rho 1 sgma 1
% load('results_data/singleDrivers_new_onlyEigen_hier_rho-00_Control_ordLapXf_1_averageOver_100_Rep_n_256_Xf_mu1_sigma1_Tf_1.mat')
% idx2keep = 1:length(outputsComponents);

%%%% other TOPO
% load('results_data/singleDrivers_ER_new_onlyEigen_hier_rho-04_Control_ordLapXf_1_averageOver_100_Rep_n_256_Xf_mu1_sigma10_Tf_1.mat')
% load('results_data/singleDrivers_BA_new_onlyEigen_hier_rho-04_Control_ordLapXf_1_averageOver_100_Rep_n_256_Xf_mu1_sigma10_Tf_1.mat')

% idx2keep = 1:length(outputsComponents);

%% keep only power of 2 x


xTick = outputsComponents(idx2keep);

lowMet = relativeErrorArraySpectralCosine2Plot(:,idx2keep);
highMet = relativeErrorArraySpectral2OriginalCosine2Plot(:,idx2keep);
mixMet = lowMet + highMet;

%% main Pic  
mkSize = 20;
lWd = 4;

fig = figure;

fig.Position(3:4) = 250*[1 1];
fig.Color = [1 1 1];


nTicks = length(xTick);

ppLapHigh = plot((xTick) , mean(highMet));
hold on;

ppLapHigh.Marker = '.';
ppLapHigh.MarkerSize = mkSize;


ppLapLow= plot((xTick) , mean(lowMet));
hold on;

ppLapLow.Marker = '.';
ppLapLow.MarkerSize = mkSize;

ppLapHigh.LineStyle = '-';



hold on;
ppLapMix = plot((xTick) , mean(mixMet ));

ppLapMix.Marker = 'none';
ppLapMix.MarkerSize = mkSize;



ppLapLow.Color ='k';
ppLapHigh.Color = 0.5*[1 1 1];
ppLapMix.Color ='k';
ppLapMix.LineStyle = '--';

ppLapHigh.LineWidth=  lWd;
ppLapLow.LineWidth=  lWd;
ppLapMix.LineWidth=  3;


ax = gca;
ax.XScale = "log";
% ax.YScale = "log";

xlim([0 258])

% ylim([0 1.1])


myLabels = {'\delta' , '\eta' , '\delta+\eta'};
leg = legend(ax , [ppLapLow , ppLapHigh ,ppLapMix], myLabels , 'Box','off'...
    ,'Location' , 'Northwest');

leg.FontSize = 15;

ax.Box = 'off';
ax.FontSize = 18;

ax.XTick = [1 4 16 64 256];

ax.XDir = "reverse";

ax.FontSize = 18;

ax.XTickLabelRotation = 0;

