clear all;

%% add path 
subdir = genpath("fig_utils");
addpath(subdir);
subdir = genpath("results_data");
addpath(subdir);
%%
load('results_data/workSpace_SNR_optimLap_optimAvg_scale64_Nd_256_100Rep.mat')

% xTick = outputsComponents;
xTick = sigmaArray;

alphaVal = 0.2; alphaVal2 = alphaVal;
cOrd = colororder;
matlabBlue = cOrd(1,:);
matlabRed = cOrd(2,:);
matlabYellow = cOrd(3,:);
matlabPurp = cOrd(4,:);

% idxShiftArrayTrigo=partSpinArray;
fig = figure;

fig.Position(3:4) = [260 300];
fig.Color = [1 1 1];

% ppAvg = plot_distribution(xTick , (clusAvgRelativeNotOptimHigh(~isnan(sum(...
% clusAvgCosineNotOptimHigh,2)),:)) , 'Color',matlabBlue,'Alpha', alphaVal);

hold on;
ppLap = plot_distribution(xTick , lapEigsCosineNotOptimHigh , 'Color',matlabRed ,'Alpha', alphaVal2);
hold on;

ppLapOrd = plot_distribution(xTick , lapEigsCosineOptimHigh, 'Color', matlabBlue,'Alpha', alphaVal);
% hold on;
% ppAvgOrd = plot_distribution(xTick , (clusAvgRelativeOptimHigh(~isnan(sum(...
% clusAvgCosineNotOptimHigh,2)),:)) , 'Color',matlabPurp,'Alpha', alphaVal);
% 

xline(1)
myLeg = {'eigen. ord.' , 'eigen. not ord.' };

myLeg = {'y^{EIG} sorted' , '\lambda sorted' };


leg = legend( [ppLapOrd , ppLap    ] , myLeg , 'Box' , 'off'...
    , 'Location','northeast' , 'FontSize',14  , 'NumColumns',1 );

leg.FontWeight = 'normal';

% leg.Position(4) = 0.4;

% xlabel('sigma')
% ylabel('cosine distance')

% xlim([0 257])
ax = gca;

ax.FontSize = 14;

% xlim([0 258])
% ylim([0 0.5])

ax.XScale = 'log';
xlim([0.01 100])
% ylim([0 0.4])

% ax.FontWeight = 'bold';


%% %% inset 
% sizePic = 0.3;
% ax2 = axes('Position',[0.25 0.5 (sizePic)  sizePic]);

% load('/Users/remy.benmessaoud/ownCloud/remy.benmessaoud/MATLAB/outputControllability/Laplacian_rpz/workSpace_optimLap_optimAvg_LaplacianNorm_Nd_256_sigma100_100Rep.mat')
% xTick = outputsComponents;

% ppAvg = plot_distribution(xTick , (clusAvgRelativeNotOptimHigh(~isnan(sum(...
% clusAvgCosineNotOptimHigh,2)),:))/2 , 'Color',matlabBlue, 'Alpha', alphaVal);
% 
% hold on;

fig = figure;
fig.Position(3:4) = [100 100];

ppLap = plot_distribution(xTick , (lapEigsCosineNotOptim) , 'Color',matlabRed,'Alpha', alphaVal2);
hold on;

ppLapOrd = plot_distribution(xTick , (lapEigsCosineOptim), 'Color', matlabBlue,'Alpha', alphaVal);
hold on;
% ppAvgOrd = plot_distribution(xTick , (clusAvgRelativeOptimHigh(~isnan(sum(...
% clusAvgCosineNotOptimHigh,2)),:)) , 'Color',matlabPurp,'Alpha', alphaVal);

% myLeg = {'eigenmaps' , 'structural clusters' , 'eigenmaps-ord.' , 'functional clusters'};
% 
% legend( [ppLap , ppAvg , ppLapOrd , ppAvgOrd ] , myLeg , 'Box' , 'off'...
%     , 'Location','north' , 'FontSize',16  , 'NumColumns',1)

% xlabel('output dimension r')
% ylabel('cosine distance')

% xlim([0 257])
ax = gca;

ax.Box = 'on';
ax.XTick = [0.01 1 100];

ax.YTick = [0.99996 1];

ax.YLim = [0.99995 1];

ax.YTickLabelRotation = 60;

ax.YTickLabel = {'0.99' , '1'};



ax.FontSize = 14;
xlim([0 258])
% ylim([0 0.5])

ax.XScale = 'log';
% ax.XScale = 'log';
% yline(0.5)
xline(1)

axis square