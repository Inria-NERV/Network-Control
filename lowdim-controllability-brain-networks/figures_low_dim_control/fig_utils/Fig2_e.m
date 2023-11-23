

%% load bigger Network
clear all

% load('results_data/ArithmeticMean_driverSameModule_r_4_8_24_und_hierarchical_normStable_rhoOpt_1.00e-04_Control_ordLapXf_averageOver_200_Rep_n_256_Xf10_Tf_1.mat')


% load('/Users/remy.benmessaoud/ownCloud/remy.benmessaoud/figures_low_dim_control/results_data/driverSameModule_und_hierarchical_normStable_rhoOpt_1e-00_Control_ordLapXf1_averageOver_100_Rep_n_256_Xf1_Tf_1.mat')

load('results_data/driverSameModule_und_hierarchical_normStable_rhoOpt_1e-00_Control_ordLapXf1_averageOver_100_Rep_n_256_Xf1_Tf_1.mat')

% same driver 
% load('results_data/driverSameModule_und_hierarchical_normStable_rhoOpt_1e-04_Control_ordLapXf1_averageOver_400_Rep_n_256_Xf10_Tf_1.mat')
% 
% 
% % change driver
% load('/Users/remy.benmessaoud/ownCloud/remy.benmessaoud/MATLAB/outputControllability/vary_target/resultsAverageWrkSpace/driverChangeModule_und_hierarchical_normStable_rhoOpt_1e-04_Control_ordLapXf1_averageOver_111_Rep_n_256_Xf10_Tf_1.mat')
% 
% 

% 
% % all drivers 
% load('/Users/remy.benmessaoud/ownCloud/remy.benmessaoud/MATLAB/outputControllability/vary_target/resultsAverageWrkSpace/allDrivers_und_hierarchical_normStable_rhoOpt_1e-04_Control_ordLapXf1_averageOver_100_Rep_n_256_Xf10_Tf_1.mat')
% 

xTicks = targetSizeArray;

% %network bigger
% load('driverSameModule_und_hierarchical_normStable_rhoOpt_1e-04_Control_ordLapXf1_averageOver_400_Rep_n_512_Xf10_Tf_1.mat')
% xTicks = networkSizeArray;
% 

%% Plotting

selectTask = [1 2 3];

nTasks = length(selectTask);



nTicks = length(xTicks);



myCord = colororder;


nCol = 2 ; lWd = 2; ftSize = 13;  legFtSize = 14;

legend1 = { 'target' , '\delta' , '\eta'  };
for k = 1:nTasks
    legend1{1,k+3} = sprintf('r = %i' , outputDimArray(k)) ;
end


%% subplot 1 
mkSz = 7;

cOrd = colororder;

fig = figure;
% fig.Position = [1981 364 300 330];

% fig.Position(3:4) = [814 687];
fig.Position(3:4) = [300 300];
fig.Color = [1 1 1];

relativeErrorArraySpectralCosine2Plot(relativeErrorArraySpectralCosine2Plot==0) = nan;
relativeErrorArraySpectral2OriginalCosine2Plot(relativeErrorArraySpectral2OriginalCosine2Plot==0) = nan;

% subplot(1,2,1)

axis square
ax = gca;

hold on ;
pp1K = plot( xTicks ,relativeErrorArraySpectralCosine2Plot(selectTask,flip(1:nTicks))' , 'Marker', '^', 'MarkerSize' , mkSz,  'Color', 'k');
hold on;
ppAvgK = plot(xTicks , relativeErrorArraySpectral2OriginalCosine2Plot(selectTask,flip(1:nTicks))' , 'Marker', '^', 'MarkerSize' , mkSz,  'Color', 'k');
hold on ;
pp = plot( xTicks ,relativeErrorArraySpectralCosine2Plot(selectTask,flip(1:nTicks))' , 'Marker', '^', 'MarkerSize' , mkSz);
hold on;
ppAvg = plot(xTicks , relativeErrorArraySpectral2OriginalCosine2Plot(selectTask,flip(1:nTicks))' , 'Marker', '^', 'MarkerSize' , mkSz);
hold on;
% ppAll = plot(xTicks , errorAllTargetsCosineOrigin2Plot(1,:)' , 'Marker', 'diamond', 'MarkerSize' , mkSz,  'Color', 'k', 'LineWidth',lWd);
ppAll = plot(xTicks , errorAllTargetsOriginCosine2Plot(1,flip(1:nTicks))' , 'Marker', 'diamond', 'MarkerSize' , mkSz,  'Color', 'k', 'LineWidth',lWd,'LineStyle',':');



% bb = 5
nTasks = length(selectTask);
for k = 1 : nTasks
    pp(k).Color = myCord(k,:);
    ppAvg(k).Color = myCord(k,:);
    pp1K(k).LineWidth = lWd  ;
    ppAvgK(k).LineWidth = lWd  ;
    ppAvgK(k).LineStyle = '--'  ;
    ppAvg(k).LineStyle = '--'  ;
    
    pp1K(k).Marker = 'none'  ;
    ppAvgK(k).Marker = 'none'  ;
    
    pp(k).LineWidth = lWd  ;
    ppAvg(k).LineWidth = lWd  ;
end


ax = gca;


ppArray = [ppAll  ,pp1K(1), ppAvgK(1)];

for kTask = 1:nTasks
    ppArray = [ppArray  ,  pp(kTask)  ];
end


leg = legend(ppArray,...
    legend1,'Location' ,'north', ...
    'FontSize' , legFtSize , 'Box' , 'off'  ,'NumColumns' , 2 );



xlim([3.8 257])


ax.Box = 'off';


axis square

ax.XTick = xTicks;
ax.XTickLabel = string((xTicks(flip(1:nTicks))));

ax.FontSize = 14;
ax.XScale = 'log';



%% plot the sum

mkSize = 9;
lWd = 2.3


fig = figure;

% fig.Position(3:4) = [814 687];
fig.Position(3:4) = 250*[1 1];
fig.Color = [1 1 1];

relativeErrorArraySpectralCosine2Plot(relativeErrorArraySpectralCosine2Plot==0) = nan;
relativeErrorArraySpectral2OriginalCosine2Plot(relativeErrorArraySpectral2OriginalCosine2Plot==0) = nan;

mixMetric = relativeErrorArraySpectralCosine2Plot + relativeErrorArraySpectral2OriginalCosine2Plot;

% subplot(1,2,1)

axis square
ax = gca;


mixMetric(2,1) = nan;

pp = plot( xTicks ,mixMetric(selectTask,(1:nTicks))' , 'LineStyle', '--', 'MarkerSize' , mkSize);
hold on;

ppAll = plot(xTicks , 2*errorAllTargetsOriginCosine2Plot(1,(1:nTicks))' , 'Marker', 'diamond', 'MarkerSize' , mkSz,  'Color', 'k', 'LineWidth',lWd,'LineStyle',':');


% bb = 5
markers = {'^','s' , 'd'};

nTasks = length(selectTask);
for k = 1 : nTasks
    pp(k).Color = 0.2*[1 1 1];
    pp(k).MarkerFaceColor = 0.2*[1 1 1];

    pp(k).LineWidth = lWd  ;
    pp(k).Marker = markers{1 , k}  ;
    
end

ax = gca;


ppArray = [];

for kTask = 1:nTasks
    ppArray = [ppArray  ,  pp(nTasks+1 -kTask)  ];
end


legend1 = { };
for k = 1:nTasks
    legend1{1,k} = sprintf('r = %i' , outputDimArray(nTasks+1 -k)) ;
end

leg = legend(ppArray,...
    legend1,'Location' ,'northwest', ...
    'FontSize' , legFtSize , 'Box' , 'off'  ,'NumColumns' , 1 );


ax.Box = 'off';


axis square

ax.XTick = xTicks;

ax.FontSize = 14;
ax.XScale = 'log';

xlim([3 280])
% ylim([0.3 1.1])

ax.XDir = "reverse";


