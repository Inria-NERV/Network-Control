

%% load bigger Network
clear all

%% add path 
subdir = genpath("fig_utils");
addpath(subdir);
subdir = genpath("results_data");
addpath(subdir);
%%
% load('results_data/allDriversInOut_r_2_8_32_und_hierarchical_normStable_rhoOpt_1e-04_Control_ordLapXf1_averageOver_400_Rep_n_256_Xf10_Tf_1.mat')

load('results_data/allDriversInOut_r_4_und_hierarchical_normStable_rhoOpt_1e-04_Control_ordLapXf1_averageOver_400_Rep_n_256_Xf10_Tf_1.mat')



%% Plotting

selectTask = [1 2 3];

selectTask = 1;


nTasks = length(selectTask);


xTicks = targetSizeArray;

%xTicks = targetSizeArray;

nTicks = length(xTicks);



myCord = colororder;


nCol = 2 ; lWd = 2; ftSize = 13;  legFtSize = 14;

legend1 = { 'target' , '\delta' , '\eta'  };
for k = 1:nTasks
    legend1{1,k+3} = sprintf('r = %i' , outputDimArray(k)) ;
end


% 

%% plot the sum

mkSz = 20;
lWd = 4;



fig = figure;

% fig.Position(3:4) = [814 687];
fig.Position(3:4) = 250*[1 1];
fig.Color = [1 1 1];

relativeErrorArraySpectralCosineIn2plot(relativeErrorArraySpectralCosineIn2plot==0) = nan;
relativeErrorArraySpectral2OriginalCosineIn2plot(relativeErrorArraySpectral2OriginalCosineIn2plot==0) = nan;

relativeErrorArraySpectralCosineOut2plot(relativeErrorArraySpectralCosineOut2plot==0) = nan;
relativeErrorArraySpectral2OriginalCosineOut2plot(relativeErrorArraySpectral2OriginalCosineOut2plot==0) = nan;


mixMetricIn = relativeErrorArraySpectralCosineIn2plot + relativeErrorArraySpectral2OriginalCosineIn2plot;
mixMetricOut = relativeErrorArraySpectralCosineOut2plot + relativeErrorArraySpectral2OriginalCosineOut2plot;

% subplot(1,2,1)

axis square
ax = gca;


mixMetric(2,1) = nan;


ppIn = plot( xTicks ,mixMetricIn(selectTask,(1:nTicks))' , 'LineStyle', '--', 'MarkerSize' , mkSz);


hold on;
ppOut = plot( xTicks ,mixMetricOut(selectTask,(1:nTicks))' , 'LineStyle', '--', 'MarkerSize' , mkSz);


% hold on;
% 
% ppAll = plot(xTicks , 2*errorAllTargetsOriginCosine2Plot(1,(1:nTicks))' , 'Marker', 'diamond', 'MarkerSize' , mkSz,  'Color', 'k', 'LineWidth',lWd,'LineStyle',':');


% bb = 5
markers = {'^','s' , 'd'};

nTasks = length(selectTask);
for k = 1 : nTasks
    % ppIn(k).Color = 0.2*[1 1 1];
    % ppIn(k).MarkerFaceColor = 0.2*[1 1 1];

    ppIn(k).LineWidth = lWd  ;
    ppIn(k).Marker = '.' ;

    ppOut(k).LineWidth = lWd  ;
    ppOut(k).Marker = '.' ;
    
end

ax = gca;


ppArray = [];

for kTask = 1:nTasks
    ppArray = [ppArray  ,  ppIn(nTasks+1 -kTask)  ];
end

ppArray = [ppIn , ppOut];

% legend1 = { };
% for k = 1:nTasks
%     legend1{1,k} = sprintf('r = %i' , outputDimArray(nTasks+1 -k)) ;
% end

leg = legend(ppArray,...
    {'drivers in' ,'drivers out'} ,'Location' ,'northwest', ...
    'FontSize' , legFtSize , 'Box' , 'off'  ,'NumColumns' , 1 );


ax.Box = 'off';


axis square

ax.XTick = xTicks;
ax.FontSize = 20;
ax.XScale = 'log';

xlim([4 265])
 ylim([0.58 1.4])

ax.XDir = "reverse";





%% subplot 1  IN
% mkSz = 7;
% 
% cOrd = colororder;
% 
% fig = figure;
% % fig.Position = [1981 364 300 330];
% 
% % fig.Position(3:4) = [814 687];
% fig.Position(3:4) = [300 300];
% fig.Color = [1 1 1];
% 
% relativeErrorArraySpectralCosine2Plot(relativeErrorArraySpectralCosine2Plot==0) = nan;
% relativeErrorArraySpectral2OriginalCosine2Plot(relativeErrorArraySpectral2OriginalCosine2Plot==0) = nan;
% 
% % subplot(1,2,1)
% 
% axis square
% ax = gca;
% 
% hold on ;
% pp1K = plot( xTicks ,relativeErrorArraySpectralCosineIn2plot(selectTask,flip(1:nTicks))' , 'Marker', '^', 'MarkerSize' , mkSz,  'Color', 'k');
% hold on;
% ppAvgK = plot(xTicks , relativeErrorArraySpectral2OriginalCosineIn2plot(selectTask,flip(1:nTicks))' , 'Marker', '^', 'MarkerSize' , mkSz,  'Color', 'k');
% hold on ;
% pp = plot( xTicks ,relativeErrorArraySpectralCosineIn2plot(selectTask,flip(1:nTicks))' , 'Marker', '^', 'MarkerSize' , mkSz);
% hold on;
% ppAvg = plot(xTicks , relativeErrorArraySpectral2OriginalCosineIn2plot(selectTask,flip(1:nTicks))' , 'Marker', '^', 'MarkerSize' , mkSz);
% hold on;
% % ppAll = plot(xTicks , errorAllTargetsCosineOrigin2Plot(1,:)' , 'Marker', 'diamond', 'MarkerSize' , mkSz,  'Color', 'k', 'LineWidth',lWd);
% % ppAll = plot(xTicks , errorAllTargetsOriginCosineIn2plot(1,flip(1:nTicks))' , 'Marker', 'diamond', 'MarkerSize' , mkSz,  'Color', 'k', 'LineWidth',lWd,'LineStyle',':');
% 
% 
% % bb = 5
% nTasks = length(selectTask);
% for k = 1 : nTasks
%     pp(k).Color = myCord(k,:);
%     ppAvg(k).Color = myCord(k,:);
%     pp1K(k).LineWidth = lWd  ;
%     ppAvgK(k).LineWidth = lWd  ;
%     ppAvgK(k).LineStyle = '--'  ;
%     ppAvg(k).LineStyle = '--'  ;
% 
%     pp1K(k).Marker = 'none'  ;
%     ppAvgK(k).Marker = 'none'  ;
% 
%     pp(k).LineWidth = lWd  ;
%     ppAvg(k).LineWidth = lWd  ;
% end
% 
% 
% ax = gca;
% 
% 
% ppArray = [ppAll  ,pp1K(1), ppAvgK(1)];
% 
% for kTask = 1:nTasks
%     ppArray = [ppArray  ,  pp(kTask)  ];
% end
% 
% 
% leg = legend(ppArray,...
%     legend1,'Location' ,'north', ...
%     'FontSize' , legFtSize , 'Box' , 'off'  ,'NumColumns' , 2 );
% 
% 
% 
% xlim([3.8 257])
% 
% 
% ax.Box = 'off';
% 
% 
% axis square
% 
% ax.XTick = xTicks;
% ax.XTickLabel = string((xTicks(flip(1:nTicks))));
% 
% ax.FontSize = 14;
% ax.XScale = 'log';
% 
% 
% %% subplot 1  OUT
% mkSz = 7;
% 
% cOrd = colororder;
% 
% fig = figure;
% % fig.Position = [1981 364 300 330];
% 
% % fig.Position(3:4) = [814 687];
% fig.Position(3:4) = [300 300];
% fig.Color = [1 1 1];
% 
% relativeErrorArraySpectralCosine2Plot(relativeErrorArraySpectralCosine2Plot==0) = nan;
% relativeErrorArraySpectral2OriginalCosine2Plot(relativeErrorArraySpectral2OriginalCosine2Plot==0) = nan;
% 
% % subplot(1,2,1)
% 
% axis square
% ax = gca;
% 
% hold on ;
% pp1K = plot( xTicks ,relativeErrorArraySpectralCosineOut2plot(selectTask,flip(1:nTicks))' , 'Marker', '^', 'MarkerSize' , mkSz,  'Color', 'k');
% hold on;
% ppAvgK = plot(xTicks , relativeErrorArraySpectral2OriginalCosineOut2plot(selectTask,flip(1:nTicks))' , 'Marker', '^', 'MarkerSize' , mkSz,  'Color', 'k');
% hold on ;
% pp = plot( xTicks ,relativeErrorArraySpectralCosineOut2plot(selectTask,flip(1:nTicks))' , 'Marker', '^', 'MarkerSize' , mkSz);
% hold on;
% ppAvg = plot(xTicks , relativeErrorArraySpectral2OriginalCosineOut2plot(selectTask,flip(1:nTicks))' , 'Marker', '^', 'MarkerSize' , mkSz);
% hold on;
% % ppAll = plot(xTicks , errorAllTargetsCosineOrigin2Plot(1,:)' , 'Marker', 'diamond', 'MarkerSize' , mkSz,  'Color', 'k', 'LineWidth',lWd);
% % ppAll = plot(xTicks , errorAllTargetsOriginCosineIn2plot(1,flip(1:nTicks))' , 'Marker', 'diamond', 'MarkerSize' , mkSz,  'Color', 'k', 'LineWidth',lWd,'LineStyle',':');
% 
% 
% % bb = 5
% nTasks = length(selectTask);
% for k = 1 : nTasks
%     pp(k).Color = myCord(k,:);
%     ppAvg(k).Color = myCord(k,:);
%     pp1K(k).LineWidth = lWd  ;
%     ppAvgK(k).LineWidth = lWd  ;
%     ppAvgK(k).LineStyle = '--'  ;
%     ppAvg(k).LineStyle = '--'  ;
% 
%     pp1K(k).Marker = 'none'  ;
%     ppAvgK(k).Marker = 'none'  ;
% 
%     pp(k).LineWidth = lWd  ;
%     ppAvg(k).LineWidth = lWd  ;
% end
% 
% 
% ax = gca;
% 
% 
% ppArray = [ppAll  ,pp1K(1), ppAvgK(1)];
% 
% for kTask = 1:nTasks
%     ppArray = [ppArray  ,  pp(kTask)  ];
% end
% 
% 
% leg = legend(ppArray,...
%     legend1,'Location' ,'north', ...
%     'FontSize' , legFtSize , 'Box' , 'off'  ,'NumColumns' , 2 );
% 
% 
% 
% xlim([3.8 257])
% 
% 
% ax.Box = 'off';
% 
% 
% axis square
% 
% ax.XTick = xTicks;
% ax.XTickLabel = string((xTicks(flip(1:nTicks))));
% 
% ax.FontSize = 14;
% ax.XScale = 'log';

