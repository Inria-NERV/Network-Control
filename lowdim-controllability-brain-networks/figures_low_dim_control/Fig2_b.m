
clear all
%% add path 
subdir = genpath("fig_utils");
addpath(subdir);
subdir = genpath("results_data");
addpath(subdir);

%% load bigger Network

load('results_data/finalArithmetic_und_hierarchical_normStable_rhoOpt_1.00e-04_Control_ordLapXf1_centralDrivers_Nd_1_8_64_averageOver_200_Rep_n_256_X0_std0_Xf_std10_Tf_1.mat')

% load('results_data/finalArithmetic_und_hierarchical_normStable_rhoOpt_1.00e-00_Control_ordLapXf1_centralDrivers_Nd_1_8_64_averageOver_100_Rep_n_256_X0_std0_Xf_std1_Tf_1.mat')

%% Plotting

selectDriversTask = [1 2 3];

NdArrayOrig = NdArray;

NdArray2plot = NdArrayOrig(selectDriversTask);

nTasks = length(selectDriversTask);

xTicks = 1:nOutputs;
% xTicks = outputsComponents/nOrig*100;

xTicks = outputsComponents; %/nOrig*100;

myCord = colororder;

letterIdx = {'a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)', 'h)'};

nCol = 2 ; lWd = 4; ftSize = 16;  legFtSize = 12;

legend1 = { 'Average' , 'Eigen'  };
for k = 1:length(NdArray2plot)
    legend1{1,k+2} = sprintf('n_d = %i' , NdArray2plot(k)) ;
end


%% Cosine
ldW = 4;
mkSz = 20;

cOrd = colororder;

fig = figure;
% fig.Position = [1981 364 300 330];

fig.Position(3:4) = 250*[1 1];
fig.Color = [1 1 1];

lapCos = (relativeErrorArraySpectralCosine2Plot(selectDriversTask,:)')/2;


statesCos = lapCos;

% axis square
ax = gca;

hold on ;
ppK = plot( xTicks ,lapCos , 'Marker', '.', 'MarkerSize' , mkSz,  'Color', 'k');
hold on;
% ppAvgK = plot(xTicks , avgCos , 'Marker', '.', 'MarkerSize' , mkSz,  'Color', 'k');
% hold on ;
pp = plot( xTicks , lapCos , 'Marker', '.', 'MarkerSize' , mkSz);
% hold on;
% ppAvg = plot(xTicks , avgCos , 'Marker', '.', 'MarkerSize' , mkSz);


hold on ;
ppAvgK = plot( xTicks ,statesCos , 'Marker', 'none', 'MarkerSize' , mkSz,  'Color', 'k' , 'LineStyle' , '--');
hold on;
% ppAvgK = plot(xTicks , avgCos , 'Marker', '.', 'MarkerSize' , mkSz,  'Color', 'k');
% hold on ;
ppAvg = plot( xTicks , statesCos , 'Marker', '.', 'MarkerSize' , mkSz);


% bb = 5
nTasks = length(selectDriversTask);
for k = 1 : nTasks
    pp(k).Color = myCord(k,:);
    ppK(k).LineWidth = lWd  ;
    ppK(k).Marker = 'none'  ;
    pp(k).LineWidth = lWd  ;

    ppAvg(k).Color = myCord(k,:);
    ppAvgK(k).LineWidth = lWd  ;
    ppAvgK(k).Marker = 'none'  ;
    ppAvg(k).LineWidth = lWd  ;

end



ax = gca;
% ax.FontWeight = 'bold';

ax.YScale = 'log';
ax.XScale = 'log';


ppArray = [ ];

for kTask = 1:nTasks
    ppArray = [ppArray  ,  pp(kTask)  ];
end

ppArray = [ppArray, ppK(1) , ppAvgK(1)];

legend1 = {    };
for k = 1:(nTasks)
    legend1{1,k} = sprintf('n_d = %i' , NdArray2plot(k)) ;
end

plot([1 256] , [0.5 0.5] , 'LineStyle','--' , 'Color' , 'k')

legend1{1,4} =  'eigenmaps'   ;

legend1{1,5} =  'clusters'  ;

leg = legend(ppArray,...
    legend1,'Location' ,'northeast', ...
    'FontSize' , legFtSize , 'Box' , 'off'  ,'NumColumns' , 1 );

leg.FontWeight = 'normal';


leg.Position(4) = 0.5;
ax.FontSize = 14;



%%%%%%%% Cosine eigen and all

cOrd = colororder;

fig = figure;
% fig.Position = [1981 364 300 330];

fig.Position(3:4) = 260*[1 1];
fig.Color = [1 1 1];

nTicks = length(xTicks);


ppArray = [ ];

for kTask = 1 : nTasks
    lapCosArray = (relativeErrorArraySpectralCosine2PlotArray(kTask,:,:));

    lapCosArray = reshape(lapCosArray , [nOutputs nRepetitions ])';


    pLapDis = plot(xTicks ,mean(lapCosArray(:,(1:nTicks))) ,'Color' , cOrd(kTask , :));

    hold on;

    pLapDis.LineWidth = ldW ;
    pLapDis.Marker = '.';
    pLapDis.MarkerSize = mkSz;
    ppArray = [ppArray    pLapDis ];

end

ax = gca;


ppArray = [ppArray ];


legend1 = {};
for k = 1:length(NdArray2plot)
    legend1{1,k} = sprintf('n_d = %i' , NdArray2plot(k)) ;
end

leg = legend(ax , ppArray,...
    legend1,'Location' ,'southeast', ...
    'FontSize' , legFtSize , 'Box' , 'off'  ,'NumColumns' , 1 );

leg.FontWeight = 'normal';


ax.Box = 'off';
ax.XTick = [1 4 16 64 256];

ax.Box = 'off';


leg.Position(4) = 0.25;
leg.FontSize = 14;
ax.FontSize = 18;

% xlim(ax,[0.9 560])
ylim(ax,[0 1])

% ax.YScale = 'linear';
% ax.YScale = 'log';

ax.XScale = 'log';
ax.XLim = [0.7 257];


ax.XDir = "reverse";



%% Cosine low and high
ldW = 2.3;
mkSz = 20;

cOrd = colororder;

fig = figure;
% fig.Position = [1981 364 300 330];

fig.Position(3:4) = 250*[1 1];
fig.Color = [1 1 1];

lapCos = (relativeErrorArraySpectralCosine2Plot(selectDriversTask,:)')/2;


statesCos = lapCos;

% axis square
ax = gca;

hold on ;
ppK = plot( xTicks ,lapCos , 'Marker', '.', 'MarkerSize' , mkSz,  'Color', 'k');
hold on;
% ppAvgK = plot(xTicks , avgCos , 'Marker', '.', 'MarkerSize' , mkSz,  'Color', 'k');
% hold on ;
pp = plot( xTicks , lapCos , 'Marker', '.', 'MarkerSize' , mkSz);
% hold on;
% ppAvg = plot(xTicks , avgCos , 'Marker', '.', 'MarkerSize' , mkSz);


hold on ;
ppAvgK = plot( xTicks ,statesCos , 'Marker', '.', 'MarkerSize' , mkSz,  'Color', 'k' , 'LineStyle' , '--');
hold on;
% ppAvgK = plot(xTicks , avgCos , 'Marker', '.', 'MarkerSize' , mkSz,  'Color', 'k');
% hold on ;
ppAvg = plot( xTicks , statesCos , 'Marker', '.', 'MarkerSize' , mkSz);


% bb = 5
nTasks = length(selectDriversTask);
for k = 1 : nTasks
    pp(k).Color = myCord(k,:);
    ppK(k).LineWidth = lWd  ;
    ppK(k).Marker = 'none'  ;
    pp(k).LineWidth = lWd  ;

    ppAvg(k).Color = myCord(k,:);
    ppAvgK(k).LineWidth = lWd  ;
    ppAvgK(k).Marker = 'none'  ;
    ppAvg(k).LineWidth = lWd  ;

end



ax = gca;
% ax.FontWeight = 'bold';

ax.YScale = 'log';
ax.XScale = 'log';


ppArray = [ ];

for kTask = 1:nTasks
    ppArray = [ppArray  ,  pp(kTask)  ];
end

ppArray = [ppArray, ppK(1) , ppAvgK(1)];

legend1 = {    };
for k = 1:(nTasks)
    legend1{1,k} = sprintf('n_d = %i' , NdArray2plot(k)) ;
end

plot([1 256] , [0.5 0.5] , 'LineStyle','--' , 'Color' , 'k')

legend1{1,4} =  'eigenmaps'   ;

legend1{1,5} =  'clusters'  ;

leg = legend(ppArray,...
    legend1,'Location' ,'northeast', ...
    'FontSize' , legFtSize , 'Box' , 'off'  ,'NumColumns' , 1 );

leg.FontWeight = 'normal';


leg.Position(4) = 0.5;
ax.FontSize = 14;



%%%%%%%% Cosine eigen and all

cOrd = colororder;

fig = figure;
% fig.Position = [1981 364 300 330];

fig.Position(3:4) = 260*[1 1];
fig.Color = [1 1 1];

nTicks = length(xTicks);


ppArray = [ ];

for kTask = 1 : nTasks
    lapCosArray = (relativeErrorArraySpectralCosine2PlotArray(kTask,:,:));

    lapCosArray = reshape(lapCosArray , [nOutputs nRepetitions ])';

    lapCosArrayHigh = (relativeErrorArraySpectral2OriginalCosine2PlotArray(kTask,:,:));

    lapCosArrayHigh = reshape(lapCosArrayHigh , [nOutputs nRepetitions ])';

    pLapDis = plot(xTicks ,mean(lapCosArray(:,(1:nTicks))) ,'Color' , cOrd(kTask , :));

    pLapDis.LineWidth = ldW ;
    pLapDis.Marker = '.';
    pLapDis.MarkerSize = mkSz;
    ppArray = [ppArray    pLapDis ];

    hold on;

    pLapDisHigh = plot(xTicks ,mean(lapCosArrayHigh(:,(1:nTicks))) ,'Color' , cOrd(kTask , :),'LineStyle','--');

    pLapDisHigh.LineWidth = ldW ;
    pLapDisHigh.Marker = 'none';
    pLapDisHigh.MarkerSize = mkSz;


    hold on;

    pLapDisMix = plot(xTicks ,mean(lapCosArrayHigh+lapCosArray) ,'Color' , cOrd(kTask , :),'LineStyle','-.');
    
    pLapDisHigh.LineWidth = ldW ;
    pLapDisHigh.Marker = '.';
    pLapDisHigh.MarkerSize = mkSz;


end

ax = gca;


ppArray = [ppArray ];


legend1 = {};
for k = 1:length(NdArray2plot)
    legend1{1,k} = sprintf('n_d = %i' , NdArray2plot(k)) ;
end

leg = legend(ax , ppArray,...
    legend1,'Location' ,'southeast', ...
    'FontSize' , legFtSize , 'Box' , 'off'  ,'NumColumns' , 1 );

leg.FontWeight = 'normal';


ax.Box = 'off';
ax.XTick = [1 4 16 64 256];

ax.Box = 'off';


leg.Position(4) = 0.25;
leg.FontSize = 14;
ax.FontSize = 14;

% xlim(ax,[0.9 560])
ylim(ax,[0 1])

% ax.YScale = 'linear';
% ax.YScale = 'log';

ax.XScale = 'log';
ax.XLim = [0.7 257];


ax.XDir = "reverse";

