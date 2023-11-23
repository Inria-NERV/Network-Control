
clear all;
load('results_data/vary_rho_wrkSpace_sym.mat')

%% plot 

fig = figure;
% fig.Position = [1981 364 300 330];

% fig.Position(3:4) = [814 687];
fig.Position(3:4) = [260 300];
fig.Color = [1 1 1];

lwdth = 2;
mkSz = 14;

% subplot(1,2,1)

% relativeErrorArraySpectralCosine(relativeErrorArraySpectralCosine==0) = nan;
% relativeErrorArrayAvgCosine2Plot(relativeErrorArrayAvgCosine2Plot==0) = nan;



hold on;
ppLap = plot(1./rhoArray , (relativeErrorArraySpectralCosine) , 'Marker', '.', 'MarkerSize' , mkSz, 'LineWidth',lWd);


% hold on;
% ppAll = plot(rhoArray , (errorAllTargetsOriginNotCosine(:,1)) ,'Color','k', 'Marker', '.', 'MarkerSize' , mkSz, 'LineWidth',lWd);


ax = gca;

ax.XScale = 'log';

ax.FontSize = 14;

% leg = legend(  {'s=8',  's=32', 's=96','s=128','s=160','s=256'},'Location' ,'north', ...
%     'FontSize' , legFtSize , 'Box' , 'off'  ,'NumColumns' ,  1 );



ax.XTick = 10.^[ -15 -10 -5 0 +5 10 15];

% ax.XLim = 10.^[-16 16];

ax.YLim = [-1 1];


xline(10^4 , 'LineStyl'  , '--')
yline(0 , 'LineStyl'  , '-')

leg = legend(  {'r=4',  'r=16', 'r=64','r=256'},'Location' ,'southeast', ...
    'FontSize' , legFtSize , 'Box' , 'off'  ,'NumColumns' ,  2 );

% axis square

%% plot target control


fig = figure;
% fig.Position = [1981 364 300 330];

% fig.Position(3:4) = [814 687];
fig.Position(3:4) = [100 100];
fig.Color = [1 1 1];

lwdth = 2;
mkSz = 8;

% subplot(1,2,1)

% relativeErrorArraySpectralCosine(relativeErrorArraySpectralCosine==0) = nan;
% relativeErrorArrayAvgCosine2Plot(relativeErrorArrayAvgCosine2Plot==0) = nan;



hold on;
ppLap = plot(1./rhoArray , (energyArraySpectral) , 'Marker', '.', 'MarkerSize' , mkSz, 'LineWidth',lWd);

ax = gca;

ax.XScale = 'log';
ax.YScale = 'log';

ax.FontSize = 12;

% leg = legend(  {'s=8',  's=32', 's=96','s=128','s=160','s=256'},'Location' ,'north', ...
%     'FontSize' , legFtSize , 'Box' , 'off'  ,'NumColumns' ,  1 );



ax.XTick = 10.^[-12 0 +12];
ax.YTick = 10.^[-30  +6];

% ax.XLim = 10.^[-16 16];

xline(10^4 , 'LineStyl'  , '--')
yline(0 , 'LineStyl'  , '-')

% leg = legend(  {'r=4',  'r=16', 'r=64','r=256'},'Location' ,'northeast', ...
%     'FontSize' , legFtSize , 'Box' , 'off'  ,'NumColumns' ,  1 );


axis square
