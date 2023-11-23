%% plot systems controllability 9nets
clear all
%% add path
subdir = genpath("fig_utils");
addpath(subdir);
subdir = genpath("results_data");
addpath(subdir);
%% load the tensor of system level controllability for all subjects
load('results_data/controllabilityTensor_9nets.mat')

rowSum = sum(controlTensor , 2);
okRows = rowSum ~= -Inf;
% avgcontrolSystems = mean(10.^controlTensor(okRows,:));
% stdcontrolSystems = std(10.^controlTensor(okRows,:));

avgcontrolSystems = mean(controlTensor(okRows,:));
stdcontrolSystems = std(controlTensor(okRows,:));

nNetworks = 9;
% load colormap schaefer
load('mySchaeferColorMap_9nets.mat')

%% %% metagraph only cortical

avgcontrolSystemsMat = reshape(avgcontrolSystems , [nNetworks , nNetworks]);

Networks9CellVert = {'Vis'; 'SomMot'; 'DorsAttn';  'SalVentAttn';...
    'Limbic'; 'Cont'; 'Default'; 'TempPar' ; 'Sub'};

Networks9CellVert = {'VIS'; 'SMN'; 'DAN';  'SVAN';...
    'LIM'; 'FPCN'; 'DMN'; 'TPJ' ; 'SUB'};
Networks9CellVertCort = {'VIS'; 'SMN'; 'DAN';  'SVAN';...
    'LIM'; 'FPCN'; 'DMN'; 'TPJ' };

%%%% remove diagonal
avgcontrolSystemsMat2plot = avgcontrolSystemsMat;
for k = 1 : nNetworks
    avgcontrolSystemsMat2plot(k,k) = nan;
    avgcontrolSystemsMat2plot(end,k) = nan;
    avgcontrolSystemsMat2plot(k,end) = nan;
end

fig = figure;
fig.Position(3:4) = 390 * [1 1];

img = imagesc((avgcontrolSystemsMat2plot(1:end-1 , 1:end-1)).^0.5);

ax = gca;
ax.XTick = 1:nNetworks-1;
ax.XTickLabel = Networks9CellVertCort;
ax.YTick = 1:nNetworks-1;
ax.YTickLabel = Networks9CellVertCort;

% ax.FontWeight = 'bold';
xlabel('driver' , 'FontSize', 14);
ylabel('target', 'FontSize', 14);
ax.XAxisLocation = 'top';
% xlabel('Controllability macro')

% colormap('jet')
% colormap(flipud((hot)))
myHot = hot();
myCmap = flipud((hot));
% myCmap = [ones(64,3);flipud((hot))];

% myCmap = createcolormap([1 1 1] , myHot(0 , :));
myCord = colororder;
myCmap = createcolormap([1 1 1] ,[1 1 1] , myCord(2,:) , rgb('FireBrick'));


myCmap = createcolormap([1 1 1] , myCord(2,:) , rgb('FireBrick'));

% %%%%% enhance contrast
% newmap = contrast(abs(log10(avgcontrolSystemsMat2plot(1:end-1 , 1:end-1))));
%
% myCmap = flipud(newmap)+ rgb('FireBrick');
% myCmap = (myCmap - min(min(myCmap))) / (max(max(myCmap)) - min(min(myCmap)));
%


colormap("jet")

colormap(flipud((pink())))

cb = colorbar;
cb.Label.String  = {'log_1_0(\lambda_m_i_n(W)) '};
cb.Label.FontSize = 16;

% cb.Position(1) = 1;
ax.FontSize = 16;

axis square

% img.AlphaData = 0.1;

img.CDataMapping = 'scaled';


ax.XGrid = 'on';
ax.YGrid = 'on';
ax.XTick = (1:nNetworks-1)-0.5;
ax.YTick = (1:nNetworks-1)-0.5;
% ax.GridLineWidth = 1;
ax.GridColor = 0.05*[1 1 1];

ax.GridAlpha = 0.6;

%% plot network circle

avgcontrolSystemsMatWithZeros = (reshape(avgcontrolSystems , [nNetworks , nNetworks]));
avgcontrolSystemsMatWithZeros = avgcontrolSystemsMatWithZeros - min(min(avgcontrolSystemsMatWithZeros));
avgcontrolSystemsMatWithZeros = avgcontrolSystemsMatWithZeros - diag(diag(avgcontrolSystemsMatWithZeros));

avgcontrolSystemsMatWithZeros = avgcontrolSystemsMatWithZeros/ max(max(avgcontrolSystemsMatWithZeros));

G = digraph(avgcontrolSystemsMatWithZeros(1:8 , 1:8)');

Networks9CellCort = {' Vis', ' SomMot',  'DorsAttn ',  'SalVentAttn ',...
    'Limbic ', 'Cont ',' Default', ' TempPar' };

Networks9CellCort =  {'', '',  '',  '','', '','', '' };


myEdgeCdata = G.Edges.Weight;
myEdgeCdata = (myEdgeCdata - min(myEdgeCdata))/( max(myEdgeCdata) -  min(myEdgeCdata));
myGrayMap = createcolormap(rgb('White') , rgb('Black')  , rgb('Black') );

% myGrayMap = flipud(gray());
%
% myGrayMap = myGrayMap(1:end-100,:);
%
% myGrayMap = 0.5*ones(200,3);


myEdgeColor =   myGrayMap(1 + floor((myEdgeCdata.^2)*(size(myGrayMap , 1)-1)) , :); % scale

outdeg = sum(avgcontrolSystemsMatWithZeros(1:8 , 1:8));
outdeg  = (outdeg - min(outdeg))/( max(outdeg) -  min(outdeg));

% myEdgesWidthData = G.Edges.Weight;
% myEdgeCdata = (myEdgeCdata - min(myEdgeCdata))/( max(myEdgeCdata) -  min(myEdgeCdata));

fig = figure;
% fig
pGcircle = plot(G , 'NodeLabel',Networks9CellCort , 'NodeColor',myColorMap(1:8,:) , 'LineWidth',0.1 + 10*myEdgeCdata.^2 ...
    , 'EdgeColor', myEdgeColor , 'MarkerSize', 10 + 10*outdeg , 'NodeFontSize', 17 , 'Layout','circle' ...
    ,'ArrowSize', 1 + 15*myEdgeCdata.^1.5 ,'EdgeAlpha',0.5 );


maxCortWeight = max(G.Edges.Weight);
ax = gca;
ax.Box = 'off';
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
ax.FontSize = 14;



reorderNodes = 1;
if reorderNodes

    %%% swap net1 and net2
    net1 = 5;
    net2 = 7;

    bb = pGcircle.XData(net2);
    pGcircle.XData(net2) = pGcircle.XData(net1) ;
    pGcircle.XData(net1)  = bb;
    %%%%%%%%%
    %%%%%%%%%
    bb = pGcircle.YData(net2);
    pGcircle.YData(net2) = pGcircle.YData(net1) ;
    pGcircle.YData(net1)  = bb;
    %%%%%%%%%
    
     %%% swap net1 and net2
    net1 = 2;
    net2 = 4;

    bb = pGcircle.XData(net2);
    pGcircle.XData(net2) = pGcircle.XData(net1) ;
    pGcircle.XData(net1)  = bb;
    %%%%%%%%%
    %%%%%%%%%
    bb = pGcircle.YData(net2);
    pGcircle.YData(net2) = pGcircle.YData(net1) ;
    pGcircle.YData(net1)  = bb;
    %%%%%%%%%
    
         %%% swap net1 and net2
    net1 = 5;
    net2 = 6;

    bb = pGcircle.XData(net2);
    pGcircle.XData(net2) = pGcircle.XData(net1) ;
    pGcircle.XData(net1)  = bb;
    %%%%%%%%%
    %%%%%%%%%
    bb = pGcircle.YData(net2);
    pGcircle.YData(net2) = pGcircle.YData(net1) ;
    pGcircle.YData(net1)  = bb;
    %%%%%%%%%

             %%% swap net1 and net2
    net1 = 5;
    net2 = 7;

    bb = pGcircle.XData(net2);
    pGcircle.XData(net2) = pGcircle.XData(net1) ;
    pGcircle.XData(net1)  = bb;
    %%%%%%%%%
    %%%%%%%%%
    bb = pGcircle.YData(net2);
    pGcircle.YData(net2) = pGcircle.YData(net1) ;
    pGcircle.YData(net1)  = bb;
    %%%%%%%%%

                 %%% swap net1 and net2
    net1 = 5;
    net2 = 8;

    bb = pGcircle.XData(net2);
    pGcircle.XData(net2) = pGcircle.XData(net1) ;
    pGcircle.XData(net1)  = bb;
    %%%%%%%%%
    %%%%%%%%%
    bb = pGcircle.YData(net2);
    pGcircle.YData(net2) = pGcircle.YData(net1) ;
    pGcircle.YData(net1)  = bb;
    %%%%%%%%%
    %%%%%%%%%
    % bb = pGcircle.XData(8);
    % pGcircle.XData(8) = pGcircle.XData(5) ;
    % pGcircle.XData(5)  = bb;
    % %%%%%%%%%
    % %%%%%%%%%
    % bb = pGcircle.YData(8);
    % pGcircle.YData(8) = pGcircle.YData(5) ;
    % pGcircle.YData(5)  = bb;
    % %%%%%%%%%
    %
    % %%%%%%%%%
    % bb = pGcircle.XData(7);
    % pGcircle.XData(7) = pGcircle.XData(6) ;
    % pGcircle.XData(6)  = bb;
    % %%%%%%%%%
    % %%%%%%%%%
    % bb = pGcircle.YData(7);
    % pGcircle.YData(7) = pGcircle.YData(6) ;
    % pGcircle.YData(6)  = bb;
    % %%%%%%%%%
end
