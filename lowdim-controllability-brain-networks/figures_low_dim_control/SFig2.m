
clear all;

%% add path
subdir = genpath("fig_utils");
addpath(subdir);
subdir = genpath("results_data");
addpath(subdir);

%% plot gramian spectrum

data_dir = 'results_data/gramian_full_brain_vary_r';

ROI_info_Table = readtable( 'results_data/Schaefer200_allinfo.csv');
networksAssignment = table2array(ROI_info_Table(:,5));
nNetworks = max(networksAssignment);


load('results_data/1000367_gram_vary_r.mat')




%% per system
mkSize = 5;
% cOrd = colororder;

load('mySchaeferColorMap_9nets.mat')

% maxdim = 50;
nSystems = 9;
n =214;

lW = 0.1;
fig = figure;
fig.Position(3:4)  = [600 400];

linearray = [];


for kNet = flip(1:   nSystems)

    xTick = (1:maxdim)+kNet/nSystems*0.8 ;%+ 0.1*rand(1,maxdim);

    xTick = (0.6:maxdim-0.4)+ kNet/nSystems*0.8 ;%+ 0.1*rand(1,maxdim);
    %xTick = (1:maxdim)+ kNet/nSystems*0.8 ;%+ 0.1*rand(1,maxdim);


    % lambdaArray = reshape(lambdaMinGramOut(:,kSys,:) , [n maxdim]);
    lambdaArray = reshape(lambdaMinGramOut(:,kNet,:) , [n,maxdim]);

    
    %%%%%%% plot negative values
    poslambdaArray = nan(size(lambdaArray));
    %     neglambdaArray = nan(size(lambdaArray));
    poslambdaArray(lambdaArray > 0) = lambdaArray(lambdaArray > 0);
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% plot negative values
    neglambdaArray = nan(size(lambdaArray));
    %     neglambdaArray = nan(size(lambdaArray));
    neglambdaArray(lambdaArray < 0) = lambdaArray(lambdaArray < 0);

    % axes(ha(1))

    subplot(2,1,1)

    hold on;

    % pp = plot( xTick , poslambdaArray(nodes,flip(1:length(xTick)))',...
    %     'Color', myColorMap(kNet , :), 'LineWidth' , lW...
    %     , 'LineStyle','none' , 'Marker','o' , 'MarkerSize',mkSize);

    pp = plot( xTick , poslambdaArray,...
        'Color', myColorMap(kNet , :), 'LineWidth' , lW...
        , 'LineStyle','none' , 'Marker','.' , 'MarkerSize',mkSize);

    % pp = plot( xTick , lambdaArray,...
    %      'LineWidth' , lW...
    %     , 'LineStyle','none' , 'Marker','o' , 'MarkerSize',mkSize);

    % linearray = [linearray pp];

    xlim([0 55])
    hold on;



    minAll = min(min(abs(lambdaArray)));
    %     plot(log10(minVal+lambdaArray) - log10(minVal+lambdaArray) )


    % axes(ha(2))
    subplot(2,1,2)

    hold on;

    % pNeg = plot(xTick , negLambda2show(nodes,flip(1:length(xTick))) ...
    %     , 'Color', myColorMap(kNet , :), 'LineWidth' , lW...
    %     , 'LineStyle','none' , 'Marker','o', 'MarkerSize',mkSize+1);
    % hold on;

    pNeg = plot(xTick , neglambdaArray ...
        , 'Color', myColorMap(kNet , :), 'LineWidth' , lW...
        , 'LineStyle','none' , 'Marker','.', 'MarkerSize',mkSize);
    xlim([0 55])
    hold on;



    % pNeg = plot(xTick , negLambda2show ...
    %     ,  'LineWidth' , lW...
    %     , 'LineStyle','none' , 'Marker','o', 'MarkerSize',mkSize);

end
% end

subplot(2,1,1)
ax = gca;
ax.YScale = 'log';
ax.FontSize = 14;

% ax.YTick = [10^-20 10^-10];
ax.XTick = [];
xline(5,'--')
ax.Box = "off";


ax.XScale = 'log';
ax.XDir = 'reverse';
ax.XLim = [0.65 maxdim+1];

ax1 = ax;


%%%%%% flip axis 2
subplot(2,1,2)
ax = gca;
ax.YScale = 'log';
hold on;
ax.FontSize = 14;

% ax.YTick = [10^-25 10^-20];
ax.YLim = [-10^-19.5 -10^-28.5];
ax.YTick = [-10^-21 -10^-28];

ax.XScale = 'log';
ax.XDir = 'reverse';
% ax.YDir = 'reverse';


ax.XTick = 5:5:45;
ax. Box = "off";

% ax.XTickLabel = string(flip((ax.XTickLabel)));
xline(5,'--') % 45 = max dim - r = 50 - 5;
% ax.XLim = [0 51];


%%% put subplots together
ax.Position(2) = 0.225;

% ax.XScale = 'log';
ax2 = ax;

%ax1.XLim = [1 55];
%ax2.XLim = [1 55];

ax2.XTick =[1 2 3  5 7 10 20 37];

ax.XLim = [0.65 maxdim+1];

ax.XTickLabel{13,1} = 'n';


%% Inset : load r after which negztive values appear
load('results_data/r_critic_9targetNets.mat')



r_mean = mean(r_critic)

r_median = median(r_critic)

r_std = std(r_critic)

%% plot histo inset
% sizePic = 0.225;
% % ax2 = axes('Position',[0.205 0.73 (sizePic- 0.05)  sizePic]);
% ax2 = axes('Position',[0.3 0.6 (sizePic)  sizePic]);
%
% ax = ax2;

fig = figure;
fig.Position(3:4) = 160 * [1 1];

histogram(r_critic , 'FaceColor',0.25*[1 1 1]);
ax = gca;

%xline(5)
ax.YScale = "log";
ax.FontSize = 20;
ax.XTick = 3:2:9;
% ax.YTick = [10^1 10^3];
ax.Box = 'off';
axis tight square
ax.XLim(1) = 2;
ax.XLim(2) = 11;
ax.XTickLabelRotation = 0;


%% fig check

% figure;
% plot(mean(lambdaArray(:,8:end)))
%
