%% script to look dependance of controllability measures with strength

load('results_data/bigTensor_results_fullBrainControl_on_6134_subjects.mat')

nSubGood = sum(isOkArray);
goodSubjects = find(isOkArray); % without Inf or nan

nSubjects = 6134;


%% loop over dimensions and subjects

MarkerFaceAlpha2 = 0.1;

ROI_info_Table = readtable( 'results_data/Schaefer200_allinfo.csv');
networksAssignment = table2array(ROI_info_Table(:,5));
nNetworks = max(networksAssignment);

% load('/Users/remy.benmessaoud/ownCloud/remy.benmessaoud/MATLAB/UK-biobank/controllability_ukb/colormap_schaefer_8nets.mat')
% myColorMap(9,:) = rgb('SaddleBrown');

load('mySchaeferColorMap_9nets.mat')

nSubGood2 = 5000;
nSubGood2 = 6123;

XX = [reshape(strengthTensor(:,goodSubjects(1:nSubGood2)) , [214 * nSubGood2 1 ])...
    reshape(averageControllabilityFullBrainTensor(:,goodSubjects(1:nSubGood2)) , [214 * nSubGood2 1 ])...
    reshape(allGramMinBasedTensor(:,goodSubjects(1:nSubGood2)) , [214 * nSubGood2 1 ])...
    reshape(lapGramMinBasedTensor(:,5,goodSubjects(1:nSubGood2)) , [214 * nSubGood2 1 ])];

[rho,pval] = corr(XX ,'Rows' , 'complete')


mkSz = 4;

fig=figure;
fig.Position(3:4) = 250 * [1 1];


x = 1 + reshape(strengthTensor(:,goodSubjects(1:nSubGood2)) , [1 214 * nSubGood2]) ;


y = reshape(10.^lapGramMinBasedTensor(:,5,goodSubjects(1:nSubGood2)) , [1 214 * nSubGood2]);
myYLim = 10.^[-16 -4];
myXLim = 10.^[2 5];

s = scatter(x, y, ...
    6*ones(1 , 214 * nSubGood2) , reshape(networksAssignment*ones(1,nSubGood2) , [1 214 * nSubGood2]) ...
    , 'filled' ,'MarkerFaceAlpha' , MarkerFaceAlpha2  , 'MarkerFaceColor','flat');

%     b = regress(y' ,[ones(214 * nSubGood2,1)  x']);
%     h = plot([min(x) max(x)] , b(1) + b(2)*[min(x) max(x)]);


ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
ax.YLim = myYLim ;
ax.XLim = myXLim;

ax.FontSize = 14;

axis square
%     h = line([min(x) max(x)] , b(1) + b(2)*[min(x) max(x)]);
%     h.Color = 'k';
%     h.LineWidth = 2;

colormap(myColorMap)

ax.XTick = 10.^[2 5];

ax.YTick = 10.^[-14 -6];


