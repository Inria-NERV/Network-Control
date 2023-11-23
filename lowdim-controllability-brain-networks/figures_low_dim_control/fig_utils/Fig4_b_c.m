
clear all;


%% add path 
subdir = genpath("fig_utils");
addpath(subdir);
subdir = genpath("results_data");
addpath(subdir);
%%

ROI_info_Table = readtable( '/results_data/Schaefer200_allinfo.csv');
networksAssignment = table2array(ROI_info_Table(:,5));
nNetworks = max(networksAssignment);

n=214;



leftIdx = [1:100 201:207];
rightIdx = [101:200 208:214];

%% load results data

% load('/Users/remy.benmessaoud/ownCloud/remy.benmessaoud/MATLAB/UK-biobank/controllability_ukb/results_node_level/bigTensor_results_target9nets_on_6134_subjects.mat')

load('results_data/bigTensor_results_target9nets_on_6134_subjects.mat')


%% ratio in out
metricInside = zeros(nNetworks, nSubjects);
metricOutside = zeros(nNetworks, nSubjects);

lookBestIn = zeros(n,nNetworks);
lookBestOut = zeros(n,nNetworks);
metricAll = zeros(nNetworks, nSubjects);


for netIdx = 1:9
    target = find(networksAssignment == netIdx);
    notTarget = setdiff(1:n , target);

    allSubjectsMetric = (abs(reshape(lapGramMinBasedTensor(:,netIdx, :), [n , nSubjects])));

    allSubjectsMetric(allSubjectsMetric==-Inf) = nan;

    metricAll(netIdx,:)=mean(allSubjectsMetric);

    metricInside(netIdx,:) = mean(allSubjectsMetric(target,:));
    metricOutside(netIdx,:) = mean(allSubjectsMetric(notTarget,:));

    lookBestIn(target,netIdx) = mean(allSubjectsMetric(target,:),2);
    lookBestOut(notTarget,netIdx) = mean(allSubjectsMetric(notTarget,:),2);
    
end

avgIn = mean(metricInside,2)';
avgOut = mean(metricOutside,2)';

stdIn = std(metricInside,0,2)';
stdOut = std(metricOutside,0,2)';

%%%%%%%%%%%%%%%%%%%%%
metricRatio = log10(metricInside./metricOutside);
% metricRatio = (metricInside-metricOutside);


% metricRatio = (metricInside-metricOutside)./(metricOutside);

%%%%% get rid of outliers
metricRatio(metricRatio==0) = nan;
metricRatio(metricRatio==Inf) = nan;
metricRatio(metricRatio==-Inf) = nan;
metricRatio(metricRatio==-Inf) = nan;

% metricRatio(metricRatio>10^6) = nan;


avgRatio = mean(metricRatio,2 , 'omitnan')';
stdRatio = std(metricRatio,0,2, 'omitnan')';


%%%%%%%%%
figure;histogram(metricRatio(1,:)) ; ax = gca;ax.XScale = 'log';
%%%%%%%%%

avgIn./avgOut

avgIn-avgOut



lookBestOut(201:end , :) = 0;

[ii , bb ]=max(lookBestOut)


%%%%%% ANOVA all
[p,tbl,stats] = anova1(metricAll');
% [p,tbl,stats] = kruskalwallis(metricAll');

[c,m,h,gnames] = multcompare(stats);


gnames = {'VIS'; 'SMN'; 'DAN';  'SVAN';...
                'LIM'; 'FPCN'; 'DMN'; 'TPJ' ; 'SUB'};

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl.("Group A") = gnames(tbl.("Group A"));
tbl.("Group B") = gnames(tbl.("Group B"))


%%%%%% ANOVA Ratio
% [p,tbl,stats] = anova1(metricRatio');
[p,tbl,stats] = kruskalwallis(metricRatio');
[c,m,h,gnames] = multcompare(stats);


gnames = {'VIS'; 'SMN'; 'DAN';  'SVAN';...
                'LIM'; 'FPCN'; 'DMN'; 'TPJ' ; 'SUB'};

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl.("Group A") = gnames(tbl.("Group A"));
tbl.("Group B") = gnames(tbl.("Group B"))


%% plot bars all drivers
barData = zeros(1 , nNetworks);
barDataErr = zeros(1 , nNetworks);

for netIdx = 1:9
    target = find(networksAssignment == netIdx);
    notTarget = setdiff(1:n , target);

    allSubjectsMetric = reshape(lapGramMinBasedTensor(:,netIdx, :), [n , nSubjects]);

    myMetric = mean(allSubjectsMetric, 2);
    barData(netIdx) = mean(myMetric);
    barDataErr(netIdx) = std(mean(allSubjectsMetric, 1));
end

[ii , kk] = sort(barData, 'descend');

load('mySchaeferColorMap_9nets.mat')
fig = figure;
fig.Position(3:4)=230*[1 0.5];

bb = bar(barData(kk));
hold on;
err = errorbar(1:nNetworks, bb.YEndPoints , nan(size(barDataErr(kk) )),barDataErr(kk) , ...
        'Color', 'k' , 'LineStyle' , 'none');

bb.FaceColor = 'flat';
ax = gca;
ax.YScale = 'log';

bb.CData = myColorMap(kk,:);
netsNames  = {'VIS'; 'SMN'; 'DAN';  'SVAN';...
                'LIM'; 'FPCN'; 'DMN'; 'TPJ' ; 'SUB'};
netsNamesOrd = {};
for kNet = 1:nNetworks
    netsNamesOrd{1,kNet} = netsNames{ kk(kNet),1};
end
ax.XTickLabel =netsNamesOrd;

% ax.YTick = 10.^[-10:-8];
ax.FontSize = 14;
% ax.YGrid = 'on';
ax.Box = "off";
% axis square
ax.XLim = [  0.3   9.7];


%% plot bars Ratio Inside /outside
barData =  avgRatio ;
barDataErr = stdRatio;


[ii , kk] = sort(barData, 'descend');

load('mySchaeferColorMap_9nets.mat')
fig = figure;
fig.Position(3:4)=230*[1 0.5];


bb = bar(barData(kk));
hold on;
err = errorbar(1:nNetworks, bb.YEndPoints , nan(size(barDataErr(kk) )),barDataErr(kk) , ...
        'Color', 'k' , 'LineStyle' , 'none');

bb.FaceColor = 'flat';
ax = gca;
% ax.YScale = 'log';

bb.CData = myColorMap(kk,:);
netsNames  = {'VIS'; 'SMN'; 'DAN';  'SVAN';...
                'LIM'; 'FPCN'; 'DMN'; 'TPJ' ; 'SUB'};
netsNamesOrd = {};
for kNet = 1:nNetworks
    netsNamesOrd{1,kNet} = netsNames{ kk(kNet),1};
end
ax.XTickLabel =netsNamesOrd;

% ax.YTick = 10.^[-10:-8];
ax.FontSize = 14;
% ax.YGrid = 'on';
ax.Box = "off";
% axis square
ax.XLim = [  0.3   9.7];

% ax.YScale = 'log';


bb.BaseValue =  0.75;
% ax.YLim(1) = bb.BaseValue;
% ax.XTickLabelRotation = 38;

