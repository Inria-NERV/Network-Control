
clear all

subdir = genpath("fig_utils");
addpath(subdir);
subdir = genpath("results_data");
addpath(subdir);

%% load tensor 
load('results_data/tensor_LeftRight.mat')

LLtensor = mean(Larray(:,1:100),2);
RLtensor = mean(Larray(:,101:200),2);

RRtensor = mean(Rarray(:,101:200),2);
LRtensor = mean(Rarray(:,1:100),2);


LL = mean(LLtensor)
RL = mean(RLtensor)
RR = mean(RRtensor)
LR = mean(LRtensor)


ROI_info_Table = readtable( 'results_data/Schaefer200_allinfo.csv');
networksAssignment = table2array(ROI_info_Table(:,5));
nNetworks = max(networksAssignment);

leftIdx = [1:100 201:207];
rightIdx = [101:200 208:214];

intraCohenHem = computeCohen_d( LLtensor,RRtensor, 'paired');
interCohenHem = computeCohen_d( LRtensor,RLtensor, 'paired');


%% satsts sing test
figure; histogram(log(abs(RLtensor)) )
hold on; histogram(log(abs(LRtensor)) );
hold on; histogram(RLtensor-LRtensor )

[p,h,stats] = signtest(LLtensor , RRtensor)
[p,h,stats] = signtest(RLtensor , LRtensor)


[h,p,ci,stats] = ttest((LLtensor) , (RRtensor))
[h,p,ci,stats] = ttest((RLtensor) , (LRtensor))


%% chech stat in the same hemisphere 

left_ipsi_vs_contra = computeCohen_d( LLtensor,LRtensor, 'paired')
right_ipsi_vs_contra = computeCohen_d( RRtensor,RLtensor, 'paired')

[p,h,stats] = signtest(LLtensor , LRtensor)

[p,h,stats] = signtest(RRtensor , RLtensor)

%% plot bar no stack but group
LLstd = std(LLtensor);
RLstd = std(RLtensor);
RRstd = std(RRtensor);
LRstd = std(LRtensor);

data = ([LL LR ; RR RL]);
datastd = ([LLstd RRstd ; LRstd RLstd]);

figure;
bb0 = bar(ones(10,4) );

bb0(1).FaceColor = 'flat';
bb0(2).FaceColor = 'flat';
bb0(3).FaceColor = 'flat';
bb0(4).FaceColor = 'flat';

bb0(1).CData = ones(10,1)* rgb('DarkRed');
bb0(2).CData = ones(10,1)* rgb('Pink');
bb0(3).CData = ones(10,1)* rgb('DarkBlue');
bb0(4).CData = ones(10,1)* rgb('LightBlue');


fig = figure;
fig.Position(3:4) = [270 270];

bb = bar( data );
bb(1).BarWidth = 0.5;
bb(2).BarWidth = 0.5;

bb(1).FaceColor = 'flat';
bb(2).FaceColor = 'flat';


ax = gca;


ylim( [5*10^-14 1*10^-11])
xlim( [0.6 2.4])

bb(1).CData = [rgb('DarkRed') ; rgb('DarkBlue')];

bb(2).CData = [rgb('Pink') ; rgb('LightBlue')];


ax.FontSize = 14;


ax.XTickLabel = {'Left' , 'Right'};


hold on;
errorbar(bb(1).XEndPoints , bb(1).YEndPoints ,nan(1,2), [LLstd RRstd] , 'LineStyle', 'none', 'Color', 'k','LineWidth',1.75)

hold on;
errorbar(bb(2).XEndPoints , bb(2).YEndPoints ,nan(1,2), [LRstd RLstd] , 'LineStyle', 'none', 'Color', 'k','LineWidth',1.75)


ax.YScale = 'log';

legend(ax , [bb0(1) , bb0(2) , bb0(3) , bb0(4)] , {'L→L' , 'L→R' , 'R→R' , 'R→L'} ,'Box' , 'off', ...
    'FontSize' , 14 , 'NumColumns', 2)

axis square
ax.Box = 'off';

ax.YGrid = 'on';

%% plot bar no stack but group intra inter
LLstd = std(LLtensor);
RLstd = std(RLtensor);
RRstd = std(RRtensor);
LRstd = std(LRtensor);

data = ([LL RR ; LR RL]);
datastd = ([LLstd LRstd ; RRstd RLstd]);

figure;
bb0 = bar(ones(10,4) );

bb0(1).FaceColor = 'flat';
bb0(2).FaceColor = 'flat';
bb0(3).FaceColor = 'flat';
bb0(4).FaceColor = 'flat';

bb0(1).CData = ones(10,1)* rgb('DarkRed');
bb0(2).CData = ones(10,1)* rgb('Pink');
bb0(3).CData = ones(10,1)* rgb('DarkBlue');
bb0(4).CData = ones(10,1)* rgb('LightBlue');


fig = figure;
fig.Position(3:4) = [270 270];

bb = bar( data );
bb(1).BarWidth = 0.66;
bb(2).BarWidth = 0.66;

bb(1).FaceColor = 'flat';
bb(2).FaceColor = 'flat';


ax = gca;


ylim( [10^-14 10^-10])
xlim( [0.6 2.4])

bb(1).CData = [rgb('DarkRed') ; rgb('Pink')];

bb(2).CData = [rgb('DarkBlue') ; rgb('LightBlue')];


ax.FontSize = 14;


ax.XTickLabel = {'Ipsi' , 'Contra'};


hold on;
errorbar(bb(1).XEndPoints , bb(1).YEndPoints ,nan(1,2), [LLstd LRstd] , 'LineStyle', 'none', 'Color', 'k','LineWidth',1)

hold on;
errorbar(bb(2).XEndPoints , bb(2).YEndPoints ,nan(1,2), [RRstd RLstd] , 'LineStyle', 'none', 'Color', 'k','LineWidth',1)


ax.YScale = 'log';

legend(ax , [bb0(1) , bb0(2) , bb0(3) , bb0(4)] , {'L→L' , 'L→R' , 'R→R' , 'R→L'} ,'Box' , 'off', ...
    'FontSize' , 16 , 'NumColumns', 2)

axis square
ax.Box = 'off';

ax.YGrid = 'off';


