
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
% figure; histogram(log(abs(RLtensor)) )
% hold on; histogram(log(abs(LRtensor)) );
% hold on; histogram(RLtensor-LRtensor )

[p,h,stats] = signtest(LLtensor , RRtensor)

[p,h,stats] = signtest(RLtensor , LRtensor)



%% do the same but per system
nNetworks = 9;
N = 6134;

LLsysArray = zeros(1,nNetworks);
RRsysArray = zeros(1,nNetworks);
LRsysArray = zeros(1,nNetworks);
RLsysArray = zeros(1,nNetworks);

LLsysArrayStd = zeros(1,nNetworks);
RRsysArrayStd = zeros(1,nNetworks);
LRsysArrayStd = zeros(1,nNetworks);
RLsysArrayStd = zeros(1,nNetworks);

intraPval = zeros(1,nNetworks);
interPval = zeros(1,nNetworks);

intraCohen = zeros(1,nNetworks);
interCohen = zeros(1,nNetworks);

intraLI = zeros(N,nNetworks);
interLI = zeros(N,nNetworks);

for kNet = 1:nNetworks

nodes = find(networksAssignment ==kNet);

LLsysArray(kNet) = mean(mean(Larray(:,setdiff(nodes , rightIdx)),2));
RRsysArray(kNet) = mean(mean(Rarray(:,setdiff(nodes , leftIdx)),2));

LRsysArray(kNet) = mean(mean(Rarray(:,setdiff(nodes , rightIdx)),2));
RLsysArray(kNet) = mean(mean(Larray(:,setdiff(nodes , leftIdx)),2));
%%%%%%%%%%%%%%
LLsysArrayStd(kNet) = std(mean(Larray(:,setdiff(nodes , rightIdx)),2));
RRsysArrayStd(kNet) = std(mean(Rarray(:,setdiff(nodes , leftIdx)),2));

LRsysArrayStd(kNet) = std(mean(Rarray(:,setdiff(nodes , rightIdx)),2));
RLsysArrayStd(kNet) = std(mean(Larray(:,setdiff(nodes , leftIdx)),2));


[p,h,stats] = signtest(log(abs(mean(Larray(:,setdiff(nodes , rightIdx)),2))) , log(abs(mean(Rarray(:,setdiff(nodes , leftIdx)),2))) );
intraPval(kNet)=p;

[p,h,stats] = signtest(log(abs(mean(Larray(:,setdiff(nodes , leftIdx)),2))) , log(abs(mean(Rarray(:,setdiff(nodes , rightIdx)),2))));
interPval(kNet)=p;


% [h,p,stats] = ttest(log(abs(mean(Larray(:,setdiff(nodes , rightIdx)),2))+eps) - log(abs(mean(Rarray(:,setdiff(nodes , leftIdx)),2))+eps));
% intraPval(kNet)=p;
% 
% [h,p,stats] = ttest(log(abs(mean(Larray(:,setdiff(nodes , leftIdx)),2))+eps) - log(abs(mean(Rarray(:,setdiff(nodes , rightIdx)),2))+eps));
% interPval(kNet)=p;

% [p,h,stats] = ranksum(mean(Larray(:,setdiff(nodes , rightIdx)),2) , mean(Rarray(:,setdiff(nodes , leftIdx)),2));
% intraPval(kNet)=p;
% 
% [p,h,stats] = ranksum(mean(Larray(:,setdiff(nodes , leftIdx)),2) , mean(Rarray(:,setdiff(nodes , rightIdx)),2));
% interPval(kNet)=p;

intraCohen(kNet) = computeCohen_d( mean(Rarray(:,setdiff(nodes , leftIdx)),2),mean(Larray(:,setdiff(nodes , rightIdx)),2), 'paired');
interCohen(kNet) = computeCohen_d(mean(Larray(:,setdiff(nodes , leftIdx)),2) , mean(Rarray(:,setdiff(nodes , rightIdx)),2), 'paired');


intraLI(:,kNet) = (mean(Rarray(:,setdiff(nodes , leftIdx)),2)-mean(Larray(:,setdiff(nodes , rightIdx)),2))...
    ./(mean(Rarray(:,setdiff(nodes , leftIdx)),2)+mean(Larray(:,setdiff(nodes , rightIdx)),2));
interLI(:,kNet) = (mean(Larray(:,setdiff(nodes , leftIdx)),2) - mean(Rarray(:,setdiff(nodes , rightIdx)),2))...
    ./(mean(Larray(:,setdiff(nodes , leftIdx)),2) + mean(Rarray(:,setdiff(nodes , rightIdx)),2));


end

barData = [LLsysArray ;RRsysArray; LRsysArray; RLsysArray];

fig = figure;
fig.Position(3:4) = [530 270];

bb = bar( barData' );
myColors = [rgb('DarkRed'); rgb('DarkBlue'); rgb('Pink');rgb('LightBlue')];

for kBar = 1:4
bb(kBar).BarWidth = .9;

bb(kBar).FaceColor = myColors(kBar , :);

end

ax = gca;
ax.YScale = 'log';
ax.Box = 'off';
 ylim([10^-16 10^-10]);
% xlim( [0.6 2.4])
ax.FontSize = 14;
ax.XTickLabel = {'VIS' , 'SMN', 'DAN', 'SVAN', 'LIM', 'FPCN', 'DMN', 'TPJ','SUB'};
%ylim([10^-14 8*10^-11])

ax.YGrid = 'on';ax.XGrid = 'on';

stdData = [LLsysArrayStd;RRsysArrayStd;LRsysArrayStd;RLsysArrayStd];
for kBar = 1:4
hold on;
errorbar(bb(kBar).XEndPoints , bb(kBar).YEndPoints ,nan(size(LLsysArrayStd)), stdData(kBar,:) , 'LineStyle', 'none', 'Color', 'k')
end

grid off
%set(gca,('xticklabelrotation',45)

%% 
kNet = 5;
nodes = find(networksAssignment ==kNet);

N = 6134;

figure; histogram(mean(Larray(:,setdiff(nodes , rightIdx)),2) )
hold on; histogram(mean(Rarray(:,setdiff(nodes , leftIdx)),2)) ;


figure; histogram(mean(Larray(:,setdiff(nodes , rightIdx)),2) -mean(Rarray(:,setdiff(nodes , leftIdx)),2))

[p,h,stats] = ranksum(mean(Larray(1:N,setdiff(nodes , rightIdx)),2) , mean(Rarray(1:N,setdiff(nodes , leftIdx)),2))
% [h,p,stats] = ttest(mean(Larray(:,setdiff(nodes , rightIdx)),2) , mean(Rarray(:,setdiff(nodes , leftIdx)),2))


%% plot cohen's d

fig = figure;
fig.Position(3:4) = [300 300];

bb = bar(  [intraCohen  ; interCohen]');
myColors = [rgb('DarkRed'); rgb('DarkBlue'); rgb('Pink');rgb('LightBlue')];

for kBar = 1:2
bb(kBar).BarWidth = 0.45;

bb(kBar).FaceColor = myColors(kBar , :);

end

ax = gca;
% ax.YScale = 'log';
ax.Box = 'off';
% ylim(10^-12 * [5 13])
% xlim( [0.6 2.4])
ax.FontSize = 14;
ax.XTickLabel = {'VIS' , 'SMN', 'DAN', 'SVAN', 'LIM', 'FPCN', 'DMN', 'TPJ','SUB'};
% ylim([10^-14 8*10^-11])

ax.YGrid = 'on';ax.XGrid = 'on';


%% plot cohen's d

fig = figure;
fig.Position(3:4) = [300 500];

bb = barh(  [intraCohen  ; interCohen]');
myColors = [rgb('DarkRed'); rgb('DarkBlue'); rgb('Pink');rgb('LightBlue')];
myColors = [rgb('Black'); rgb('Gray'); rgb('Pink');rgb('LightBlue')];

for kBar = 1:2
bb(kBar).BarWidth = 0.45;

bb(kBar).FaceColor = myColors(kBar , :);

end

ax = gca;
% ax.YScale = 'log';
ax.Box = 'off';
% ylim(10^-12 * [5 13])
% xlim( [0.6 2.4])
ax.FontSize = 14;
ax.YTickLabel = {'VIS' , 'SMN', 'DAN', 'SVAN', 'LIM', 'FPCN', 'DMN', 'TPJ','SUB'};
% ylim([10^-14 8*10^-11])

ax.YGrid = 'on';ax.XGrid = 'on';

legend(ax  , {'intra' , 'inter'} ,'Box' , 'off', ...
    'FontSize' , 14 , 'NumColumns', 2)


ax.XTick = [-0.5 -0.2 0.2 0.5]