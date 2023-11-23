%% plot systems controllability 9nets
clear all

%% load the tensor of system level controllability for all subjects
load('results_data/controllabilityTensor_9nets.mat')
% load('results_data/controllabilityTensor_9nets_noLog.mat')

rowSum = sum(controlTensor , 2);
okRows = rowSum ~= -Inf;
% avgcontrolSystems = mean(10.^controlTensor(okRows,:));
% stdcontrolSystems = std(10.^controlTensor(okRows,:));

avgcontrolSystems = mean(controlTensor(okRows,:));
stdcontrolSystems = std(controlTensor(okRows,:));

nNetworks = 9;
% load colormap schaefer
load('mySchaeferColorMap_9nets.mat')

Networks9CellVert = {'VIS'; 'SMN'; 'DAN';  'SVAN';...
                'LIM'; 'FPCN'; 'DMN'; 'TPJ' ; 'SUB'};

%% look driverness and targetness per subject
nSubjects= 6134;

drivernessArray = nan(nSubjects,nNetworks);
targetnessArray = nan(nSubjects,nNetworks);

for kSub = 1 : nSubjects
    kSub
    if okRows(kSub)
        metagraph = reshape(10.^controlTensor(kSub,:) , [nNetworks , nNetworks]);
        drivernessArray(kSub,:) = sum(metagraph) ;
        targetnessArray(kSub,:) = sum(metagraph , 2)';
    end
end


%% plotting asym driver target
avgDriverness = mean(drivernessArray , 'omitnan');
stdDriverness = std(drivernessArray , 'omitnan');

avgTargetness = mean(targetnessArray , 'omitnan');
stdTargetness = std(targetnessArray , 'omitnan');


driverTargetAsymArray = drivernessArray - targetnessArray;

% driverTargetAsymArray = log(abs(drivernessArray)./abs(targetnessArray));


asymMetric = mean(driverTargetAsymArray, 'omitnan');
asymMetricErr = std(driverTargetAsymArray, 'omitnan');

[mm , ii] = sort(asymMetric);

fig = figure;
fig.Position(3:4) =[ 280  360];

% scatter(driverness , targetness , 160*ones(1,nNetworks) , 5*ones(1,nNetworks) , 'filled' , 'Marker', 'o')
% hold on;
% scatter(driverness , targetness , 130*ones(1,nNetworks) , 1:nNetworks , 'filled' , 'Marker', 'o')
% hold on;

bb  = barh(asymMetric(ii));
bb.FaceColor = "flat";

myLabels = {};
for kNet = 1:nNetworks
    myLabels{1 , kNet} = Networks9CellVert{ii(kNet),1};

    bb.CData(kNet , :) = myColorMap(ii(kNet) , :);
end

ax = gca;
ax.Box = 'off';
ax.YTickLabel = myLabels;
% ax.YTickLabel = {'','','','','','','','','',};

% ax.XTick = -1:0.4:1;

ax.XGrid = 'on';
ax.YGrid = 'on';
% xlim([-0.6 0.6])
% axis equal square
% xlabel('log_1_0( driverness )')
% ylabel('log_1_0( targetness )')
% ax.FontWeight = 'bold';
ax.FontSize = 16;

% only show half bars
asymErrPos = nan(size(asymMetricErr(ii)));
asymErrNeg = nan(size(asymMetricErr(ii)));

asymErrNeg(1:5) = asymMetricErr(ii(1:5));
asymErrPos(6:9) = asymMetricErr(ii(6:9));


asymMetricErrNeg = asymMetricErr;
asymMetricErrNeg(asymMetric>0) = nan;

asymMetricErrPos = asymMetricErr;
asymMetricErrPos(asymMetric<0) = nan;


hold on;
% err = errorbar(asymMetric(ii) , bb.XEndPoints , asymErrNeg, asymErrPos ,'horizontal', ...
%     'k' , 'LineStyle' , 'none');
err = errorbar(asymMetric(ii) , bb.XEndPoints , asymMetricErrNeg(ii) , asymMetricErrPos(ii) ,'horizontal', ...
    'k' , 'LineStyle' , 'none');

err.LineWidth = 0.8;
% % ax.XScale = 'log';
% xlim(10^-8 * [-0.22 0.22])
% ax.XTick = 10^-8 * [-0.15 : 0.05: 0.15];

% ax.XTick = 10^-8 * [-0.2 : 0.1: 0.2];



%% plotting asym driver target horizental

targetnessArray(targetnessArray<0) = nan;
targetnessArray(targetnessArray==0) = nan;

avgDriverness = mean(drivernessArray , 'omitnan');
stdDriverness = std(drivernessArray , 'omitnan');

avgTargetness = mean(targetnessArray , 'omitnan');
stdTargetness = std(targetnessArray , 'omitnan');


driverTargetAsymArray = ((drivernessArray)) - ((targetnessArray));


% driverTargetAsymArray = log(abs(drivernessArray)) - log(abs(targetnessArray));

% driverTargetAsymArray = log(abs(drivernessArray)./abs(targetnessArray));
% 
% driverTargetAsymArray = log(abs(drivernessArray)) - log(abs(targetnessArray));


asymMetric = mean(driverTargetAsymArray, 'omitnan');
asymMetricErr = std(driverTargetAsymArray, 'omitnan');

[mm , ii] = sort(asymMetric , 'descend');

fig = figure;
fig.Position(3:4) =[ 280  200];

% scatter(driverness , targetness , 160*ones(1,nNetworks) , 5*ones(1,nNetworks) , 'filled' , 'Marker', 'o')
% hold on;
% scatter(driverness , targetness , 130*ones(1,nNetworks) , 1:nNetworks , 'filled' , 'Marker', 'o')
% hold on;

bb  = bar(asymMetric(ii));
bb.FaceColor = "flat";

myLabels = {};
for kNet = 1:nNetworks
    myLabels{1 , kNet} = Networks9CellVert{ii(kNet),1};

    bb.CData(kNet , :) = myColorMap(ii(kNet) , :);
end

ax = gca;
ax.Box = 'off';
ax.XTickLabel = myLabels;
% ax.YTickLabel = {'','','','','','','','','',};

% ax.XTick = -1:0.4:1;

ax.XGrid = 'off';
ax.YGrid = 'off';
% xlim([-0.6 0.6])
% axis equal square
% xlabel('log_1_0( driverness )')
% ylabel('log_1_0( targetness )')
% ax.FontWeight = 'bold';
ax.FontSize = 16;

% only show half bars
asymErrPos = nan(size(asymMetricErr(ii)));
asymErrNeg = nan(size(asymMetricErr(ii)));

asymErrNeg(1:5) = asymMetricErr(ii(1:5));
asymErrPos(6:9) = asymMetricErr(ii(6:9));


asymMetricErrNeg = asymMetricErr;
asymMetricErrNeg(asymMetric>0) = nan;

asymMetricErrPos = asymMetricErr;
asymMetricErrPos(asymMetric<0) = nan;


hold on;
% err = errorbar(asymMetric(ii) , bb.XEndPoints , asymErrNeg, asymErrPos ,'horizontal', ...
%     'k' , 'LineStyle' , 'none');
err = errorbar( 1:nNetworks, bb.YEndPoints , asymMetricErrNeg(ii) , asymMetricErrPos(ii) , ...
    'k' , 'LineStyle' , 'none');

err.LineWidth = 0.8;
% ax.YScale = 'log';

% ax.YLim=[-8 9];
ax.XTickLabelRotation=45;
% xlim(10^-8 * [-0.22 0.22])
% ax.XTick = 10^-8 * [-0.15 : 0.05: 0.15];

% ax.XTick = 10^-8 * [-0.2 : 0.1: 0.2];




%% plotting asym driver target horizental split in two

targetnessArray(targetnessArray<0) = nan;
targetnessArray(targetnessArray==0) = nan;

avgDriverness = mean(drivernessArray , 'omitnan');
stdDriverness = std(drivernessArray , 'omitnan');

avgTargetness = mean(targetnessArray , 'omitnan');
stdTargetness = std(targetnessArray , 'omitnan');


driverTargetAsymArray = ((drivernessArray)) - ((targetnessArray));


%% the stat
%%%%%% ANOVA all
[p,tbl,stats] = anova1(driverTargetAsymArray);
[c,m,h,gnames] = multcompare(stats);


gnames = {'VIS'; 'SMN'; 'DAN';  'SVAN';...
                'LIM'; 'FPCN'; 'DMN'; 'TPJ' ; 'SUB'};

tbl = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl.("Group A") = gnames(tbl.("Group A"));
tbl.("Group B") = gnames(tbl.("Group B"))

%% the plot 
% driverTargetAsymArray = log(abs(drivernessArray)) - log(abs(targetnessArray));

% driverTargetAsymArray = log(abs(drivernessArray)./abs(targetnessArray));
% 
% driverTargetAsymArray = log(abs(drivernessArray)) - log(abs(targetnessArray));


asymMetric = mean(driverTargetAsymArray, 'omitnan');
asymMetricErr = std(driverTargetAsymArray, 'omitnan');


asymMetricPos = asymMetric;
asymMetricPos(asymMetric < 0 ) = nan;
asymMetricErrPos = asymMetricErr;
asymMetricErrPos(asymMetricErr < 0 ) = nan;


asymMetricNeg = asymMetric;
asymMetricNeg(asymMetric > 0 ) = nan;
asymMetricErrNeg = asymMetricErr;
asymMetricErrNeg(asymMetricErr > 0 ) = nan;



[mm , ii] = sort(asymMetric , 'descend');

fig = figure;
fig.Position(3:4) =[ 450  300];


% scatter(driverness , targetness , 160*ones(1,nNetworks) , 5*ones(1,nNetworks) , 'filled' , 'Marker', 'o')
% hold on;
% scatter(driverness , targetness , 130*ones(1,nNetworks) , 1:nNetworks , 'filled' , 'Marker', 'o')
% hold on;

%%%%% pos 
subplot(2,1,1)
bb  = bar(asymMetricPos(ii));
bb.FaceColor = "flat";

myLabels = {};
for kNet = 1:nNetworks
    myLabels{1 , kNet} = Networks9CellVert{ii(kNet),1};

    bb.CData(kNet , :) = myColorMap(ii(kNet) , :);
end

ax = gca;
ax.Box = 'off';
ax.XTickLabel = {};
% ax.YTickLabel = {'','','','','','','','','',};

% ax.XTick = -1:0.4:1;

ax.XGrid = 'off';
ax.YGrid = 'off';
% xlim([-0.6 0.6])
% axis equal square
% xlabel('log_1_0( driverness )')
% ylabel('log_1_0( targetness )')
% ax.FontWeight = 'bold';
ax.FontSize = 16;

% only show half bars
asymErrPos = nan(size(asymMetricErr(ii)));
asymErrNeg = nan(size(asymMetricErr(ii)));

asymErrNeg(1:5) = asymMetricErr(ii(1:5));
asymErrPos(6:9) = asymMetricErr(ii(6:9));


asymMetricErrNeg = asymMetricErr;
asymMetricErrNeg(asymMetric>0) = nan;

asymMetricErrPos = asymMetricErr;
asymMetricErrPos(asymMetric<0) = nan;


hold on;
% err = errorbar(asymMetric(ii) , bb.XEndPoints , asymErrNeg, asymErrPos ,'horizontal', ...
%     'k' , 'LineStyle' , 'none');
err = errorbar( 1:nNetworks, bb.YEndPoints , asymMetricErrNeg(ii) , asymMetricErrPos(ii) , ...
    'k' , 'LineStyle' , 'none');

err.LineWidth = 0.8;
% ax.YScale = 'log';

% ax.YLim=[-8 9];
ax.XTickLabelRotation=45;
% xlim(10^-8 * [-0.22 0.22])
% ax.XTick = 10^-8 * [-0.15 : 0.05: 0.15];

% ax.XTick = 10^-8 * [-0.2 : 0.1: 0.2];
ax.YScale = "log";

bb.BaseValue = 10^-11;
ax1 = ax;
ax.YTick = [10^-10  10^-9];
ax.YLim = [10^-10.5  10^-8.5];

xlim([0.25 9.75])
%%%%%% do neg 
subplot(2,1,2)
bb  = bar(asymMetricNeg(ii));
bb.FaceColor = "flat";

myLabels = {};
for kNet = 1:nNetworks
    myLabels{1 , kNet} = Networks9CellVert{ii(kNet),1};

    bb.CData(kNet , :) = myColorMap(ii(kNet) , :);
end

ax = gca;
ax.Box = 'off';
ax.XTickLabel = {};
% ax.YTickLabel = {'','','','','','','','','',};

% ax.XTick = -1:0.4:1;

ax.XGrid = 'off';
ax.YGrid = 'off';
% xlim([-0.6 0.6])
% axis equal square
% xlabel('log_1_0( driverness )')
% ylabel('log_1_0( targetness )')
% ax.FontWeight = 'bold';
ax.FontSize = 16;

% only show half bars
asymErrPos = nan(size(asymMetricErr(ii)));
asymErrNeg = nan(size(asymMetricErr(ii)));

asymErrNeg(1:5) = asymMetricErr(ii(1:5));
asymErrPos(6:9) = asymMetricErr(ii(6:9));


asymMetricErrNeg = asymMetricErr;
asymMetricErrNeg(asymMetric>0) = nan;

asymMetricErrPos = asymMetricErr;
asymMetricErrPos(asymMetric<0) = nan;


hold on;
% err = errorbar(asymMetric(ii) , bb.XEndPoints , asymErrNeg, asymErrPos ,'horizontal', ...
%     'k' , 'LineStyle' , 'none');
err = errorbar( 1:nNetworks, bb.YEndPoints , asymMetricErrNeg(ii) , asymMetricErrPos(ii) , ...
    'k' , 'LineStyle' , 'none');

err.LineWidth = 0.8;
% ax.YScale = 'log';

% ax.YLim=[-8 9];
ax.XTickLabelRotation=45;
% xlim(10^-8 * [-0.22 0.22])
% ax.XTick = 10^-8 * [-0.15 : 0.05: 0.15];

% ax.XTick = 10^-8 * [-0.2 : 0.1: 0.2];
ax.YScale = "log";

ax.Position(2) = 0.23;

ax.YTick = [-10^-9  -10^-11];
ax.YLim = -[10^-8  10^-12];


ax2 = ax;

ax1.Position(2) = ax2.Position(2) + ax2.Position(4);

ax1.Position(1) = ax1.Position(1) + 0.05;
ax2.Position(1) = ax2.Position(1) + 0.05;
ax.XTickLabel = myLabels;
xlim([0.25 9.75])
