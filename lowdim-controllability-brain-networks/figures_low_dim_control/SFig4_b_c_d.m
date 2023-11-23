%% plot systems controllability 18nets
clear all

%% load the tensor of system level controllability for all subjects
Networks18Cell = {'Vis-Cent', 'Vis-Peri', 'SMN-A','SMN-B',  'DAN-A',  'DAN-B'...
    ,  'SVAN-A','SVAN-B', 'LIM-B', 'LIM-A', ...
    'FPCN-A', 'FPCN-B', 'FPCN-C', 'DMN-A','DMN-B','DMN-C', 'TPJ' , 'Sub'};

load('results_data/controllabilityTensor18nets.mat')

rowSum = sum(controlTensor , 2);
okRows = rowSum ~= -Inf;
avgcontrolSystems = mean(10.^controlTensor(okRows,:));
stdcontrolSystems = std(10.^controlTensor(okRows,:));

nNetworks = 18;
% load colormap schaefer
load('mySchaeferColorMap_18nets.mat')

%%%% this is the metagraph


%% %% metagraph only cortical

avgcontrolSystemsMat = reshape(avgcontrolSystems , [nNetworks , nNetworks]);

Networks9CellVert  = {'Vis-Cent';'Vis-Peri'; 'SMN-A';'SMN-B';  'DAN-A';  'DAN-B'...
    ;  'SVAN-A';'SVAN-B'; 'LIM-B'; 'LIM-A'; ...
    'FPCN-A'; 'FPCN-B'; 'FPCN-C';'DMN-A';'DMN-B';'DMN-C'; 'TPJ' ; 'Sub'};


Networks9CellVertCort = {'Vis-Cent', 'Vis-Peri', 'SMN-A','SMN-B',  'DAN-A',  'DAN-B'...
    ,  'SVAN-A','SVAN-B', 'LIM-B', 'LIM-A', ...
    'FPCN-A', 'FPCN-B', 'FPCN-C', 'DMN-A','DMN-B','DMN-C', 'TPJ' };

%%%% remove diagonal
avgcontrolSystemsMat2plot = avgcontrolSystemsMat;
for k = 1 : nNetworks
    avgcontrolSystemsMat2plot(k,k) = nan;
    avgcontrolSystemsMat2plot(end,k) = nan;
    avgcontrolSystemsMat2plot(k,end) = nan;
end

fig = figure;
fig.Position(3:4) = 390 * [1 1];

img = imagesc(log10(avgcontrolSystemsMat2plot(1:end-1 , 1:end-1)).^1);
% img = imagesc(log10(avgcontrolSystemsMat2plot(1:end-1 , 1:end-1)));

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


% colormap("sky")

cb = colorbar;
% cb.Label.String  = {'log_1_0(\lambda_m_i_n(W)) '};
cb.Label.FontSize = 16;

% cb.Position(1) = 1;
ax.FontSize = 13.5;

axis square

% img.AlphaData = 0.1;

img.CDataMapping = 'scaled';

ax.XGrid = 'on';
ax.YGrid = 'on';
ax.XTick = (1:17)-0.5;
ax.YTick = (1:17)-0.5;
% ax.GridLineWidth = 1;
ax.GridColor = 0.05*[1 1 1];

ax.XTickLabelRotation = 60;
ax.GridAlpha = 0.6;






%% look driverness and targetness per subject
nSubjects= 6134;

drivernessArray = nan(nSubjects,nNetworks);
targetnessArray = nan(nSubjects,nNetworks);

for kSub = 1 : nSubjects
    kSub
    if okRows(kSub)
        metagraph = reshape(controlTensor(kSub,:) , [nNetworks , nNetworks]);
        drivernessArray(kSub,:) = sum(10.^metagraph) ;
        targetnessArray(kSub,:) = sum(10.^metagraph , 2)';
    end
end




%% plotting asym driver target
avgDriverness = mean(drivernessArray , 'omitmissing');
stdDriverness = std(drivernessArray , 'omitmissing');

avgTargetness = mean(targetnessArray , 'omitmissing');
stdTargetness = std(targetnessArray , 'omitmissing');


asymMetric = mean(drivernessArray - targetnessArray, 'omitmissing');
asymMetricErr = std(drivernessArray - targetnessArray, 'omitmissing');

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
ax.YTick = 1:nNetworks;
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

hold on;
err = errorbar(asymMetric(ii) , bb.XEndPoints , asymMetricErr(ii) ,'horizontal', ...
    'k' , 'LineStyle' , 'none');


err.LineWidth = 0.8;
% ax.XScale = 'log';


xlim(10^-8 * [-0.4 0.4])
ax.XTick = 10^-8 * [-0.4 :0.2 :0.4];



%% plotting asym driver target horizental split in two

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
fig.Position(3:4) =[ 440  300];


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


ax1 = ax;
ax.YTick = [10^-10  10^-8];
ax.YLim = [10^-11.5  10^-8];
bb.BaseValue = 10^-11.5;

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
% xlim(10^-8 * [-0.22 0.22])
% ax.XTick = 10^-8 * [-0.15 : 0.05: 0.15];

% ax.XTick = 10^-8 * [-0.2 : 0.1: 0.2];
ax.YScale = "log";

ax.Position(2) = 0.23;

ax.YTick = [-10^-8  -10^-10];
ax.YLim = -[10^-8  10^-11.5];


ax2 = ax;

ax1.Position(2) = ax2.Position(2) + ax2.Position(4);

ax.XTick = 1:nNetworks;

ax.XTickLabel = myLabels;
ax.XTickLabelRotation=65;


%% plotting sum driver target
avgDriverness = mean(drivernessArray , 'omitmissing');
stdDriverness = std(drivernessArray , 'omitmissing');

avgTargetness = mean(targetnessArray , 'omitmissing');
stdTargetness = std(targetnessArray , 'omitmissing');


totalDegree = mean(drivernessArray + targetnessArray, 'omitmissing');
totalDegreeErr = std(drivernessArray + targetnessArray, 'omitmissing');

[mm , ii] = sort(totalDegree);

fig = figure;
fig.Position(3:4) =[ 280  360];

% scatter(driverness , targetness , 160*ones(1,nNetworks) , 5*ones(1,nNetworks) , 'filled' , 'Marker', 'o')
% hold on;
% scatter(driverness , targetness , 130*ones(1,nNetworks) , 1:nNetworks , 'filled' , 'Marker', 'o')
% hold on;

bb  = barh(totalDegree(ii) );
bb.FaceColor = "flat";

% bb.BaseValue = min(totalDegree) - 1;

myLabels = {};
for kNet = 1:nNetworks
    myLabels{1 , kNet} = Networks9CellVert{ii(kNet),1};

    bb.CData(kNet , :) = myColorMap(ii(kNet) , :);
end

ax = gca;
ax.Box = 'off';
ax.YTick = 1:nNetworks;
ax.YTickLabel = myLabels;
% ax.YTickLabel = {'','','','','','','','','',};

ax.XTick = 10.^[-9 -7];

% ax.XGrid = 'on';
% ax.YGrid = 'on';
% xlim([-0.6 0.6])
% axis equal square
% xlabel('log_1_0( driverness )')
% ylabel('log_1_0( targetness )')
% ax.FontWeight = 'bold';
ax.FontSize = 16;

hold on;
err = errorbar(totalDegree(ii) , bb.XEndPoints ,nan(size(totalDegreeErr(ii))), totalDegree(ii) ,'horizontal', ...
    'k' , 'LineStyle' , 'none');


err.LineWidth = 1.2;

ax.XScale = 'log';
bb.BaseValue = bb.BaseValue * 0.6;


