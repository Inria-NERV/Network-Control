%% plot lateralization of 9 nets quantitative approach
clear all
%% add path 
subdir = genpath("fig_utils");
addpath(subdir);
subdir = genpath("results_data");
addpath(subdir);


%%
ROI_info_Table = readtable( 'results_data/Schaefer200_allinfo.csv');
networksAssignment = table2array(ROI_info_Table(:,5));
nNetworks = max(networksAssignment);


load('results_data/lateralizationOf9netsNoLog_tensor_results_on_6134_subjects_lowdim5.mat')


%% plot bars Inside only r= 5 worst  and size rescale fig
Networks9CellVert = {'VIS'; 'SMN'; 'DAN';  'SVAN';...
                'LIM'; 'FPCN'; 'DMN'; 'TPJ' ; 'SUB'};

%%%% lateralization of size

% figure; 
cc = histogram(networksAssignment , 'Visible','off');
nodesPersys = cc.Values;
% figure;histogram(networksAssignment(101:200))
sizeRight = [12 18 11 15 8 19 13 4 7]./nodesPersys;


% figure;histogram(networksAssignment(1:100))
sizeLeft = [12 16 11 11 6 18 24 2 7]./nodesPersys;


sizeLateralization = (sizeRight - sizeLeft)./(sizeRight + sizeLeft);

% allMetricsMean = [mean(strengthLateralization);mean(avgContLateralization);...
%     mean(worstCaseContLateralization);mean(lowDimContLateralization)]';
% 
% allMetricsStd = [std(strengthLateralization);std(avgContLateralization);...
%     std(worstCaseContLateralization);std(lowDimContLateralization)]';


allMetricsMean = [sizeLateralization;
    mean(0.5*worstCaseContLateralization , 'omitnan');mean(0.5*lowDimContLateralization , 'omitnan')]';

allMetricsStd = [nan(1,nNetworks);...
    std(0.5*worstCaseContLateralization , 'omitnan');std(0.5*lowDimContLateralization , 'omitnan')]';

%%%% take MSE
allMetricsStd = allMetricsStd/sqrt(nSubjects);

% allMetricsMean(8,1) = -0.26;

%%%%%%%% sort 
[bb , ii] = sort(allMetricsMean(:,3));
Networks9CellVertOrd = {};
for k = 1:nNetworks 
    Networks9CellVertOrd{k,1} = Networks9CellVert{ii(k),1};
end
%%%%%%%%%%%%%%%%
fig = figure;
% fig.Position(3:4) = [250 300];
fig.Position(3:4)=250*[1 1];

% myBars = barh(lateralizationCoordInside ) ; %- lateralizationCoordInsideStrength);
myBars = barh(allMetricsMean(ii,:) , 'BarWidth' , 0.7 );


allMetricsStdRight = allMetricsStd;
allMetricsStdLeft = allMetricsStd;

allMetricsStdRight(2,3)=nan;allMetricsStdRight(6,3)=nan;allMetricsStdRight(7,3)=nan;
allMetricsStdLeft(1,3)=nan;allMetricsStdLeft(3,3)=nan;allMetricsStdLeft(4,3)=nan;allMetricsStdLeft(5,3)=nan;allMetricsStdLeft(8,3)=nan;allMetricsStdLeft(9,3)=nan;

allMetricsStdRight(:,2)=nan(9,1);allMetricsStdRight(8,2)=allMetricsStd(8,2);
allMetricsStdLeft(8,2)=nan;

for kMet = 1:3
    hold on;
    err = errorbar(allMetricsMean((ii),kMet) , myBars(kMet).XEndPoints , allMetricsStdLeft(ii,kMet), allMetricsStdRight(ii,kMet) ...
        ,'horizontal', ...
        'k' , 'LineStyle' , 'none');
end
% myBars(1).FaceColor = 

cOrd = colororder;
barColor = [rgb('Black');0.3*[1 1 1];rgb('LightGray') ;rgb('SteelBlue')];
barColor = [0.78*[1 1 1];cOrd(1,:);cOrd(2,:);rgb('SteelBlue')];

for kDim = 1: 3
    myBars(kDim).FaceColor = 'flat' ;
    myBars(kDim).CData = barColor(kDim , :);
    myBars(kDim).BarWidth = 0.8;
end

ax = gca;
ax.YTickLabel = {'','','','','','','','',''};
% ax.YTickLabelRotation = 45;
% xlabel('controllability lateralization (R-L)')

% ax.XLim = [-0.25 0.25];
% xline(xLims(1) ,"Label", "LEFT" , 'LabelHorizontalAlignment', 'right' , 'LabelVerticalAlignment', 'bottom')
% xline(xLims(2) ,"Label", "RIGHT" , 'LabelHorizontalAlignment', 'left' , 'LabelVerticalAlignment', 'bottom')

% ax.FontWeight = 'bold';
ax.Box = 'off';
ax.FontSize = 16;
% ax.FontSmoothing = 'on';
% axis square
ax.XTick = [-1 :0.5 : 1];
ax.XLim = [-0.6  1];

% ax.XGrid = 'on';
% ax.YGrid = 'on';

ax.YTickLabel = Networks9CellVertOrd;
% ylim([0.5 10])
leg = legend({'#nodes' , 'worst' , 'low dim.'}, 'Box' , 'off' , ...
    'FontSize' , 14 , 'NumColumns',3 , 'Location' , 'North');

leg.Position(1) = 0.01;
% xlim([-0.7 0.7])


%% same with ordering

Networks9CellVert = {'VIS'; 'SMN'; 'DAN';  'SVAN';...
                'LIM'; 'FPCN'; 'DMN'; 'TPJ' ; 'SUB'};

%%%% lateralization of size

% figure; 
cc = histogram(networksAssignment , 'Visible','off');
nodesPersys = cc.Values;
% figure;histogram(networksAssignment(101:200))
sizeRight = [12 18 11 15 8 19 13 4 7]./nodesPersys;


% figure;histogram(networksAssignment(1:100))
sizeLeft = [12 16 11 11 6 18 24 2 7]./nodesPersys;


sizeLateralization = (sizeRight - sizeLeft)./(sizeRight + sizeLeft);

% allMetricsMean = [mean(strengthLateralization);mean(avgContLateralization);...
%     mean(worstCaseContLateralization);mean(lowDimContLateralization)]';
% 
% allMetricsStd = [std(strengthLateralization);std(avgContLateralization);...
%     std(worstCaseContLateralization);std(lowDimContLateralization)]';


allMetricsMeanOrig = [sizeLateralization;
    mean(0.5*worstCaseContLateralization , 'omitnan');...
    mean(0.5*lowDimContLateralization , 'omitnan')]';


allMetricsMean = allMetricsMeanOrig;
allMetricsMean(allMetricsMean>0) = log(allMetricsMean(allMetricsMean>0));
allMetricsMean(allMetricsMean<0) = -log(-allMetricsMean(allMetricsMean<0));


%allMetricsMean = log(1+allMetricsMeanOrig);


allMetricsStd = [nan(1,nNetworks);...
    std(0.5*worstCaseContLateralization , 'omitnan');...
    std(0.5*lowDimContLateralization , 'omitnan')]';

allMetricsStd=allMetricsStd/sqrt(6134);

%%%%% sort 

[bb , ii] = sort(allMetricsMean(:,3));
Networks9CellVertOrd = {};
for k = 1:nNetworks 
    Networks9CellVertOrd{k,1} = Networks9CellVert{ii(k),1};
end


fig = figure;
% fig.Position(3:4) = [250 300];
fig.Position(3:4)=250*[1 1];

% myBars = barh(lateralizationCoordInside ) ; %- lateralizationCoordInsideStrength);
myBars = barh(allMetricsMean(ii,:) , 'BarWidth' , 0.7 );


allMetricsStdRight = allMetricsStd;
allMetricsStdLeft = allMetricsStd;

allMetricsStdRight(2,3)=nan;allMetricsStdRight(6,3)=nan;allMetricsStdRight(7,3)=nan;
allMetricsStdLeft(1,3)=nan;allMetricsStdLeft(3,3)=nan;allMetricsStdLeft(4,3)=nan;allMetricsStdLeft(5,3)=nan;allMetricsStdLeft(8,3)=nan;allMetricsStdLeft(9,3)=nan;

allMetricsStdRight(:,2)=nan(9,1);allMetricsStdRight(8,2)=allMetricsStd(8,2);
allMetricsStdLeft(8,2)=nan;

for kMet = 1:3
    hold on;
    xEnd = myBars(kMet).XEndPoints ;
    barLeft = allMetricsStdLeft(:,kMet);
    barRight = allMetricsStdRight(:,kMet);

    err = errorbar(allMetricsMean((ii),kMet) , xEnd , barLeft(ii), barRight(ii) ...
        ,'horizontal', ...
        'k' , 'LineStyle' , 'none');
end
% myBars(1).FaceColor = 

cOrd = colororder;
barColor = [rgb('Black');0.3*[1 1 1];rgb('LightGray') ;rgb('SteelBlue')];
barColor = [0.78*[1 1 1];cOrd(1,:);cOrd(2,:);rgb('SteelBlue')];

for kDim = 1: 3
    myBars(kDim).FaceColor = 'flat' ;
    myBars(kDim).CData = barColor(kDim , :);
    myBars(kDim).BarWidth = 0.8;
end

ax = gca;
ax.YTickLabel = {'','','','','','','','',''};
% ax.YTickLabelRotation = 45;
% xlabel('controllability lateralization (R-L)')

% ax.XLim = [-0.25 0.25];
% xline(xLims(1) ,"Label", "LEFT" , 'LabelHorizontalAlignment', 'right' , 'LabelVerticalAlignment', 'bottom')
% xline(xLims(2) ,"Label", "RIGHT" , 'LabelHorizontalAlignment', 'left' , 'LabelVerticalAlignment', 'bottom')

% ax.FontWeight = 'bold';
ax.Box = 'off';
ax.FontSize = 16;
% ax.FontSmoothing = 'on';
% axis square
% ax.XTick = [-1 :0.5 : 1];


ax.XGrid = 'on';
ax.YGrid = 'on';

ax.YTickLabel = Networks9CellVertOrd;
% ylim([0.5 10])
leg = legend({'#nodes' , 'worst' , 'low dim.'}, 'Box' , 'off' , ...
    'FontSize' , 14 , 'NumColumns',3 , 'Location' , 'North');

leg.Position(1) = 0.01;
% xlim([-0.7 0.7])

% ax.XLim(1) = -ax.XLim(2);


%% controllability is not only due to network asymetry

% X = (0.5*lowDimContLateralization(:,ii));
% Y = ones(nSubjects,1)*sizeLateralization(ii);
%
% Z = X-Y;
% %%%%%% ANOVA Ratio
% [p,tbl,stats] = anova1(Z);
% [c,m,h,gnames] = multcompare(stats);
%
%
% gnames = {'VIS'; 'SMN'; 'DAN';  'SVAN';...
%                 'LIM'; 'FPCN'; 'DMN'; 'TPJ' ; 'SUB'};
%
% tbl = array2table(c,"VariableNames", ...
%     ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
% tbl.("Group A") = gnames(tbl.("Group A"));
% tbl.("Group B") = gnames(tbl.("Group B"))

% pVal = zeros(1,nNetworks);
% for kNet = 1:nNetworks
%     X = (0.5*lowDimContLateralization(:,(kNet))) ;
%     [p,h,stats] = signtest(X , sizeLateralization((kNet)));
%     pVal(kNet) = p;
% end

%% regression between the network left-right  asymmetry and the mean lateralization of the slef-regulation control
X = allMetricsMean(:,1) ; % size lateralization
Y = allMetricsMean(:,2) ; % self regulation lateralization

[b,bint,r,rint,stats] = regress(Y,[ones(nNetworks,1) X]);

% stats are R2 , F, p , error variance
stats
