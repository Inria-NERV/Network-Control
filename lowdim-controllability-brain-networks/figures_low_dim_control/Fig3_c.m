
clear all;
%% add path
subdir = genpath("fig_utils");
addpath(subdir);
subdir = genpath("results_data");
addpath(subdir);

%% load matrix and info table


ROI_info_Table = readtable( 'results_data/Schaefer200_allinfo.csv');
networksAssignment = table2array(ROI_info_Table(:,5));
nNetworks = max(networksAssignment);

n = 214;
coord = table2array(ROI_info_Table(:,7:9));
coord = coord - ones(n,1)* mean(coord);


color = table2array(ROI_info_Table(:,13 : 15))/255;


alphaBrain = 0.2;
alphaBrain = 1;

range = [0 1];

matrix = ones(n,n);
% matrix(1,1)

target =1:n;

matrixTarget = zeros(size(matrix));
matrixTarget(target,target) = matrix(target,target) ;

leftIdx = [1:100 201:207];
rightIdx = [101:200 208:214];

lab = [];
for k = 1:n
    lab = [lab ; ''];
end

lL = min(matrixTarget(matrixTarget ~= 0));
Aviz  = zeros(n,n);
Aviz(target , target)= log10(1 + matrixTarget / lL);
n = size(Aviz,1);

%%%%%%%%%%
fracViz = 0.005;
thViz = min(min(Aviz)) +  fracViz * (max(max(Aviz)) - min(min(Aviz)));
Aviz=double(Aviz < thViz).*Aviz;
Aviz = (Aviz + Aviz')/2;
str = sum(Aviz);
% subplot(1,2,2)
G = graph(Aviz);

Gright = graph(Aviz(rightIdx , rightIdx));
Gleft = graph(Aviz(leftIdx , leftIdx));

%% load results for subject
load('results_data/bigTensor_withNegative_results_fullBrainControl_on_6134_subjects.mat')

avgControllability_groupAvg = mean((averageControllabilityFullBrainTensor),2);

allGramMinBasedTensor(allGramMinBasedTensor==-Inf) = nan;

% % worstControllability_groupAvg = mean(10.^(allGramMinBasedTensor),2 ,'omitmissing' );
% worstControllability_groupAvg = mean((allGramMinBasedTensor),2 ,'omitnan' );% abs?

lowDimTensor = reshape((lapGramMinBasedTensor(:,5,:)) , [n , 6134]);

% lowDimTensor(lowDimTensor==-Inf) = nan;
% lowControllability_groupAvg = mean((lowDimTensor),2,'omitnan' );
% 

%% look for significance difference between high and low
% figure;
% histogram(allGramMinBasedTensor(1,:))

% maxHigh = max((allGramMinBasedTensor));
% minHigh = min((allGramMinBasedTensor));
%
% highTensorCentered = (exp(allGramMinBasedTensor) - ones(n,1)*minHigh) ./(ones(n,1)*(maxHigh-minHigh));

Z_high = zscore(log(abs(allGramMinBasedTensor)),0,1); % abs?
Z_low = zscore(log(abs(lowDimTensor)),0,1);


worstControllability_groupAvg = mean((Z_high),2 ,'omitnan' );% abs?
lowControllability_groupAvg = mean((Z_low),2,'omitnan' );

figure;histogram(worstControllability_groupAvg);
hold on;
histogram(lowControllability_groupAvg);

ax = gca;
% Z_high = (allGramMinBasedTensor);
% Z_low = (lowDimTensor);

% figure;
% histogram(Z_high(1,:))

tVariable = nan(1,n);
pBonferoni = 0.0000000001/n;
pVals = zeros(1,n);
cohenD = zeros(1,n);

for k =1:n
    % [h,p,ci,stats] = ttest( Z_low(k,:) , Z_high(k,:));
    [p,h,stats] = signtest(Z_low(k,:) , Z_high(k,:));

    d = computeCohen_d(Z_low(k,:) , Z_high(k,:), 'paired');

    cohenD(k)=d;
    if p<pBonferoni && abs(d)>0.5
        tVariable(k) = stats.zval;
        pVals(k) = p;
    end
end

figure;
histogram(tVariable,100)

%% loop over low dim 4 views  r=1,5,all
figure;
histogram(tVariable,100)

transparent = 1;

if transparent
    epsPlot = 100;
    alphaBrain = 0.2;
else
    epsPlot = 100;
    alphaBrain = 1;
end

range = [0 1];
% subplotIdx = [1 2 7 8;...
%     3 4 9 10;...
%     5 6 11 12];

subplotIdx = [1 2 3 4;...
    5 6 7 8;...
    9 10 11 12];

leftIdx = [1:100 201:207];
rightIdx = [101:200 208:214];
lowDimArray = [1:6];

fig = figure;
fig.Position(3:4) = 500*[ 0.7 1.5];
[ha, pos] = tight_subplot(6,2,[.005 .005],[.01 .01],[.01 .01]);

color = table2array(ROI_info_Table(:,13 : 15))/255;


for kDim = 1:3

    if kDim ==1
        controlMetricChosen = ((worstControllability_groupAvg));
        myTitle = 'worst-case';
    elseif kDim ==2
        controlMetricChosen = ((lowControllability_groupAvg));
        myTitle = sprintf( 'low dim.');
    else

        %%%%%% look significance difference
        controlMetricChosen = abs(tVariable);

        myTitle = sprintf( 'worst VS low dim. : t-stat' );
    end


    
    controlMetricChosen(isnan(controlMetricChosen))=1;

    sizeMetricLin = (controlMetricChosen - min(controlMetricChosen)) / (max(controlMetricChosen) - min(controlMetricChosen));
    sizeMetric = sizeMetricLin.^1.5;
    markerSize = 1 + 9*(sizeMetric - min(sizeMetric))/ (max(sizeMetric - min(sizeMetric)));



    if kDim <3
        %%%%% loop over views

        for kView = 1:4
            %         subplot(6,6, subplotIdx(netIdx , kView))
            axes(ha(subplotIdx(kDim , kView)))

            if kView == 1
                h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'lh');
                hold on;
                plot((Gleft) , 'MarkerSize' , markerSize(leftIdx), 'NodeColor' , color(leftIdx , :) ,  'XData' , coord(leftIdx,1)- epsPlot...
                    , 'NodeLabel' , lab ,'YData' , coord(leftIdx,2) , 'ZData' , coord(leftIdx,3) , 'EdgeAlpha' , 0.5 , 'EdgeColor' , 'k' , 'LineWidth', 2 )

                view([-90 0])
            elseif kView == 2
                h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'rh');
                hold on;
                plot((Gleft) , 'MarkerSize' , markerSize(rightIdx), 'NodeColor' , color(rightIdx , :) ,  'XData' , coord(rightIdx,1)+ epsPlot...
                    , 'NodeLabel' , lab ,'YData' , coord(rightIdx,2) , 'ZData' , coord(rightIdx,3) , 'EdgeAlpha' , 0.5 , 'EdgeColor' , 'k' , 'LineWidth', 2 )

                view([90 0])
            elseif kView == 3
                h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'lh');
                hold on;
                plot((Gleft) , 'MarkerSize' , markerSize(leftIdx), 'NodeColor' , color(leftIdx , :) ,  'XData' , coord(leftIdx,1) + epsPlot...
                    , 'NodeLabel' , lab ,'YData' , coord(leftIdx,2) , 'ZData' , coord(leftIdx,3) , 'EdgeAlpha' , 0.5 , 'EdgeColor' , 'k' , 'LineWidth', 2 )

                view([90 0])
            elseif kView == 4
                h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'rh');
                hold on;
                plot((Gleft) , 'MarkerSize' , markerSize(rightIdx), 'NodeColor' , color(rightIdx , :) ,  'XData' , coord(rightIdx,1)- epsPlot...
                    , 'NodeLabel' , lab ,'YData' , coord(rightIdx,2) , 'ZData' , coord(rightIdx,3) , 'EdgeAlpha' , 0.5 , 'EdgeColor' , 'k' , 'LineWidth', 2 )

                view([-90 0])
            end

            axis equal tight off
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        greyPart = 3.7;
        mkSizeScatter = 90;
        emptyNode = 2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif kDim ==3
        facecolor = color;

        nodesSigPos = find(and(~isnan(tVariable), tVariable>0));
        nodesSigNeg = find(and(~isnan(tVariable), tVariable<0));
        mkSize = 9;

        %%%%%%%%%%
        for kView = 1:4
            %         subplot(6,6, subplotIdx(netIdx , kView))
            axes(ha(subplotIdx(kDim , kView)))

            if kView == 1
                h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'lh');
                hold on;


                scatter3(coord(intersect(leftIdx , nodesSigPos),1)- epsPlot , coord(intersect(leftIdx , nodesSigPos),2)...
                    , coord(intersect(leftIdx , nodesSigPos),3),...
                    mkSizeScatter ,  color(intersect(leftIdx , nodesSigPos) , :),'filled','MarkerEdgeColor','flat'  )

                %%%%%%%%% do the neg nodes%%%%%%%%%%%%%
                
                scatter3(coord(intersect(leftIdx , nodesSigNeg),1)- epsPlot , coord(intersect(leftIdx , nodesSigNeg),2)...
                    , coord(intersect(leftIdx , nodesSigNeg),3),...
                    mkSizeScatter ,  color(intersect(leftIdx , nodesSigNeg) , :),'MarkerEdgeColor','flat' ,'LineWidth' , emptyNode)


                %%%%%%%%%%%%%
                view([-90 0])
            elseif kView == 2
                h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'rh');
                hold on;


                scatter3(coord(intersect(rightIdx , nodesSigPos),1)+ epsPlot , coord(intersect(rightIdx , nodesSigPos),2)...
                    , coord(intersect(rightIdx , nodesSigPos),3),...
                    mkSizeScatter ,  color(intersect(rightIdx , nodesSigPos) , :),'filled','MarkerEdgeColor','flat'  )

                %%%%%%%%% do the neg nodes%%%%%%%%%%%%%
                
                scatter3(coord(intersect(rightIdx , nodesSigNeg),1)+ epsPlot , coord(intersect(rightIdx , nodesSigNeg),2)...
                    , coord(intersect(rightIdx , nodesSigNeg),3),...
                    mkSizeScatter ,  color(intersect(rightIdx , nodesSigNeg) , :),'MarkerEdgeColor','flat' ,'LineWidth' , emptyNode)


                %%%%%%%%%%%%%

                view([90 0])
            elseif kView == 3
                h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'lh');
                hold on;


                                scatter3(coord(intersect(leftIdx , nodesSigPos),1)+ epsPlot , coord(intersect(leftIdx , nodesSigPos),2)...
                    , coord(intersect(leftIdx , nodesSigPos),3),...
                    mkSizeScatter ,  color(intersect(leftIdx , nodesSigPos) , :),'filled','MarkerEdgeColor','flat'  )

                %%%%%%%%% do the neg nodes%%%%%%%%%%%%%
                
                scatter3(coord(intersect(leftIdx , nodesSigNeg),1)+ epsPlot , coord(intersect(leftIdx , nodesSigNeg),2)...
                    , coord(intersect(leftIdx , nodesSigNeg),3),...
                    mkSizeScatter ,  color(intersect(leftIdx , nodesSigNeg) , :),'MarkerEdgeColor','flat' ,'LineWidth' , emptyNode)



                view([90 0])
            elseif kView == 4
                h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'rh');
                hold on;

                scatter3(coord(intersect(rightIdx , nodesSigPos),1)- epsPlot , coord(intersect(rightIdx , nodesSigPos),2)...
                    , coord(intersect(rightIdx , nodesSigPos),3),...
                    mkSizeScatter ,  color(intersect(rightIdx , nodesSigPos) , :),'filled','MarkerEdgeColor','flat'  )

                %%%%%%%%% do the neg nodes%%%%%%%%%%%%%
                
                scatter3(coord(intersect(rightIdx , nodesSigNeg),1)- epsPlot , coord(intersect(rightIdx , nodesSigNeg),2)...
                    , coord(intersect(rightIdx , nodesSigNeg),3),...
                    mkSizeScatter ,  color(intersect(rightIdx , nodesSigNeg) , :),'MarkerEdgeColor','flat' ,'LineWidth' , emptyNode)


                %%%%%%%%%%%%%

                view([-90 0])
            end

            axis equal tight off
        end
    end

end

% exportgraphics(fig , fullfile('UK-biobank/controllability_ukb/results_node_level/control_fullBrain_Figures' , ...
%   sprintf('4views_avg6k_control_fullBrain.png' )), 'Resolution' , 300);



%% colorbar
% figure;
% imagesc(-max(abs(tVariable)) : 1:max(abs(tVariable)));
%
% cc = colorbar;
%
% myColromap = jet();
%
% off = 25;
% myColromap(128-off : 128+off , :) = 0.5*ones(2*off+1 , 3);
%
% colormap(myColromap)
%
% cc.FontSize = 14;
% cc.Position(4) = 0.5;
%
% cc.Position(3) = 0.03;
% cc.Position(1) = 0.92;``




sigROI_id_pos = find(and(~isnan(tVariable) , tVariable>0));
sigROI_pos = ROI_info_Table(sigROI_id_pos,[1 6])


sigROI_id_neg = find(and(~isnan(tVariable) , tVariable<0));
sigROI_neg = ROI_info_Table(sigROI_id_neg,[1 6])



cohenDPos = cohenD(sigROI_id_pos)'

cohenDNeg = cohenD(sigROI_id_neg)'


