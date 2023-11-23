clear all

load('results_data/single_driver_onlyEigen_Xf10/workSpace_singleDriver_onlyEigen_rho4_ord0_Xf10_n_256_run_6.mat')

%% add path 
subdir = genpath("fig_utils");
addpath(subdir);
subdir = genpath("results_data");
addpath(subdir);

%% 
% plot one graph several eigen maps

minNodeSize = 4 ;
markerSizeFact = 4;

lab = [];
for k = 1:n
    lab = [lab ; ''];
end



fig = figure;
% fig.Position(3:4) = [814 687];
fig.Position(3:4) = [1057 542];
fig.Color = [1 1 1];
[ha, pos] = tight_subplot(2,3,[.035 .02],[.1 .01],[.01 .01]);


%%%%%%%%%%%%%
matrix = graphs{1,1};
% matrix = matrix(target , target);
if issparse(matrix)
    matrix = full(matrix);
end


% cut the graph
% matrix(1:n/2 , n/2+1 : n) = 0;
% matrix(n/2+1 : n , 1:n/2) = 0;


isSym = norm(matrix - matrix') < 10^(-6);
if isSym
    G = graph(matrix);
else
    G = digraph(matrix');
end


mode = 1:6;
% mode = 251:256;

% %%%% normalize Laplacian?
Ggsp = gsp_graph(matrix);
Ggsp = gsp_create_laplacian(Ggsp, 'combinatorial');
% Ggsp = gsp_create_laplacian(Ggsp, 'normalized');
targetlaplac = full(Ggsp.L);

[VnotSorted,lambdaMatNotSorted] = eig( targetlaplac );
lambdaNotSorted = diag(lambdaMatNotSorted);

[lambda , idxAsc] = sort(lambdaNotSorted , 'ascend'); % for Laplacian
V = VnotSorted(: , idxAsc);



%%%%% plot
for kMode = (1:6)

    ax = ha(kMode);
    axes(ax)

    markerSizeMet =  sum(matrix);

    markerSize = minNodeSize + markerSizeFact*(markerSizeMet - min(markerSizeMet))/ (max(markerSizeMet) - min(markerSizeMet));

    myColorMap = jet();


    colorMetricBrut = V(:,mode(kMode));
    % findIndices
    if min(colorMetricBrut) ~= max(colorMetricBrut)
        colorMetric = (colorMetricBrut - min(colorMetricBrut)) / ( max(colorMetricBrut) - min(colorMetricBrut));
    else
        colorMetric = myMet;
    end
    myMin = min(colorMetricBrut);
    myMax = max(colorMetricBrut);
    absMax = max(abs([myMin myMax]));
    colorMetric = colorMetricBrut/ absMax ;
    %     colorMapIdx = floor(254*(colorMetric)) + 1 ; % if jet
    colorMapIdx = floor(127*(colorMetric)) + 128 ; % if jet Abs
    % colorMapIdx = floor(31*(colorMetric)) + 32 ; % if polarmap


    color = myColorMap(colorMapIdx , :);

    % plot black first
    plot( G , 'MarkerSize' , markerSize+4.2 , 'NodeColor' , 0.3*[1 1 1] , 'NodeLabel' , lab ,...
        'EdgeAlpha' , 0.05 , 'EdgeColor' , 'k' )% , 'Layout' , 'subspace' , 'Dimension' , dimViz(kScale))

    hold on;

    % then plot color
    plot( G , 'MarkerSize' , markerSize, 'NodeColor' , color , 'NodeLabel' , lab ,...
        'EdgeAlpha' , 0.05 , 'EdgeColor' , 'k' )% , 'Layout' , 'subspace' , 'Dimension' , dimViz(kScale))
    hold on;
    %     axis equal tight off
    % set(ax, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
    % xlabel(ax, sprintf('#modules %i \n' , kMode) , 'FontWeight', 'bold' , 'FontSize' , 14)
    %     polarmap();
    %     ax = gca;
    ax.Box = 'off';

    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';

    ax.Title.String = sprintf('eigenmap V_%i \n' ,  kMode);

    ax.Title.VerticalAlignment = 'top';
    %     ax.Title.HorizontalAlignment = 'right';
    ax.Title.Margin = 4;
    ax.Title.FontSize = 13.2;
end


