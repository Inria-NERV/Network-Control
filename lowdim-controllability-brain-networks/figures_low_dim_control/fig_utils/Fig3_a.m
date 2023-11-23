clear all;
%% add path 
subdir = genpath("fig_utils");
addpath(subdir);
subdir = genpath("results_data");
addpath(subdir);

%% load matrix and info table
% matrixStruct = load('UK-biobank/controllability_ukb/results_node_level/workSpaces/1000367_matrix.mat');
% matrix = matrixStruct.matrix;

n = 214;
matrix = ones(n,n);

ROI_info_Table = readtable( 'results_data/Schaefer200_allinfo.csv');
networksAssignment = table2array(ROI_info_Table(:,5));
nNetworks = max(networksAssignment);

color = table2array(ROI_info_Table(:,13:15));
color = color/max(max(color));

coord = table2array(ROI_info_Table(:,7:9));
coord = coord - ones(n,1)* mean(coord);

Networks9Cell = {'Vis', 'SomMot',  'DorsAttn',  'SalVentAttn',...
    'Limbic', 'Cont','Default', 'TempPar' , 'Sub'};
legCell = {'VIS' , 'SMN', 'DAN', 'SVAN', 'LIM', 'FPCN', 'DMN', 'TPJ', 'SUB'};
legCell8 = {'VIS' , 'SMN', 'DAN', 'SVAN', 'LIM', 'FPCN', 'DMN', 'TPJ'};
leftIdx = [1:100 201:207];
rightIdx = [101:200 208:214];

load('mySchaeferColorMap_9nets.mat')


figure;
bb = bar(1:9 , zeros(9,9));
for kNet = 1 : 9
    bb(kNet).FaceColor = 'flat';
    bb(kNet).CData = myColorMap(kNet , :);
end
% colormap(myColorMap)
leg = legend(legCell , 'Box' , 'off' , 'NumColumns',3);
leg.FontSize = 16;

figure;
bb = bar(1:8 , zeros(8,8));
for kNet = 1 : 8
    bb(kNet).FaceColor = 'flat';
    bb(kNet).CData = myColorMap(kNet , :);
end
% colormap(myColorMap)
leg = legend(legCell , 'Box' , 'off' , 'NumColumns',4);
leg.FontSize = 16;


alphaBrain = 0.36;
range = [0 1];


leftIdx = [1:100 201:207];
rightIdx = [101:200 208:214];


strength = sum(matrix);

lab = [];
for k = 1:n
    lab = [lab ; ''];
end

target = 1:n;
lL = min(matrix(matrix ~= 0));
Aviz  = zeros(n,n);
Aviz(target , target)= log10(1 + matrix / lL);
n = size(Aviz,1);

%%%%%%%%%%
fracViz = 0.05;
thViz = min(min(Aviz)) +  fracViz * (max(max(Aviz)) - min(min(Aviz)));
Aviz=double(Aviz < thViz).*Aviz;
Aviz = (Aviz + Aviz')/2;
str = sum(Aviz);
% subplot(1,2,2)

Gright = graph(Aviz(rightIdx , rightIdx));
Gleft = graph(Aviz(leftIdx , leftIdx));


%% plot scahefer parcellation all in one plot 5 views
alphaColor2 = 1;

% color = myColorMap(networksAssignment , :);

fig = figure;
fig.Position = [2200 607 580 420];

[ha, pos] = tight_subplot(1,5,[.01 .01],[.2 .01],[.01 .01]);

rangeNoir = [0 0.01];

for kView = 1:5
    %         subplot(6,6, subplotIdx(netIdx , kView))
    axes(ha(kView))

    if kView == 1


        h_brain = my_simpleBrainSurfaceHemisphere(rangeNoir , 1 , alphaBrain , 'lh');

        for netIdx = 1:9


            nodeColors = nan(n,3);
            target = find(networksAssignment == netIdx);
            nodeColors(target,:) = color(target , :);

            % networkCmap = cMapArray{1,netIdx};

            % target = find(networksAssignment == netIdx);


            %     myMetric = double(networksAssignment == netIdx);

            myMetric = zeros(n,1);
            myMetric(target) = strength(target);

            myMetric = myMetric/max(myMetric);

            hold on;
            % h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_OnlyTarget2(range  , 1 , 1 , 'lh'  , myMetric ,...
            %     cMapArray{1 , netIdx} , target );

            h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaColor2 , 'lh' , nodeColors  ) ;

        end
        view([-90 0])

    elseif kView == 2
        h_brain = my_simpleBrainSurfaceHemisphere(rangeNoir , 1 , alphaBrain , 'rh');
        for netIdx = 1:9
            nodeColors = nan(n,3);
            target = find(networksAssignment == netIdx);
            nodeColors(target,:) = color(target , :);

            % networkCmap = cMapArray{1,netIdx};
            %
            % target = find(networksAssignment == netIdx);

            %     myMetric = double(networksAssignment == netIdx);

            myMetric = zeros(n,1);
            myMetric(target) = strength(target);

            myMetric = myMetric/max(myMetric);

            hold on;
            % h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_OnlyTarget2(range  , 1 , 1 , 'rh'  , myMetric ,...
            %     cMapArray{1 , netIdx} , target );
            h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaColor2 , 'rh' , nodeColors  ) ;


        end

        view([90 0])
    elseif kView == 3
        h_brain = my_simpleBrainSurface(rangeNoir , 1 , alphaBrain );


        for netIdx = 1:9
            nodeColors = nan(n,3);
            target = find(networksAssignment == netIdx);
            nodeColors(target,:) = color(target , :);
            myMetric = zeros(n,1);
            myMetric(target) = strength(target);

            myMetric = myMetric/max(myMetric);

            hold on;
            h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaColor2 , 'lh' , nodeColors  ) ;

            % hl = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_OnlyTarget2(range  , 1 , 1 , 'lh'  , myMetric ,...
            %     cMapArray{1 , netIdx} , target );
            hold on;
            % hr = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_OnlyTarget2(range  , 1 , 1 , 'rh'  , myMetric ,...
            %     cMapArray{1 , netIdx} , target );
            h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaColor2 , 'rh' , nodeColors  ) ;

        end

        view([0 90])

    elseif kView == 4
        h_brain = my_simpleBrainSurfaceHemisphere(rangeNoir , 1 , alphaBrain , 'lh');
        for netIdx = 1:9

            nodeColors = nan(n,3);
            target = find(networksAssignment == netIdx);
            nodeColors(target,:) = color(target , :);

            myMetric = zeros(n,1);
            myMetric(target) = strength(target);

            myMetric = myMetric/max(myMetric);

            hold on;
            % h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_OnlyTarget2(range  , 1 , 1 , 'lh'  , myMetric ,...
            %     cMapArray{1 , netIdx} , target );
            h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaColor2 , 'lh' , nodeColors  ) ;

        end
        view([90 0])
    elseif kView == 5
        h_brain = my_simpleBrainSurfaceHemisphere(rangeNoir , 1 , alphaBrain , 'rh');
        for netIdx = 1:9

            nodeColors = nan(n,3);
            target = find(networksAssignment == netIdx);
            nodeColors(target,:) = color(target , :);

            myMetric = zeros(n,1);
            myMetric(target) = strength(target);

            myMetric = myMetric/max(myMetric);

            hold on;
            % h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_OnlyTarget2(range  , 1 , 1 , 'rh'  , myMetric ,...
            %     cMapArray{1 , netIdx} , target );
            h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaColor2 , 'rh' , nodeColors  ) ;
        end
        view([-90 0])
    end

    axis equal tight off
end


%     exportgraphics(fig , sprintf('UK-biobank/parcellations_code/fig_RSN_surface/RSN_%s.png' , Networks9Cell{netIdx}) , 'Resolution' , 300)
kk = 5;







%% plot scahefer parcellation all in one plot 4 views
alphaColor2 = 1;

% color = myColorMap(networksAssignment , :);

fig = figure;
fig.Position = [2200 607 580 420];

[ha, pos] = tight_subplot(2,2,[.01 .01],[.2 .01],[.01 .01]);

rangeNoir = [0 0.01];

for kView = 1:4
    %         subplot(6,6, subplotIdx(netIdx , kView))
    axes(ha(kView))

    if kView == 1


        h_brain = my_simpleBrainSurfaceHemisphere(rangeNoir , 1 , alphaBrain , 'lh');

        for netIdx = 1:9


            nodeColors = nan(n,3);
            target = find(networksAssignment == netIdx);
            nodeColors(target,:) = color(target , :);

            % networkCmap = cMapArray{1,netIdx};

            % target = find(networksAssignment == netIdx);


            %     myMetric = double(networksAssignment == netIdx);

            myMetric = zeros(n,1);
            myMetric(target) = strength(target);

            myMetric = myMetric/max(myMetric);

            hold on;
            % h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_OnlyTarget2(range  , 1 , 1 , 'lh'  , myMetric ,...
            %     cMapArray{1 , netIdx} , target );

            h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaColor2 , 'lh' , nodeColors  ) ;

        end
        view([-90 0])

    elseif kView == 2
        h_brain = my_simpleBrainSurfaceHemisphere(rangeNoir , 1 , alphaBrain , 'rh');
        for netIdx = 1:9
            nodeColors = nan(n,3);
            target = find(networksAssignment == netIdx);
            nodeColors(target,:) = color(target , :);

            % networkCmap = cMapArray{1,netIdx};
            %
            % target = find(networksAssignment == netIdx);

            %     myMetric = double(networksAssignment == netIdx);

            myMetric = zeros(n,1);
            myMetric(target) = strength(target);

            myMetric = myMetric/max(myMetric);

            hold on;
            % h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_OnlyTarget2(range  , 1 , 1 , 'rh'  , myMetric ,...
            %     cMapArray{1 , netIdx} , target );
            h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaColor2 , 'rh' , nodeColors  ) ;


        end

        view([90 0])

    elseif kView == 3
        h_brain = my_simpleBrainSurfaceHemisphere(rangeNoir , 1 , alphaBrain , 'lh');
        for netIdx = 1:9

            nodeColors = nan(n,3);
            target = find(networksAssignment == netIdx);
            nodeColors(target,:) = color(target , :);

            myMetric = zeros(n,1);
            myMetric(target) = strength(target);

            myMetric = myMetric/max(myMetric);

            hold on;
            % h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_OnlyTarget2(range  , 1 , 1 , 'lh'  , myMetric ,...
            %     cMapArray{1 , netIdx} , target );
            h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaColor2 , 'lh' , nodeColors  ) ;

        end
        view([90 0])
    elseif kView == 4
        h_brain = my_simpleBrainSurfaceHemisphere(rangeNoir , 1 , alphaBrain , 'rh');
        for netIdx = 1:9

            nodeColors = nan(n,3);
            target = find(networksAssignment == netIdx);
            nodeColors(target,:) = color(target , :);

            myMetric = zeros(n,1);
            myMetric(target) = strength(target);

            myMetric = myMetric/max(myMetric);

            hold on;
            % h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_OnlyTarget2(range  , 1 , 1 , 'rh'  , myMetric ,...
            %     cMapArray{1 , netIdx} , target );
            h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaColor2 , 'rh' , nodeColors  ) ;
        end
        view([-90 0])
    end

    axis equal tight off
end


%     exportgraphics(fig , sprintf('UK-biobank/parcellations_code/fig_RSN_surface/RSN_%s.png' , Networks9Cell{netIdx}) , 'Resolution' , 300)
kk = 5;






%% plot scahefer parcellation all in one plot 2 views
alphaColor2 = 1;

% color = myColorMap(networksAssignment , :);

fig = figure;
fig.Position = [2200 607 580 420];

[ha, pos] = tight_subplot(1,2,[.01 .01],[.2 .01],[.01 .01]);

rangeNoir = [0 0.01];

for kView = 1:2
    %         subplot(6,6, subplotIdx(netIdx , kView))
    axes(ha(kView))

    if kView == 1

        h_brain = my_simpleBrainSurfaceHemisphere(rangeNoir , 1 , alphaBrain , 'lh');

        for netIdx = 1:9


            nodeColors = nan(n,3);
            target = find(networksAssignment == netIdx);
            nodeColors(target,:) = color(target , :);

            % networkCmap = cMapArray{1,netIdx};

            % target = find(networksAssignment == netIdx);


            %     myMetric = double(networksAssignment == netIdx);

            myMetric = zeros(n,1);
            myMetric(target) = strength(target);

            myMetric = myMetric/max(myMetric);

            hold on;
            % h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_OnlyTarget2(range  , 1 , 1 , 'lh'  , myMetric ,...
            %     cMapArray{1 , netIdx} , target );

            h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaColor2 , 'lh' , nodeColors  ) ;

        end
        view([-90 0])

    elseif kView == 2
        h_brain = my_simpleBrainSurfaceHemisphere(rangeNoir , 1 , alphaBrain , 'lh');
        for netIdx = 1:9

            nodeColors = nan(n,3);
            target = find(networksAssignment == netIdx);
            nodeColors(target,:) = color(target , :);

            myMetric = zeros(n,1);
            myMetric(target) = strength(target);

            myMetric = myMetric/max(myMetric);

            hold on;
            % h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_OnlyTarget2(range  , 1 , 1 , 'lh'  , myMetric ,...
            %     cMapArray{1 , netIdx} , target );
            h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaColor2 , 'lh' , nodeColors  ) ;

        end
        view([90 0])
    end
    axis equal tight off
end


%     exportgraphics(fig , sprintf('UK-biobank/parcellations_code/fig_RSN_surface/RSN_%s.png' , Networks9Cell{netIdx}) , 'Resolution' , 300)
kk = 5;




%% plot scahefer parcellation all in one view
alphaColor2 = 1;

% color = myColorMap(networksAssignment , :);

fig = figure;
fig.Position = [2200 607 580 420];


rangeNoir = [0 0.01];


% h_brain = my_simpleBrainSurface(rangeNoir , 1 , alphaBrain );
h_brain = my_simpleBrainSurface(range , 1 , alphaBrain );


for netIdx = 1:9
    nodeColors = nan(n,3);
    target = find(networksAssignment == netIdx);
    nodeColors(target,:) = color(target , :);
    myMetric = zeros(n,1);
    myMetric(target) = strength(target);

    myMetric = myMetric/max(myMetric);

    hold on;
    h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaColor2 , 'lh' , nodeColors  ) ;

    % hl = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_OnlyTarget2(range  , 1 , 1 , 'lh'  , myMetric ,...
    %     cMapArray{1 , netIdx} , target );
    hold on;
    % hr = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_OnlyTarget2(range  , 1 , 1 , 'rh'  , myMetric ,...
    %     cMapArray{1 , netIdx} , target );
    h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaColor2 , 'rh' , nodeColors  ) ;

end

view([0 90])


axis equal tight off

kk = 5;

%% plot scahefer parcellation single networks



for netIdx = 1:9

    nodeColors = nan(n,3);
    target = find(networksAssignment == netIdx);
    nodeColors(target,:) = color(target , :);


    fig = figure;
    [ha, pos] = tight_subplot(2,2,[.01 .01],[.2 .01],[.01 .01]);


    alphaBrain = 1;
    alphaColor2 = 1;

    for kView = 1:4
        %         subplot(6,6, subplotIdx(netIdx , kView))
        axes(ha(kView))



        if kView == 1
            h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'lh');
            hold on;
            h1 = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaColor2 , 'lh' , nodeColors  ) ;
            hold on;
            view([-90 0])
        elseif kView == 2
            h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'rh');
            hold on;
            h1 = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaColor2 , 'rh' , nodeColors  ) ;
            hold on;
            view([90 0])
        elseif kView == 3
            h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'lh');
            hold on;
            h1 = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaColor2 , 'lh' , nodeColors  ) ;
            hold on;
            view([90 0])
        elseif kView == 4
            h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'rh');
            hold on;
            h1 = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaColor2 , 'rh' , nodeColors  ) ;
            hold on;
            view([-90 0])
        end

        axis equal tight off
    end
    % exportgraphics(fig , sprintf('UK-biobank/parcellations_code/fig_RSN_surface/RSN_%s_noPattern.png' , Networks9Cell{netIdx}) , 'Resolution' , 300)

end


kk = 5;
