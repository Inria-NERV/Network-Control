
clear all;


%% add path
subdir = genpath("fig_utils");
addpath(subdir);
subdir = genpath("results_data");
addpath(subdir);

%% load data


ROI_info_Table = readtable( 'results_data/Schaefer200_allinfo.csv');
networksAssignment = table2array(ROI_info_Table(:,5));
nNetworks = max(networksAssignment);

n=214;


Networks9Cell = {'Vis', 'SomMot',  'DorsAttn',  'SalVentAttn',...
    'Limbic', 'Cont','Default', 'TempPar' , 'Sub'};

leftIdx = [1:100 201:207];
rightIdx = [101:200 208:214];


load('results_data/bigTensor_results_target9nets_on_6134_subjects.mat')



alphaBrain = 0.36;
range = [0 1];
%% dealing  with ax indices
subplotIdx = [1 2 7 8;...
    3 4 9 10;...
    5 6 11 12;...
    13 14 19 20;...
    15 16 21 22;...
    17 18 23 24;...
    25 26 31 32;...
    27 28 33 34;...
    29 30 35 36];

leftIdx = [1:100 201:207];
rightIdx = [101:200 208:214];



%% Color bar
cOrd=colororder;
matlabYellow = cOrd(3,:);
matlabRed = cOrd(2,:);

myYellow = rgb("Yellow");

myYellow= 1/255*[252 234 0];

myjet = jet();
myGrayCorrMap = createcolormap(0.2*[1 1 1] ,0.5*[1 1 1], myjet(end - 72,:));

myHot = hot(100);
% myGrayCorrMap = createcolormap(0.6*[1 1 1],0.6*[1 1 1], myHot(60 , :) , myHot(end - 100 , :) , myHot(end - 50 , :));


myGrayCorrMap = [0.75*ones(2,3) ; createcolormap(10,0.65*[1 1 1],myHot(60,:)) ;myHot(60:end-30,:)];
% myGrayCorrMap = createcolormap(0.4*[1 1 1] ,0.5*[1 1 1], myHot(80 , :) , myHot(end - 50 , :));

myGrayCorrMap = [0.65*ones(1,3) ; createcolormap(100,0.65*[1 1 1],myHot(end-30,:)) ;flipud(myHot(60:end-30,:))];

% this is my best
myGrayCorrMap = [createcolormap(60,0.75*[1 1 1],myYellow) ;createcolormap(40,myYellow,matlabRed,rgb('FireBrick'))];



% myGrayCorrMap = [createcolormap(100,0.75*[1 1 1],myHot(end-29,:)) ; flipud(myHot(end-120 : end-30,:))];

myGrayCorrMap = [createcolormap(120,0.75*[1 1 1],myHot(end-29,:)) ; flipud(myHot(end-80 : end-30,:))];

% myGrayCorrMap = [createcolormap(120,0.75*[1 1 1],myHot(end-72,:)) ; (myHot(end-70 : end-10,:))];


% myGrayCorrMap = [createcolormap(80,0.75*[1 1 1],matlabYellow) ];

% myGrayCorrMap = [createcolormap(100,0.65*[1 1 1],myHot(end-50,:)); ];


myHot = hot(200);
myGrayCorrMap = [createcolormap(260,0.75*[1 1 1],myHot(end-57,:)) ; flipud(myHot(end-165 : end-55,:))];


myCMap = myGrayCorrMap;

% myCMap = hot();
% myCMap = pink();
%%%% plot color bar
figure;
imagesc([ -1 : 0.1 : 1]);
cb = colorbar;
colormap(myCMap)
cb.FontSize = 16;



%% plot 9 views  without asterisk %%%%%%%%%%%%%%
coord = table2array(ROI_info_Table(:,7:9));
coord = coord - ones(n,1)* mean(coord);

epsX = 0;
coord(:,1) = coord(:,1)-epsX;

epsY = 20;
coord(:,2) = coord(:,2)-epsY;

epsZ = -15;
coord(:,3) = coord(:,3)-epsZ;


range = [0 0.5];

myFig = figure;
axBound = gca;

fig = figure;
fig.Position(3:4) = [800 718];

alphaBrain =  1;

[ha, pos] = tight_subplot(6,6,[.0008 .01],[.2 .01],[.1 .01]);

% metricInside = zeros(nNetworks, nSubjects);
% metricOutside = zeros(nNetworks, nSubjects);

for netIdx = 1: nNetworks
    % for netIdx = 2:2

    target = find(networksAssignment == netIdx);
    notTarget = setdiff(1:n , target);

    allSubjectsMetric = reshape(lapGramMinBasedTensor(:,netIdx, :), [n , nSubjects]);

    % metricInside(netIdx,:) = mean(allSubjectsMetric(target,:));
    % metricOutside(netIdx,:) = mean(allSubjectsMetric(notTarget,:));

    myMetric = mean(abs(allSubjectsMetric), 2);


    % myMetric = log10(myMetric).^3;

    myMetric = (myMetric).^0.2;

    myMetric = (myMetric -min(myMetric)) / (max(myMetric) -min(myMetric));


    %myMetric = (myMetric).^0.6;


    idxMetric = myMetric;

    idxMetric = 1 + floor((size(myCMap , 1)-1) * idxMetric);

    nodeColors = myCMap(idxMetric , :);


    for kView = 1:4


        %%%%%% deal with range
        if (kView==1 ||kView==3)
            metric2plot = myMetric(leftIdx);
        else
            metric2plot = myMetric(rightIdx);
        end
        range = [0 (max(metric2plot) -min(metric2plot))];

        range = [0 1];

        %%%%%% first get the boundary
        axes(axBound)
        %%%% color just target
        nodeColorsBound = ones(n,3);
        nodeColorsBound(target , :) = zeros(length(target),3);

        % h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'lh');
        % hold on;
        if kView == 1
            hBound = my_simpleBrainSurfaceHemisphere_cData_colorNoMotif_noShade(range  , 1 , 1 , 'lh' , nodeColorsBound  ) ;
        elseif kView == 2
            hBound = my_simpleBrainSurfaceHemisphere_cData_colorNoMotif_noShade(range  , 1 , 1 , 'rh' , nodeColorsBound  ) ;
        elseif kView == 3
            hBound = my_simpleBrainSurfaceHemisphere_cData_colorNoMotif_noShade(range  , 1 , 1 , 'lh' , nodeColorsBound  ) ;
        elseif kView == 4
            hBound = my_simpleBrainSurfaceHemisphere_cData_colorNoMotif_noShade(range  , 1 , 1 , 'rh' , nodeColorsBound  ) ;
        end
        faces = hBound.Faces;
        nFaces = size(faces,1);

        myCdata = hBound.CData;
        numBlack = reshape(sum(myCdata, 1) , [ nFaces 3]);

        myFacesOfInterest = or(numBlack(:,1)==1, numBlack(:,1)==2);

        allverticiesOfInterest = faces(myFacesOfInterest);
        allverticiesOfInterest = reshape(allverticiesOfInterest, [numel(allverticiesOfInterest),1]);
        verticesOfInterest = unique(allverticiesOfInterest);

        %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% then dio the plot
        axes(ha(subplotIdx(netIdx , kView)))

        if kView == 1
            % h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'lh');
            % hold on;
            % h1 = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaBrain , 'lh' , nodeColors  ) ;
            % hold on;
            h1 = my_simpleBrainSurfaceHemisphere_cData_colorNoMotif_rangeHemisph(range  , 1 , alphaBrain , 'lh' , nodeColors  ) ;
            hold on;

            view([-90 0])
            % if netIdx==6
            %     view([-90 20])
            % end

        elseif kView == 2
            % h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'rh');
            % hold on;
            % h1 = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaBrain , 'rh' , nodeColors  ) ;
            % hold on;
            h1 = my_simpleBrainSurfaceHemisphere_cData_colorNoMotif_rangeHemisph(range  , 1 , alphaBrain , 'rh' , nodeColors  ) ;
            hold on;

            view([90 0])
            % if netIdx==6
            %     view([90 -20])
            % end
        elseif kView == 3
            % h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'lh');
            % hold on;
            % h1 = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaBrain , 'lh' , nodeColors  ) ;
            % hold on;

            h1 = my_simpleBrainSurfaceHemisphere_cData_colorNoMotif_rangeHemisph(range  , 1 , alphaBrain , 'lh' , nodeColors  ) ;
            hold on;

            view([90 0])
        elseif kView == 4
            % h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'rh');
            % hold on;
            % h1 = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaBrain , 'rh' , nodeColors  ) ;
            % hold on;
            h1 = my_simpleBrainSurfaceHemisphere_cData_colorNoMotif_rangeHemisph(range  , 1 , alphaBrain , 'rh' , nodeColors  ) ;
            hold on;

            view([-90 0])
        end
        h1.FaceVertexCData(verticesOfInterest , :) = zeros(length(verticesOfInterest),3);

        axis equal tight off

    end

    k = netIdx
end




%% plot 9 views  with best out driver %%%%%%%%%%%%%%
coord = table2array(ROI_info_Table(:,7:9));
coord = coord - ones(n,1)* mean(coord);

epsX = 0;
coord(:,1) = coord(:,1)-epsX;

epsY = 20;
coord(:,2) = coord(:,2)-epsY;

epsZ = -15;
coord(:,3) = coord(:,3)-epsZ;


range = [0 0.5];

myFig = figure;
axBound = gca;

fig = figure;
fig.Position(3:4) = [800 718];

alphaBrain =  1;

[ha, pos] = tight_subplot(6,6,[.0008 .01],[.2 .01],[.1 .01]);

% metricInside = zeros(nNetworks, nSubjects);
% metricOutside = zeros(nNetworks, nSubjects);

for netIdx = 1: nNetworks
% for netIdx = 2:2

    target = find(networksAssignment == netIdx);
    notTarget = setdiff(1:n , target);

    allSubjectsMetric = reshape(lapGramMinBasedTensor(:,netIdx, :), [n , nSubjects]);

    % metricInside(netIdx,:) = mean(allSubjectsMetric(target,:));
    % metricOutside(netIdx,:) = mean(allSubjectsMetric(notTarget,:));

    myMetric = mean(abs(allSubjectsMetric), 2);


    % myMetric = log10(myMetric).^3;

    myMetric = (myMetric).^0.2;

    myMetric = (myMetric -min(myMetric)) / (max(myMetric) -min(myMetric));


    %myMetric = (myMetric).^0.6;


    idxMetric = myMetric;

    idxMetric = 1 + floor((size(myCMap , 1)-1) * idxMetric);

    nodeColors = myCMap(idxMetric , :);


    for kView = 1:4


        %%%%%% deal with range
        if (kView==1 ||kView==3)
            metric2plot = myMetric(leftIdx);
        else
            metric2plot = myMetric(rightIdx);
        end
        range = [0 (max(metric2plot) -min(metric2plot))];

        range = [0 1];

        %%%%%% first get the boundary
        axes(axBound)
        %%%% color just target
        nodeColorsBound = ones(n,3);
        nodeColorsBound(target , :) = zeros(length(target),3);

        % h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'lh');
        % hold on;
        if kView == 1
            hBound = my_simpleBrainSurfaceHemisphere_cData_colorNoMotif_noShade(range  , 1 , 1 , 'lh' , nodeColorsBound  ) ;
        elseif kView == 2
            hBound = my_simpleBrainSurfaceHemisphere_cData_colorNoMotif_noShade(range  , 1 , 1 , 'rh' , nodeColorsBound  ) ;
        elseif kView == 3
            hBound = my_simpleBrainSurfaceHemisphere_cData_colorNoMotif_noShade(range  , 1 , 1 , 'lh' , nodeColorsBound  ) ;
        elseif kView == 4
            hBound = my_simpleBrainSurfaceHemisphere_cData_colorNoMotif_noShade(range  , 1 , 1 , 'rh' , nodeColorsBound  ) ;
        end
        faces = hBound.Faces;
        nFaces = size(faces,1);

        myCdata = hBound.CData;
        numBlack = reshape(sum(myCdata, 1) , [ nFaces 3]);

        myFacesOfInterest = or(numBlack(:,1)==1, numBlack(:,1)==2);

        allverticiesOfInterest = faces(myFacesOfInterest);
        allverticiesOfInterest = reshape(allverticiesOfInterest, [numel(allverticiesOfInterest),1]);
        verticesOfInterest = unique(allverticiesOfInterest);

        %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% then dio the plot
        axes(ha(subplotIdx(netIdx , kView)))

        if kView == 1
            % h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'lh');
            % hold on;
            % h1 = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaBrain , 'lh' , nodeColors  ) ;
            % hold on;
            h1 = my_simpleBrainSurfaceHemisphere_cData_colorNoMotif_rangeHemisph(range  , 1 , alphaBrain , 'lh' , nodeColors  ) ;
            hold on;

            view([-90 0])
            % if netIdx==6
            %     view([-90 20])
            % end

        elseif kView == 2
            % h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'rh');
            % hold on;
            % h1 = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaBrain , 'rh' , nodeColors  ) ;
            % hold on;
            h1 = my_simpleBrainSurfaceHemisphere_cData_colorNoMotif_rangeHemisph(range  , 1 , alphaBrain , 'rh' , nodeColors  ) ;
            hold on;

            view([90 0])
            % if netIdx==6
            %     view([90 -20])
            % end
        elseif kView == 3
            % h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'lh');
            % hold on;
            % h1 = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaBrain , 'lh' , nodeColors  ) ;
            % hold on;

            h1 = my_simpleBrainSurfaceHemisphere_cData_colorNoMotif_rangeHemisph(range  , 1 , alphaBrain , 'lh' , nodeColors  ) ;
            hold on;

            view([90 0])
        elseif kView == 4
            % h_brain = my_simpleBrainSurfaceHemisphere(range , 1 , alphaBrain , 'rh');
            % hold on;
            % h1 = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_colorNoMotif(range  , 1 , alphaBrain , 'rh' , nodeColors  ) ;
            % hold on;
            h1 = my_simpleBrainSurfaceHemisphere_cData_colorNoMotif_rangeHemisph(range  , 1 , alphaBrain , 'rh' , nodeColors  ) ;
            hold on;

            view([-90 0])
        end
        h1.FaceVertexCData(verticesOfInterest , :) = zeros(length(verticesOfInterest),3);

        axis equal tight off

    
% end
    %%%%%%%%%%%%%%
    %%% get best driver
% for netIdx = 1 : 9

    myMetricJustOut=myMetric;
    myMetricJustOut(target)=0;
    myMetricJustOut(201:214)=0;

    [ddbest , iibest] = max(myMetricJustOut)

    inLeft = isempty(intersect(rightIdx , iibest));

    mkSize = 50;
    mkSize2 = 31;

    offSetPlot = 60;


    ftSize = 42;
    myText = '*';

    view1 = [0 1 0 0 0 1 1 0 0];
    view2 = [0 0 1 0 0 0 0 1 1];
    view3 = [0 0 0 1 0 0 0 0 0];
    view4 = [1 0 0 0 1 0 0 0 0];
    %%%%%%%%%%%%%% Same but with text
    if inLeft

        % view1
        if view1(netIdx) && kView==1
            axes(ha(subplotIdx(netIdx , kView)))
            hold on;
            % text(coord(iibest,1)-offSetPlot ,coord(iibest,2),coord(iibest,3),myText,'FontSize',ftSize)
            % hold on;
            % text(coord(iibest,1)-offSetPlot ,coord(iibest,2),coord(iibest,3),'*')
            plot3(coord(iibest,1)-offSetPlot ,coord(iibest,2),coord(iibest,3) , 'Marker','.' , 'Color',  'k', 'MarkerSize' , mkSize)
            hold on;
            plot3(coord(iibest,1)-offSetPlot ,coord(iibest,2),coord(iibest,3) , 'Marker','.' , 'Color',  'w', 'MarkerSize' , mkSize2)
hold on;
        end

                % view3
        if view3(netIdx) && kView==3
            axes(ha(subplotIdx(netIdx , kView)))
            hold on;
            % text(coord(iibest,1)-offSetPlot ,coord(iibest,2),coord(iibest,3),myText,'FontSize',ftSize)
            % hold on;
            % text(coord(iibest,1)-offSetPlot ,coord(iibest,2),coord(iibest,3),'*')
            plot3(coord(iibest,1)+offSetPlot ,coord(iibest,2),coord(iibest,3) , 'Marker','.' , 'Color',  'k', 'MarkerSize' , mkSize)
            hold on;
            plot3(coord(iibest,1)+offSetPlot ,coord(iibest,2),coord(iibest,3) , 'Marker','.' , 'Color',  'w', 'MarkerSize' , mkSize2)
hold on;
        end

    else

                % view2
        if view2(netIdx) && kView==2
            axes(ha(subplotIdx(netIdx , kView)))
            hold on;
            % text(coord(iibest,1)-offSetPlot ,coord(iibest,2),coord(iibest,3),myText,'FontSize',ftSize)
            % hold on;
            % text(coord(iibest,1)-offSetPlot ,coord(iibest,2),coord(iibest,3),'*')
            plot3(coord(iibest,1)+offSetPlot ,coord(iibest,2),coord(iibest,3) , 'Marker','.' , 'Color',  'k', 'MarkerSize' , mkSize)
            hold on;
            plot3(coord(iibest,1)+offSetPlot ,coord(iibest,2),coord(iibest,3) , 'Marker','.' , 'Color',  'w', 'MarkerSize' , mkSize2)
hold on;
        end

        % view4
        if view4(netIdx) && kView==4
            axes(ha(subplotIdx(netIdx , kView)))
            hold on;
            % text(coord(iibest,1)-offSetPlot ,coord(iibest,2),coord(iibest,3),myText,'FontSize',ftSize)
            % hold on;
            % text(coord(iibest,1)-offSetPlot ,coord(iibest,2),coord(iibest,3),'*')
            plot3(coord(iibest,1)-offSetPlot ,coord(iibest,2),coord(iibest,3) , 'Marker','.' , 'Color',  'k', 'MarkerSize' , mkSize)
            hold on;
            plot3(coord(iibest,1)-offSetPlot ,coord(iibest,2),coord(iibest,3) , 'Marker','.' , 'Color',  'w', 'MarkerSize' , mkSize2)
hold on;
        end
    end
    
    end
    k = netIdx
end





