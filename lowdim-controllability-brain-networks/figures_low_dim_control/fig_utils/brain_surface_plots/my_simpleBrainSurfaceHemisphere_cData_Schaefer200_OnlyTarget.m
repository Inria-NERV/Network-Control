function h = my_simpleBrainSurfaceHemisphere_cData_Schaefer200_OnlyTarget(range  , centerPatch , alphaBrain , hemisphere , cData ,...
    mycolormapTarget , target )
% simpleBrainSurface
% Simple function to render a brain in 3D in MNI coordinates.
% Addapted from surfPlot of the MRtools collection from:
% http://mrtools.mgh.harvard.edu/index.php/Main_Page
%
% Initial surface data comes from FreeSurfer template which apparently comes from:
% Fischl, B., Sereno, M. I., Tootell, R. B.H. and Dale, A. M. (1999),
% High-resolution intersubject averaging and a coordinate system for the
% cortical surface. Hum. Brain Mapp., 8: 272284.
% doi: 10.1002/(SICI)1097-0193(1999)8:4<272::AID-HBM10>3.0.CO;2-4
%
% INPUT:
% par = specs object and specification object
% range = 2D vector graycolor range [darkest brightest]
%           { default [0.1 0.7] }
% OUTPUT:
% h = handle to the surface patch.
% SIDEEFFECTS:
% A brain surface plot is generated.
%
%
% Original:
%%% Written by Aaron P. Schultz - aschultz@martinos.org
%%%
%%% Copyright (C) 2014,  Aaron P. Schultz
%%%
%%% Supported in part by the NIH funded Harvard Aging Brain Study (P01AG036694) and NIH R01-AG027435
%%%
%%% This program is free software: you can redistribute it and/or modify
%%% it under the terms of the GNU General Public License as published by
%%% the Free Software Foundation, either version 3 of the License, or
%%% any later version.
%%%
%%% This program is distributed in the hope that it will be useful,
%%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%% GNU General Public License for more details.
%%%
% Copyright (C) 2015, Robert Rein
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    range = [0.1 0.7];
end



% get vertex data
load('simple_brain_surface.mat');


% get hemisphere data
load(['simple_brain_surface_hemisphere_' hemisphere '.mat']);

hemisphereVerticesIdx = hemisphereSurface.hemisphereVerticesIdx;
hemisphereFaces = hemisphereSurface.hemisphereFaces;
hemisphereVertices = brain.vertices(hemisphereVerticesIdx , :);

%% center the vertices
if centerPatch
%     brain.vertices = brain.vertices - ones(size(brain.vertices , 1) ,1)* mean(brain.vertices);
    hemisphereVertices = hemisphereVertices- ones(size(hemisphereVertices , 1) ,1)* mean(brain.vertices);
end
% figure background color
set(gcf,'renderer','opengl');


tmp_shading_colorHem = brain.shading_pre(hemisphereVerticesIdx) * diff(range);
tmp_shading_colorHem = tmp_shading_colorHem - min(tmp_shading_colorHem) + range(1);
shading_colorHem = repmat(tmp_shading_colorHem,1,3);


%% deal with color
if strcmp(hemisphere , 'lh')
load('mapVertex2parcel_left.mat')
end
if strcmp(hemisphere , 'rh')
load('mapVertex2parcel_right.mat')
end

% cData is a 1*214 vector 
cDataNorm = (cData - min(cData))/(max(cData) - min(cData));

nVertices = size(hemisphereVertices , 1);

newVert = zeros(1 , nVertices);
vertexShowInd = zeros(nVertices ,1);
count = 0;
for kVert = 1 : nVertices
    node = vert2parcelMap(kVert);

    if ~isempty(find(target == node))
        count = count+1;
        vertexShowInd(kVert,: )= 1;
        newVert(kVert) = count;
    end
%     newVert(kVert) = count;
end
vertexShowInd = logical(vertexShowInd);

vertexColor = zeros(nVertices ,3);
for kVert = 1 : nVertices
    node = vert2parcelMap(kVert);
    if ~isempty(find(target == node))
        vertexColor(kVert,: )= mycolormapTarget( 1+ floor(cDataNorm(node)*255), :);
    else
        vertexColor(kVert ,:)= mycolormapTarget( 1+ floor(cDataNorm(node)*255), :);
    end
end

nFaces = size(hemisphereFaces,1);
newFaces = [];
for kFace = 1 : nVertices
        if vertexShowInd(hemisphereFaces(kFace,1)) && vertexShowInd(hemisphereFaces(kFace,2)) && vertexShowInd(hemisphereFaces(kFace,3))

%     if vertexShowInd(hemisphereFaces(kFace,1)) || vertexShowInd(hemisphereFaces(kFace,2)) || vertexShowInd(hemisphereFaces(kFace,3))


        %     if sum(double([vertexShowInd(hemisphereFaces(kFace,1)) ...
        %             vertexShowInd(hemisphereFaces(kFace,2)) ...
        %             vertexShowInd(hemisphereFaces(kFace,3))])) >1
        newFaces = [newFaces;newVert(hemisphereFaces(kFace,:))];

    end
end

newVertVect = hemisphereVertices(vertexShowInd , :);
newVertFind = find(newVert);


%% deal with wholes
% find the faces that are in between the big perimeter

%% plot patch
% h = patch('vertices',hemisphereVertices(vertexShowInd , :), ...
%     'faces', newVert(hemisphereFaces(vertexShowInd,:)) , ...
%     'FaceVertexCdata',vertexColor(vertexShowInd,:));
h = patch('vertices',newVertVect, ...
    'faces', newFaces , ...
    'FaceVertexCdata',vertexColor(vertexShowInd,:));

set(h,'edgecolor','interp','facecolor','interp');
set(gca,'dataaspectratio',ones(1,3),'visible','off');
axis tight;
view(159,6);
alpha(alphaBrain)


h.FaceLighting = 'gouraud';
% h.FaceLighting = 'flat';
shading interp



