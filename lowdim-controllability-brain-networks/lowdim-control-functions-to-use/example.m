%% example script controllability

%% load connectome
tableName = ...
'0000000_connectome_edgelist_schaefer200-yeo17_rmap_s1dild.csv';

matrixSize = 216;
edgeList = table2array(readtable(tableName));
i = edgeList(:,2) +1;
j = edgeList(:,3) +1;
w = edgeList(:,4);
% build matrix
matrix = full(sparse(i,j,w , matrixSize  , matrixSize));

% add upper part
matrix = matrix+ matrix';

%%% get only real parcels
% additional parcel created by multiAtlasTransferTool
% ignore FreeSurfer medial Wall 
% lh.Background+FreeSurfer_Defined_Medial_Wall.label

realIdx = [(2 : 100) (102 :216)]; % this valid for Scaefer200 % change as needed

matrix = matrix(realIdx , realIdx);
n = size(matrix,2);
%%%%%%% the resulting network is undirected with zero diagonal

%%%% example replace connectome with random graph 
matrix = rand(n,n); % very dense network!!


%%% !!!!!!! One can load a network with all different ways possible as long
%%% as the variable "matrix" is a matlab array of size nxn

%% calculate control centrality metric for a specific target 
target = 1:12; % correspond to left Visual network in schaefer

[standardCC , lowdimCC]=...
low_dimensional_control_centrality(matrix , target);


%% calculate the capacity of each driver to achieve a specific state transfer
% define states:

x0 = zeros(n,1);x0([1:12 101:112] , 1) = 1; % VIS active
xf = zeros(n,1);xf([13:28 113:140] , 1) = 1; % SMN active


[tarnsferCapacity , lowdim_transferCapacity] = ...
    controlMetric_forStateTransfer(matrix , x0 , xf);




