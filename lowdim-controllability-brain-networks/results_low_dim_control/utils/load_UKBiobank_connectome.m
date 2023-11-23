function matrix = load_UKBiobank_connectome(list_of_edges_filename , matrixSize)
edgeList = table2array(readtable(list_of_edges_filename));
i = edgeList(:,2) +1;
j = edgeList(:,3) +1;
w = edgeList(:,4);
% build matrix
matrix = full(sparse(i,j,w , matrixSize  , matrixSize));

% add upper part
matrix = matrix+ matrix';


%% get only real parcels
% additional parcel created by multiAtlasTransferTool
% ignore FreeSurfer medial Wall 
% lh.Background+FreeSurfer_Defined_Medial_Wall.label

realIdx = [(2 : 100) (102 :216)]; % this valid for Scaefer200 % change as needed

matrix = matrix(realIdx , realIdx);

end