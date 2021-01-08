function [Nd,drivernodes,Nconfigs,ld,C] = myExactControllability(A,varargin)
% Computes the exact controllability of matrix A of arbitrary structure. 
% Input: adjacency matrix A, where entry in row i, column j, indicates a
% connection from i to j. 
% Use second argument (..,'plotting',0) to disable generating network figure.
% Default is to plot.
% Returns Nd = number of driver nodes, drivernodes = list of driver node IDs, 
% Nconfigs = number of different configurations that satisfy the requirement of Nd. 
% Each time the code is run, a different set of configuration of driver nodes 
% will be returned. The total number of driver nodes remains constant.
% Ref: Exact Controllability of complex networks by Yuan et al, Nat Comm 2013

% Code written by:
% Tapan P Patel, PhD, tapan.p.patel@gmail.com
% University of Pennsylvania
% Please cite: Patel TP et al, Automated quantification of neuronal 
% networks and single-cell calcium dynamics using calcium imaging.
% J Neurosci Methods, 2015.




% Examples 1, 2 and 3 are three different network configurations illustrated
% in Figure 1 of Liu et al, Controllability of complex networks, Nature
% 2011

% Example 1:
% Directed path, 4 nodes. Node 1 connects to node 2, node 2 connects to
% node 3, node 3 connects to node 4
% A = [0 1 0 0;
%      0 0 1 0;
%      0 0 0 1;
%      0 0 0 0];
% [Nd, drivernodes,Nconfig] = ExactControllability(A,'plotting',1);
% The above command returns Nd = 1, i.e. total of 1 driver node,
% drivernodes = 1, i.e. node #1 is the driver node
% Nconfig = 1, i.e. there is only 1 configuration of driver node

% Example 2:
% Directed star. 4 nodes, where node 1 connects to all other 3 nodes.
% A = [0 1 1 1;
%     0 0 0 0;
%     0 0 0 0;
%     0 0 0 0;];
% [Nd, drivernodes,Nconfig] = ExactControllability(A,'plotting',1);
% In this example, there are total of 3 driver nodes, Nd = 3.
% The identity of driver nodes are provided in the vector drivernodes.
% There are total of 3 distinct configurations of selecting drivernodes,
% Nconfig = 3. Each time the algorithm is run, a different set of
% drivernodes will be returned. The possibilites are: {1,3,4}, {1,2,3},
% {1,2,4}

% Example 3:
% A = [0 1 1 1 1 1;
%     0 0 0 0 0 0;
%     0 0 0 0 0 0;
%     0 0 0 0 0 0;
%     0 0 0 0 0 0;
%     0 0 0 0 1 0;];
% [Nd, drivernodes,Nconfig] = ExactControllability(A,'plotting',1);
% In this example, there are Nd = 4 driver nodes and 4 distinct
% configurations. The possiblities of drivernodes = {1,2,3,4}, {1,3,4,6},
% {1,2,3,6}, {1,2,4,6}

% Example 4:
% Let's try a more complex network. Let's build a scale-free directed
% network with 300 nodes, using a seed network size of 10 nodes and average
% connectivity 10. Please see the accompanying BAgraph_dir.m file for more
% info on the algorithm. The BAgraph_dir algorithm may take some time to
% generate an adjacency matrix. To save time, you may load a sample
% adjacency matrix by executing: load('scalefree_adjacency.mat');
% A = BAgraph_dir(300,10,10);
% A scale-free network should have a power-law degree distribution. Verify
% by 
% figure; hist(sum(A),20); xlabel('Degree'); ylabel('Count');
% [Nd, drivernodes,Nconfig] = ExactControllability(A,'plotting',1);
% There are Nd = 33 driver nodes and only 1 configuration. The identity of
% driver nodes is provided in the variable drivernodes.
% You will note that the drivernodes have very low in-degrees
% sum(A(:,drivernodes)

% Example 5:
% Lastly, let's try a random undirected symmetric Erdos Reyni network of 300 nodes with 0.01
% probability of connection between any pair-wise nodes
% A = erdos_reyni(300,0.01);
% [Nd, drivernodes,Nconfig] = ExactControllability(full(A),'plotting',1);
% There are Nd = 22 driver nodes and 16 different ways to choose 22 driver
% nodes. Note: as the network becomes less sparse, i.e. as the probability
% of connection incresease from 0.01, the number of driver nodes decreases.


%%

% In Yuan's algorithm, a(i,j) refers to a connection from j to i rather than 
% a connection from i to j. So the input needs to be transposed.

% Multiplicity of an eigenvalue i is the number of times lambda_i is repeated 
% (tolerance 1e-8). lambda_i that corresponds to the maximum multiplicity is 
% then used to compute a geomtric multiplicity.
% C is the column reduced matrix whose dependent rows  are the identities of 
% the driver nodes. There can be more than 1 set of driver nodes consistent 
% with the requirement N_D total driver nodes.

plotting = 0;
% ------------------------------------------------------------------------------
% parse varargin
% ------------------------------------------------------------------------------
if nargin > 2
    for i = 1:2:length(varargin)-1
        if isnumeric(varargin{i+1})
            eval([varargin{i} '= [' num2str(varargin{i+1}) '];']);
        else
            eval([varargin{i} '=' char(39) varargin{i+1} char(39) ';']);
        end
    end
end
%A = double(A);
%A = A';
lambda = eig(A);
%multiplicity = -1;
%lambda_M = 0;
%% find maximum multiplicity eigenvalue

%%%%% This is how it was done wrong %%M%%%%%%%
% for i=1:length(lambda)
% 	x = length(find(abs(lambda(i)-lambda)<1e-8));
% 	if(~isempty(x) && x>multiplicity)
% 		multiplicity = x;
% 		lambda_M = lambda(i);
% 	end
% end
N = size(A,1);
geoMul = zeros(N,1);
lamAlreadySeen = NaN(N,1);

for i=1:N
    lam = lambda(i);
	if sum(abs(lamAlreadySeen - lam*ones(N,1))< 10^-8) == 0
        lamAlreadySeen(i) = lam;
		geoMul(i) = N - rank(lam*eye(N)-A);
	end
end

[mu_M,I] = max(geoMul);
maxIndex = find(geoMul == mu_M);
numMaxGeoMul = length(maxIndex);

lambda_M = lamAlreadySeen(maxIndex(end));
Nd = mu_M;

%%
% mu_M is the number of minimum driver nodes
% lambda_M is used to identify the set of driver nodes
M = A- lambda_M*eye(N);
[R,p , perms] = myrref(M');
C=R';
%% Linearly dependent rows of C are the driver nodes. However, 
% there can be more than 1 configuration of the set of driver nodes. 
%Let's enumerate the possibilities. Row i is linearly independent 
%if C(i,i) ~=0 and all other entries in that row are 0. Similarly row i 
%is linearly dependent if either all rows are 0 or there are more 
%than 1 nonzero entries in the row. Configurations of driver nodes
%would include all the dependent rows and a node from each of the 
%independent row if the multiplicity at row i is >1.

li = cell(N,1);
ld = [];
for i=1:N
	if(nnz(C(i,:))==1)
		idx = find(C(i,:));
		x = li{idx};
		x = [x i];
		li{idx} = x;
	else
		ld = [ld i];
	end
end	
% Driver node configurations:
DriverConfigs = cell(N,1);
Nconfigs = 1;
for i=1:N
	if(length(li{i})>1)
	% Keep 1 and the rest are dependent rows which make them drivers
	idx = li{i};
	DriverConfigs(i) = {nchoosek(idx,length(idx)-1)};
	% Number of possible configurations is the product of the size of each DriverConfigs{i}
	Nconfigs = Nconfigs*size(DriverConfigs{i},1);
	end
end
% Give a permutation of driver node configuration.
drivernodes = ld;
for i=1:N
	if(~isempty(DriverConfigs{i}))
		x = DriverConfigs{i};
		k = randsample(size(x,1),1);
		drivernodes = [drivernodes x(k,:)];
	end
end

%% Plot
% if(plotting)
% xy = rand(size(A,1),2);
% gplotdc(A,xy);
% % Color driver nodes in blue
% hold on
% plot(xy(drivernodes,1),xy(drivernodes,2),'bo','MarkerFaceColor','b');
% title(sprintf('%d driver nodes identified\nNodes in blue are driver nodes.',Nd));
end
