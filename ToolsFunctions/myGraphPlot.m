function myGraphPlot(n , A , B, Drivers, lambda , symOpt , layout )
if nargin < 7
    layout = 'force';
end

''' -------- Controllability PreAnalysis ---------------''';
nodeNames = {};
for k=1:n
    nodeNames{k} = sprintf('%i',k);
end
if symOpt
    G = graph(A', nodeNames);
else
    G = digraph(A' , nodeNames);
end

supEdge = 0;
numD = sum(Drivers > 0);
for dd = find(Drivers==1)
    actuators = find(B(:,dd));
    if length(actuators) ==1
        G = addedge(G , sprintf('U%i' , dd) ,sprintf('%i',dd) , 1);
        supEdge = supEdge +1;
    else
        for acts = actuators'
            G = addedge(G , sprintf('U%i',dd) , sprintf('%i',acts) , B(acts,dd));
            supEdge = supEdge +1;
        end
    end
end
figure
if n < 30
    pp =plot(G ,'EdgeLabel',G.Edges.Weight , 'MarkerSize' , 10 ...
    , 'NodeFontWeight' , 'bold' , 'NodeFontSize' , 10 ,'Layout',layout);
else
    pp =plot(G , 'MarkerSize' , 10 ...
    , 'NodeFontWeight' , 'bold' , 'NodeFontSize' , 10 ,'Layout',layout);
end
%layout(G,'layered')
%legend('Drivers')
pp.EdgeCData = [zeros(1,G.numedges -supEdge) ones(1,supEdge)];
pp.NodeCData = [zeros(1,G.numnodes - numD) ones(1,numD)];
colormap jet



if sum(find(imag(lambda)))==0
    figure
    plot(lambda , zeros(n,1) , '*' ,'Color','k', 'MarkerSize' , 12)
    xline(0)
    yline(0)
    xlim([min(lambda) 1])
    %         ylim([-2 2])
    title('Sp(A)','FontWeight','bold')
else
    figure
    plot((lambda) , '*' ,'Color','k', 'MarkerSize' , 12)
    xline(0)
    yline(0)
    xlim([-2 2])
    ylim([-2 2])
    title('Sp(A)','FontWeight','bold')
end
end