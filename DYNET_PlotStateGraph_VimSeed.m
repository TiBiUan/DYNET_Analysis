%% This function plots, nicely, the information from my TP matrix
% A is the transition probability matrix: it should have size K + 2, with
% the first row/column the baseline state and the last one the unassigned
% state
function [] = DYNET_PlotStateGraph_VimSeed(FC,alpha,Label,CM,Network,Hemisphere)

    A = zeros(length(FC)+1,length(FC)+1);
    A(1,2:end) = FC;

    A_bis = A;

    T = prctile(abs(FC),alpha);

    A(abs(A) < T) = 0;
    
    % Constructs my directional graph
    G = graph((A+A')/2);
    
    % Number of regions
    n_ROIs = size(A,1);
    
    % Plots the initial version of my graph
    figure;
    p2 = plot(G,'NodeFontSize',16,'NodeLabelColor',[0.3,0.7,0.9],'NodeFontWeight','bold');

    % How much do I increase the angle by every time
    Phi = 2*pi/n_ROIs;

    p2.XData(1) = 0;
    p2.YData(1) = 0;
    
    for c = 1:n_ROIs-1
        p2.XData(c+1) = 2*cos(Phi*(c+1));
        p2.YData(c+1) = 2*sin(Phi*(c+1));
    end

    MS = 20*sum(abs(A_bis))';
    MS(1) = 20;
    
    p2.MarkerSize = MS;

    NC(1,:) = [0,0,0];
    
    for r = 1:n_ROIs-1
        NC(r+1,:) = CM(Network(r),:);
    end
    
    p2.NodeColor = NC;
    
    M{1} = 'o';
    
    for r = 2:n_ROIs
        if Hemisphere(r-1) == 1
            M{r} = 'square';
        elseif Hemisphere(r-1) == 2
            M{r} = 'diamond';
        else
            M{r} = 'pentagram';
        end
    end
    
    p2.Marker = M;

    
    % I want my edges to have their color and thickness proportional to 
    % their weight 
    n_edges = size(G.Edges.EndNodes,1);
    
    Weights = G.Edges.Weight;
    
    CMB = cbrewer('seq','Blues',10000);
    CMR = cbrewer('seq','Reds',10000);
    
    Weights_Signs = sign(Weights);
    Weights = abs(Weights);
    Weights = (Weights)/max(Weights);

    % Red or blue depending on the sign
    for e = 1:n_edges
        if Weights_Signs(e) >= 0
            tmp_col2(e,:) = CMR(floor(Weights(e)*10000),:);
        elseif Weights_Signs(e) < 0
            tmp_col2(e,:) = CMB(floor(Weights(e)*10000),:);
        end
    end
    
    % I want edge color proportional to the weight
    p2.EdgeColor = tmp_col2;

    p2.LineWidth = realmin+6*Weights;%2 + 10*Weights;
    
    p2.NodeLabel = Label;
    p2.EdgeAlpha = 0.5;

    axis off
    set(gcf,'color','w');
    set(gca,'Box','off');
end