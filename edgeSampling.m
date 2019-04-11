%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clint Olsen
% ECEN 5322: Higher-Dimensional Datasets
% Final Project: Assignment 5 and 6
% Edge Sampling Algorithm 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - edgeSampling.m: Graph subsampling algorithm 2
% - Inputs: Adjacency Matrix A, prob. dist p, num of sub
% sampled edges k
% - Outputs: Subsampled graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function As = edgeSampling(A, p, k)

    % Subsampled Adjacency Matrix
    As = zeros(size(A,1), size(A,2));
    Es = zeros(1, k);
    
    edgeIndex = 0;
    % Generate Collection of Edges (only through upper triangular)
    for i=1:size(A,1)
        for j=i:size(A,1)
            if A(i,j) > 0
                edges(edgeIndex+1).i = i;
                edges(edgeIndex+1).j = j;
                edges(edgeIndex+1).index = edgeIndex;
                edges(edgeIndex+1).weight = A(i,j);
                edgeWeights(edgeIndex+1) = A(i,j);
                edgeIndex = edgeIndex + 1;
            end
        end
    end
    
    if p == "uniform"
        Es = randperm(size(A,1), k);
        
        for i=1:k
            % Place selected Edges and Vertices into Adj. Matrix
            As(edges(Es(i)).i, edges(Es(i)).j) = edges(Es(i)).weight;            
        end
    end

    if p == "weighted"
        edgeWeights = edgeWeights ./ sum(edgeWeights);
        
        for i=1:k
            Es(i) = randsample(1:size(edgeWeights, 2), 1, true, edgeWeights);
            edgeWeights(Es(i)) = 0; % Alters to without replacement
            
            % Place selected Edges and Vertices into Adj. Matrix
            As(edges(Es(i)).i, edges(Es(i)).j) = edges(Es(i)).weight;
        end        
    end
end
% 
% function [Vs, As] = InducedGraphSampling(G, p, k)
% %INDUCED Summary of this function goes here
% %   Detailed explanation goes here
%     Vs = zeros(1, k);
%     
%     %Select K vertices based on P from G(E,V)
%     [x, y] = size(G);
%     
%     if p == 'uniform'
%         %uniform sampling without replacement
%         Vs = randperm(x, k);
%         
%     elseif p == 'proportional'
%         %calculate weights based on vertex degree's
%         Vweights = sum(G);
%         Volume = sum(Vweights);
%         w = Vweights./Volume;
%         
%         %weighted sampling without replacement
%         n = 1:x;
%         Vs = zeros(1,k);
%         for i = 1:k
%             Vs(i) = randsample(n,1,true,w);
%             w(n == J(i)) = 0;
%         end        
%     end
%     
%     %Select Edges 
%             
%     
% end