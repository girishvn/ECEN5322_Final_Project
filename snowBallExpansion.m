%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Girish Narayanswamy
% ECEN 5322: Higher-Dimensional Datasets
% Final Project: Assignment 5 and 6
% Snow Ball Expansion Algorithm 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - snowBallExpansion.m: Graph subsampling algorithm 5
% - Inputs: Adjacency Matrix A,  num of sub sampled vertices ms
% - Outputs: Subsampled graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [As] = snowBallExpansion(A,ms)

ms = round(ms);

% Construct sampled A matrix 
[m,n] = size(A); 
As = zeros(m,n); % zero init As adjacency mtx of the subset 

v = randi(m); % choose radnom starting vertex
S = v; % contained vertices
sizeS = 1; % number of contained vertices


while sizeS < ms
    
    diffMax = 0; % value to maximize to add new node
    vMax = 0;
    
    NS = mod(find(A(S,:)' > 0),m); % neighborhood of N
    NS(NS == 0) = m; % make the indexes work still
    NSUS = union(NS,S); % get union of N(S) and S
    
    for i = 1:length(NS) % iterate through neighborhood of S
        
        v = NS(i); % vertex in S
        Nv = find(A(v,:) > 0); % neighborhood of v
        
        diffNew = length(setdiff(Nv,NSUS)); % difference between v neighborhood and current set Nhood
        if diffNew > diffMax % if higher difference save new node v
            diffMax = diffNew;
            vMax = v; 
        end
    end
    
    % if there exists a node with new vertices
    if vMax > 0
        % populate edges for new v with existing nodes in S
        for i = 1:length(S)
            u = S(i);
            As(vMax,u) = A(vMax,u);
            As(u,vMax) = A(u,vMax);
        end

       S = [S, vMax]; % add vertex to vertex list
       sizeS = sizeS + 1;
    
    else % no new vertices in neighborhood cluster
        
        % find a new random unvisited node
        all_v = 1:m;
        remaining_v = setdiff(all_v, S);
        idx = randi(length(remaining_v));
        vMax = remaining_v(idx);

        % populate edges for new v with existing nodes in S
        for i = 1:length(S)
            u = S(i);
            As(vMax,u) = A(vMax,u);
            As(u,vMax) = A(u,vMax);
        end

        S = [S, vMax]; % add vertex to vertex list
        sizeS = sizeS + 1;
    end
end
end

