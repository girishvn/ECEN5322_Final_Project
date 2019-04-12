%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Girish Narayanswamy
% ECEN 5322: Higher-Dimensional Datasets
% Final Project: Assignment 5 and 6
% Edge Sampling Algorithm 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - metropolisHastingsRW.m: Graph subsampling algorithm 3
% - Inputs: Adjacency Matrix A,  num of sub sampled edges ms
% - Outputs: Subsampled graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [As] = metropolisHastingsRW(A,ms)

ms = round(ms);

% Construct sampled A matrix 
[m,n] = size(A); 
As = zeros(m,n); % zero init As adjacency mtx of the subset 
sizeS = 0; % number of edges sampled

v = randi([1 m]); % Random starting vertex
while sizeS < ms
    
    Nv = find(A(v,:) > 0); % neighbor hood of v
    Nv = intersect(Nv,find(As(v,:) == 0)); % Nv limited to unvisited nodes
    if isempty(Nv)
        v = randi([1 m]); % Random starting vertex
        continue;
    end
    
    w = Nv(randi([1 length(Nv)])); % choose a random neighbor of v
    
    p = unifrnd(0,1); % choose uniform-random probability [0 1] 
    degv = sum(A(v,:));
    degw = sum(A(w,:));
    if p <= degv/degw
        As(v,w) = A(v,w); % add connection to new adjacency matric 
        As(w,v) = A(w,v);
        
        v = w; % move to new node
        sizeS = sizeS + 1; % add an edge
    end  
end
end

