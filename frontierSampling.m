%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Girish Narayanswamy
% ECEN 5322: Higher-Dimensional Datasets
% Final Project: Assignment 5 and 6
% Frontier Sampling Algorithm 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - frontierSampling.m: Graph subsampling algorithm 4
% - Inputs: Adjacency Matrix A,  num of sub sampled edges ms
% - Outputs: Subsampled graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [As] = frontierSampling(A, ms)

ms = round(ms);

% Construct sampled A matrix 
[m,n] = size(A); 
As = zeros(m,n); % zero init As adjacency mtx of the subset 

s = randi(m);
L = randperm(m,s); % random vertex seed of size s
k = 0; % iterator

while k < ms
    
    % check to see if L must be reshuffled again
    for i = 1:length(L)
        u_test = L(i);
        Nu_test = find(A(u_test,:) > 0);
        Nu_test = intersect(Nu_test,find(As(u_test,:) == 0));
        
        if ~isempty(Nu_test) 
            break % possible accessible nodes still exist
        end
        
        if i == length(L) % no more possible nodes to visit, reseed L
            L = randperm(m,s);
        end
    end
    
    degLTotal = sum(sum(A(:,L))); % total degree weights of u in L
    pu = sum(A(:,L))/degLTotal; % probability of each vertex u
    
    if isempty(find(pu) > 0)
        continue; % if no new edges to add repeat iteration and generate new u
    end
    
    idx = randsample(s,1,true,pu); % choose a random index in L with probabilities pu
    u = L(idx); % find node number using idx generated above
    
    Nu = find(A(u,:) > 0); % neighbor hood of u
    Nu = intersect(Nu,find(As(u,:) == 0)); % Nu limited to unvisited nodes
    if isempty(Nu)
        continue; % if no new edges to add repeat iteration and generate new u
    end 
    
    v = Nu(randi(length(Nu)));
    
    As(u,v) = A(u,v); % add edge to As matrix
    As(v,u) = A(v,u);
    
    k = k + 1; % iterate number of edges
    L(idx) = v; % replace u with outgoing vertex v 
end
end

