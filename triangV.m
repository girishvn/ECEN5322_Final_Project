%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clint Olsen
% ECEN 5322: Higher-Dimensional Datasets
% Final Project: Assignment 3 and 4
% Clustering Coefficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - numTriang.m: Determines the number of triangle clusters
% a given vertex is in
% - Inputs: Adjacency Matrix A and vertex v
% - Outputs: Number of triangle v is a part of
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function numTriang = triangV(A, v)

    numTriang = 0;
    for j=1:size(A,1) % Origin node
        if A(v,j) > 0 % Loop until neighbor is found
            for k=1:size(A,1)
               if A(j,k) > 0 && k ~= v % Neighbor of neighbor found
                   if A(k,v) > 0 % Neighbor of neighbor is connected to origin
                       numTriang = numTriang + 1; % Indicate v is part of a triangle
                   end
               end
            end
        end
    end
    
    numTriang = (numTriang / 2);              
end