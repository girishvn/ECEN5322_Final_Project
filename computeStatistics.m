%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clint Olsen
% ECEN 5322: Higher-Dimensional Datasets
% Final Project: Assignment 3 and 4
% Adjacency Matrix Statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - computeStatistics.m: Computes stats for passed Adj. Matrix
% - Inputs: Adjacency Matrix A
% - Outputs: Displays Table of stat values, vector of degrees,
% and vector of clusteringCoeffs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [degreeDistribution, clusteringCoefficient] = computeStatistics(adjacencyMatrix)

    % Statistics to gather:
    graphSize = 0;
    numEdges = 0;
    volume = 0;
    density = 0;
    degreeDistribution = zeros(1,size(adjacencyMatrix,1));
    averageDegree = 0;
    clusteringCoefficient = zeros(1,size(adjacencyMatrix,1));
    avgCC = 0;
    
    % Size
    graphSize = size(adjacencyMatrix,1);
    
    %Total number of edges
    numEdges = size(find(triu(adjacencyMatrix) > 0), 1);
    
    %Volume (sum of all edge weights)
    volume = sum(sum(triu(adjacencyMatrix)));
    
    %Density (2m/(n(n-1))) where m is number of edges and n is num of nodes
    density = (2*numEdges)/(graphSize*(graphSize-1));
    
    %Degree (sum of all edge weights connected to vertex v)
    for i=1:size(adjacencyMatrix,1)
        degreeDistribution(i) = sum(adjacencyMatrix(i,:));
    end
    
    %Average Degree
    averageDegree = sum(degreeDistribution) / graphSize;
    
    %Clustering Coefficient
    for i=1:size(adjacencyMatrix,1)
        deltaV = triangV(adjacencyMatrix, i);
        if degreeDistribution(i) < 2
            clusteringCoefficient(i) = 0;
        else
            clusteringCoefficient(i) = (2*deltaV)/(degreeDistribution(i)*(degreeDistribution(i)-1));
        end
    end
    
    % Average Clustering Coefficient
    avgCC = sum(clusteringCoefficient) / graphSize;
    
    % Table output of data
    table(graphSize, numEdges, volume, density, averageDegree, avgCC)
end