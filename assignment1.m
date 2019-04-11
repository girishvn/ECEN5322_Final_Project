%%% Questions %%%
% How to compare orig and subsampled?
% Subsampled adj matrix sparse but same size?
% What should log-log plots look like?


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clint Olsen
% ECEN 5322: Higher-Dimensional Datasets
% Final Project: Assignment 1 and 2
% Temporally Aggregated Adjacency Matrix Generation
% Creates and saves adjacency matrix from list of 
% datasets in face2faceData and coPresenceData vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lyon Hospital and 2009 French Society for Hospital Hygiene Conf. Data
face2faceData = ["tij_LH10.dat", "tij_SFHH.dat"];
coPresenceData = ["tij_pres_LH10.dat", "tij_pres_SFHH.dat"];
metaData = ["sorted_indices_LH10.mat", "sorted_indices_SFHH.mat"];
fAdjOutput = ["f2f_adj_LH10.mat", "f2f_adj_SFHH.mat"];
cAdjOutput = ["cp_adj_LH10.mat", "cp_adj_SFHH.mat"];

% Loop through each data set (face to face and co-presence)
for set=1:size(face2faceData,2)
    % Load in dataset and metadata set with indices for re-indexing based on 
    % number of nodes in network.
    f2fDataSet = load(face2faceData(set));
    cpDataSet = load(coPresenceData(set));
    sortedIndices = importdata(metaData(set));

    for i = 2:3 %Re-index column 2 and 3 (i and j) between 1 and # of nodes
        for j=1:size(f2fDataSet,1)
            % Face to face data
            reIndex = find(sortedIndices == f2fDataSet(j,i));
            reIndex = reIndex(1);
            f2fDataSet(j,i) = reIndex;
        end
        for j=1:size(cpDataSet,1) 
            % Co-Presence data
            reIndex = find(sortedIndices == cpDataSet(j,i));
            reIndex = reIndex(1);
            cpDataSet(j,i) = reIndex;
        end
    end

    % Loop through re-indexed columns to populate the agg. adjacency matrix
    % By adding a edge to the adjacency matrix for every instance in time
    % that i and j are in contact
    fAdjacencyMatrix = zeros(size(sortedIndices,1), size(sortedIndices,1));
    cAdjacencyMatrix = zeros(size(sortedIndices,1), size(sortedIndices,1));
    for i=1:size(f2fDataSet,1)
        fAdjacencyMatrix(f2fDataSet(i,2), f2fDataSet(i,3)) = fAdjacencyMatrix(f2fDataSet(i,2), f2fDataSet(i,3)) + 1;
        fAdjacencyMatrix(f2fDataSet(i,3), f2fDataSet(i,2)) = fAdjacencyMatrix(f2fDataSet(i,3), f2fDataSet(i,2)) + 1;
    end
    for i=1:size(cpDataSet,1)
        cAdjacencyMatrix(cpDataSet(i,2), cpDataSet(i,3)) = cAdjacencyMatrix(cpDataSet(i,2), cpDataSet(i,3)) + 1;
        cAdjacencyMatrix(cpDataSet(i,3), cpDataSet(i,2)) = cAdjacencyMatrix(cpDataSet(i,3), cpDataSet(i,2)) + 1;
    end
    
    % Save adjacency matrices to current directory
    data = fAdjacencyMatrix;
    save(fAdjOutput(set), 'data')
    data = cAdjacencyMatrix;
    save(cAdjOutput(set), 'data')
end

















% 
% sortedTest = sortrows(test,2);
% uniqueValues = 1;
% prev = 0;
% agg = zeros(size(test,1)*2,1);
% agg = sortedTest(1,2);
% agg(1:size(test,1),1) = sortedTest(:,2);
% agg(size(test,1)+1:size(test,1)*2,1) = sortedTest(:,3);
% agg = sort(agg);
% unique = zeros(81,1);
% 
% % Reindex according to indices given
% prev = agg(1,1);
% for i=1:size(test,1)*2
%     current = agg(i,1);
%     if(current ~= prev)
%         unique(uniqueValues,1) = prev;
%         prev = current;
%         uniqueValues = uniqueValues + 1;
%     end
% end
