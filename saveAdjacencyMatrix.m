%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clint Olsen
% ECEN 5322: Higher-Dimensional Datasets
% Final Project: Assignment 1 and 2
% Temporally Aggregated Adjacency Matrix Generation
% Creates and saves adjacency matrix from list of 
% datasets in face2faceData and coPresenceData vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - saveAdjacencyMatrix.m: Saves a face-to-face and co-presence
% adjacency matrix to the current directory
%
% - Inputs: face to face file to load (.dat), co-presence file (.dat)
% sorted index file from meta data (.mat), f2f adj. matrix file name
% to save to (.mat), and co presence file name to save to (.mat).
%
% - Outputs: Saved files in current directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveAdjacencyMatrix(f2fFile, cpFile, sortedIndicesFile, f2fSaveFileName, cpSaveFileName)
    f2fDataSet = load(f2fFile);
    cpDataSet = load(cpFile);
    sortedIndices = importdata(sortedIndicesFile);

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
    save(f2fSaveFileName, 'data')
    data = cAdjacencyMatrix;
    save(cpSaveFileName, 'data')

end