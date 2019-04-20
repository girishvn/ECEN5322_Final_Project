face2faceData = ["tij_InVS13.dat", "tij_InVS15.dat", "tij_LH10.dat", "tij_LyonSchool.dat", "tij_SFHH.dat", "tij_Thiers13.dat"];
coPresenceData = ["tij_pres_InVS13.dat", "tij_pres_InVS15.dat", "tij_pres_LH10.dat", "tij_pres_LyonSchool.dat", "tij_pres_SFHH.dat", "tij_pres_Thiers13.dat"];
metaData = ["sorted_indices_lnVS13.mat", "sorted_indices_lnVS15.mat", "sorted_indices_LH10.mat", "sorted_indices_LyonSchool.mat",   "sorted_indices_SFHH.mat", "sorted_indices_Thiers13.mat"];
fAdjOutput = ["f2f_adj_InVS13.mat", "f2f_adj_InVS15.mat", "f2f_adj_LH10.mat", "f2f_adj_LyonSchool.mat", "f2f_adj_SFHH.mat", "f2f_adj_Thiers13.mat"];
cAdjOutput = ["cp_adj_InVS13.mat", "cp_adj_InVS15.mat", "cp_adj_LH10.mat", "cp_adj_LyonSchool.mat", "cp_adj_SFHH.mat", "cp_adj_Thiers13.mat"];

%%
%Q1-2 Adjacency Matrix
for set=1:size(face2faceData,2)
    saveAdjacencyMatrix(face2faceData(set), coPresenceData(set), metaData(set), fAdjOutput(set), cAdjOutput(set));
end
%%
%Q3-4 Statistical Bar Graphs
bar_graphSize = [];
bar_numEdges = [];
bar_volume = [];
bar_density = [];
bar_degreeDistribution = [];
bar_clusteringCoefficient = [];
bar_averageDegree = [];
bar_avgCC = [];

adj_pairs = ["f2f_adj_InVS13.mat", "cp_adj_InVS13.mat", ...
                 "f2f_adj_InVS15.mat", "cp_adj_InVS15.mat", ...
                 "f2f_adj_LH10.mat", "cp_adj_LH10.mat", ...
                 "f2f_adj_LyonSchool.mat", "cp_adj_LyonSchool.mat", ...
                 "f2f_adj_SFHH.mat", "cp_adj_SFHH.mat", ...
                 "f2f_adj_Thiers13.mat", "cp_adj_Thiers13.mat"];
             
for adjMatrix=1:size(adj_pairs, 2)
    curr_adjMatrix = load(adj_pairs(adjMatrix));
    data_adjMatrix = curr_adjMatrix.data;
    [graphSize, numEdges, volume, density, degreeDistribution, clusteringCoefficient, averageDegree, avgCC] = computeStatistics(data_adjMatrix, false);
    bar_graphSize = [bar_graphSize, graphSize];
    bar_numEdges = [bar_numEdges, numEdges];
    bar_volume = [bar_volume, volume];
    bar_density = [bar_density, density];
    bar_degreeDistribution = [bar_degreeDistribution, degreeDistribution];
    bar_clusteringCoefficient = [bar_clusteringCoefficient, clusteringCoefficient];
    bar_averageDegree = [bar_averageDegree, averageDegree];
    bar_avgCC = [bar_avgCC, avgCC];   
end

adj_pairs_names = ["f2f-adj-InVS13", "cp-adj-InVS13", ...
                 "f2f-adj-InVS15", "cp-adj-InVS15", ...
                 "f2f-adj-LH10", "cp-adj-LH10", ...
                 "f2f-adj-LyonSchool", "cp-adj-LyonSchool", ...
                 "f2f-adj_SFHH", "cp-adj-SFHH", ...
                 "f2f-adj_Thiers13", "cp-adj-Thiers13"];
             

             
c = categorical(adj_pairs_names);
c = reordercats(c, {'f2f-adj-InVS13', 'cp-adj-InVS13', ...
                 'f2f-adj-InVS15', 'cp-adj-InVS15', ...
                 'f2f-adj-LH10', 'cp-adj-LH10', ...
                 'f2f-adj-LyonSchool', 'cp-adj-LyonSchoomat', ...
                 'f2f-adj_SFHH', 'cp-adj-SFHH', ...
                 'f2f-adj_Thiers13', 'cp-adj-Thiers13'});

subplot(2,3,1)
bar(c, bar_graphSize)
xlabel('Adjacency Data Sets') 
ylabel('Graph Size') 
title(sprintf("Face to Face & Co-Presence\n Graph Sizes"))

subplot(2,3,2)
bar(c, bar_numEdges)
xlabel('Adjacency Data Sets') 
ylabel('Number of Edges') 
title(sprintf("Face to Face & Co-Presence\n Number of Edges"))

subplot(2,3,3)
bar(c, bar_volume)
xlabel('Adjacency Data Sets') 
ylabel('Volume') 
title(sprintf("Face to Face & Co-Presence\n Volumes"))

subplot(2,3,4)
bar(c, bar_density)
xlabel('Adjacency Data Sets') 
ylabel('Density') 
title(sprintf("Face to Face & Co-Presence\n Densities"))

subplot(2,3,5)
bar(c, bar_averageDegree)
xlabel('Adjacency Data Sets') 
ylabel('Average Degree') 
title(sprintf("Face to Face & Co-Presence\n Average Degrees"))

subplot(2,3,6)
bar(c, bar_avgCC)
xlabel('Adjacency Data Sets') 
ylabel('Average Clustering Coefficient') 
title(sprintf("Face to Face & Co-Presence\n Average Clustering Coefficient"))


%%
%Q5-6 Statistical Line graphs

%Initialize Structures to Parse and Fill Data for Plots
adjData = ["f2f_adj_InVS13.mat", "f2f_adj_InVS15.mat", "f2f_adj_LH10.mat", "f2f_adj_LyonSchool.mat", "f2f_adj_SFHH.mat", "f2f_adj_Thiers13.mat", ...
           "cp_adj_InVS13.mat", "cp_adj_InVS15.mat", "cp_adj_LH10.mat", "cp_adj_LyonSchool.mat", "cp_adj_SFHH.mat", "cp_adj_Thiers13.mat"];
adjData_names = ["f2f-adj-InVS13", "f2f-adj-InVS15", "f2f-adj-LH10", "f2f-adj-LyonSchool", "f2f-adj-SFHH", "f2f-adj-Thiers13", ...
           "cp-adj-InVS13", "cp-adj-InVS15", "cp-adj-LH10", "cp-adj-LyonSchool", "cp-adj-SFHH", "cp-adj-Thiers13"];
frequencies = [0.5, 0.6, 0.7, 0.8, 0.9];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%ALGORITHM 1%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Go Through Every Data Set and Populate Master Stat Bank With It's
%Statistics After Being Reduced By Algorithm 1
masterStatBank = [];
for adjMatrix=1:size(adjData, 2)
    %Get all combination of Stats for a Matrix
    curr_adjMatrix = load(adjData(adjMatrix));
    data_adjMatrix = curr_adjMatrix.data;
    dataStatBank = zeros(6, 5);
    i = 1;
    %Get total num of vertices
    [totalVertices, t2, t3, t4, t5, t6, t7, t8] = computeStatistics(data_adjMatrix, false);
    
    %Find Statistics For Every Kind of Frequencies
    for f=1:size(frequencies, 2)
        reduced = InducedGraphSampling(data_adjMatrix, "uniform", totalVertices*frequencies(f));
        [s1, s2, s3, s4, s5, s6, s7, s8] = computeStatistics(reduced, true);
        dataStatBank(1, i) = s1;
        dataStatBank(2, i) = s2;
        dataStatBank(3, i) = s3;
        dataStatBank(4, i) = s4;
        dataStatBank(5, i) = s7;
        dataStatBank(6, i) = s8;
        i = i + 1;
    end
    %Append the Master Stat Bank With All Freq-Stats For the Dataset Ran Through Algorithm 1
    masterStatBank = [masterStatBank; dataStatBank];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate Line Graph     %
%For Every Stat vs. Freq %
%for Every Data Set      %
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
%Create Plot for Graph Size
subplot(2,3,1)
for i = 0:11
   plot(frequencies, masterStatBank(1+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
lgd = legend;
lgd.NumColumns = 4;
xlabel('Frequencies') 
ylabel('Graph Size') 
title(sprintf("Algorithm 1 - Induced Graph Sampling(Uniform) \n Graph Sizes")) 

%Create Plot for Num Edges
subplot(2,3,2)
for i = 0:11
   plot(frequencies, masterStatBank(2+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Number of Edges') 
title(sprintf("Algorithm 1 - Induced Graph Sampling(Uniform) \n Number of Edges"))  

%Create Plot for Volumes
subplot(2,3,3)
for i = 0:11
   plot(frequencies, masterStatBank(3+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Volume') 
title(sprintf("Algorithm 1 - Induced Graph Sampling(Uniform) \n Volume")) 

%Create Plot for Densities
subplot(2,3,4)
for i = 0:11
   plot(frequencies, masterStatBank(4+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Densities') 
title(sprintf("Algorithm 1 - Induced Graph Sampling(Uniform) \n Densities")) 

%Create Plot for Average Degrees
subplot(2,3,5)
for i = 0:11
   plot(frequencies, masterStatBank(5+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Average Degree') 
title(sprintf("Algorithm 1 - Induced Graph Sampling(Uniform) \n Average Degrees")) 

%Create Plot for Average Coefficients
subplot(2,3,6)
for i = 0:11
   plot(frequencies, masterStatBank(6+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Average Coefficients') 
title(sprintf("Algorithm 1 - Induced Graph Sampling(Uniform) \n Average Coefficients")) 

disp("Algorithm 1 Done")

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%ALGORITHM 1.2%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Go Through Every Data Set and Populate Master Stat Bank With It's
%Statistics After Being Reduced By Algorithm 1.2
masterStatBank = [];
for adjMatrix=1:size(adjData, 2)
    %Get all combination of Stats for a Matrix
    curr_adjMatrix = load(adjData(adjMatrix));
    data_adjMatrix = curr_adjMatrix.data;
    dataStatBank = zeros(6, 5);
    i = 1;
  
    %Get total num of vertices
    [totalVertices, t2, t3, t4, t5, t6, t7, t8] = computeStatistics(data_adjMatrix, false);

    %Find Statistics For Every Kind of Frequencies
    for f=1:size(frequencies, 2)
        reduced = InducedGraphSampling(data_adjMatrix, "weighted", totalVertices*frequencies(f));
        [s1, s2, s3, s4, s5, s6, s7, s8] = computeStatistics(reduced, true);
        dataStatBank(1, i) = s1;
        dataStatBank(2, i) = s2;
        dataStatBank(3, i) = s3;
        dataStatBank(4, i) = s4;
        dataStatBank(5, i) = s7;
        dataStatBank(6, i) = s8;
        i = i + 1;
    end
    %Append the Master Stat Bank With All Freq-Stats For the Dataset Ran Through Algorithm 1
    masterStatBank = [masterStatBank; dataStatBank];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate Line Graph     %
%for Every Stat vs. Freq %
%for Every Data Set      %
%for Algorithm 1.2       %
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
%Create Plot for Graph Size
subplot(2,3,1)
for i = 0:11
   plot(frequencies, masterStatBank(1+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
lgd = legend;
lgd.NumColumns = 4;
xlabel('Frequencies') 
ylabel('Graph Size') 
title(sprintf("Algorithm 1 - Induced Graph Sampling(Weighted) \n Graph Sizes")) 

%Create Plot for Num Edges
subplot(2,3,2)
for i = 0:11
   plot(frequencies, masterStatBank(2+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Number of Edges') 
title(sprintf("Algorithm 1 - Induced Graph Sampling(Weighted) \n Number of Edges"))  

%Create Plot for Volumes
subplot(2,3,3)
for i = 0:11
   plot(frequencies, masterStatBank(3+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Volume') 
title(sprintf("Algorithm 1 - Induced Graph Sampling(Weighted) \n Volume")) 

%Create Plot for Densities
subplot(2,3,4)
for i = 0:11
   plot(frequencies, masterStatBank(4+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Densities') 
title(sprintf("Algorithm 1 - Induced Graph Sampling(Weighted) \n Densities")) 

%Create Plot for Average Degrees
subplot(2,3,5)
for i = 0:11
   plot(frequencies, masterStatBank(5+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Average Degree') 
title(sprintf("Algorithm 1 - Induced Graph Sampling(Weighted) \n Average Degrees")) 

%Create Plot for Average Coefficients
subplot(2,3,6)
for i = 0:11
   plot(frequencies, masterStatBank(6+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Average Coefficients') 
title(sprintf("Algorithm 1 - Induced Graph Sampling(Weighted) \n Average Coefficients")) 

disp("Algorithm 1.2 Done")

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%ALGORITHM 2%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Go Through Every Data Set and Populate Master Stat Bank With It's
%Statistics After Being Reduced By Algorithm 1.2
masterStatBank = [];
for adjMatrix=1:size(adjData, 2)
    %Get all combination of Stats for a Matrix
    curr_adjMatrix = load(adjData(adjMatrix));
    data_adjMatrix = curr_adjMatrix.data;
    dataStatBank = zeros(6, 5);
    i = 1;
    
    %Get total num of edges
    [t1, totalEdges, t3, t4, t5, t6, t7, t8] = computeStatistics(data_adjMatrix, false);
    %Find Statistics For Every Kind of Frequencies
    for f=1:size(frequencies, 2)
        reduced = edgeSampling(data_adjMatrix, "uniform", totalEdges*frequencies(f));
        [s1, s2, s3, s4, s5, s6, s7, s8] = computeStatistics(reduced, true);
        dataStatBank(1, i) = s1;
        dataStatBank(2, i) = s2;
        dataStatBank(3, i) = s3;
        dataStatBank(4, i) = s4;
        dataStatBank(5, i) = s7;
        dataStatBank(6, i) = s8;
        i = i + 1;
    end
    %Append the Master Stat Bank With All Freq-Stats For the Dataset Ran Through Algorithm 1
    masterStatBank = [masterStatBank; dataStatBank];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate Line Graph     %
%for Every Stat vs. Freq %
%for Every Data Set      %
%for Algorithm 2         %
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
%Create Plot for Graph Size
subplot(2,3,1)
for i = 0:11
   plot(frequencies, masterStatBank(1+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
lgd = legend;
lgd.NumColumns = 4;
xlabel('Frequencies') 
ylabel('Graph Size') 
title(sprintf("Algorithm 2 - Edge Sampling(Uniform) \n Graph Sizes")) 

%Create Plot for Num Edges
subplot(2,3,2)
for i = 0:11
   plot(frequencies, masterStatBank(2+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Number of Edges') 
title(sprintf("Algorithm 2 - Edge Sampling(Uniform) \n Number of Edges"))  

%Create Plot for Volumes
subplot(2,3,3)
for i = 0:11
   plot(frequencies, masterStatBank(3+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Volume') 
title(sprintf("Algorithm 2 - Edge Sampling(Uniform) \n Volume")) 

%Create Plot for Densities
subplot(2,3,4)
for i = 0:11
   plot(frequencies, masterStatBank(4+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Densities') 
title(sprintf("Algorithm 2 - Edge Sampling(Uniform) \n Densities")) 

%Create Plot for Average Degrees
subplot(2,3,5)
for i = 0:11
   plot(frequencies, masterStatBank(5+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Average Degree') 
title(sprintf("Algorithm 2 - Edge Sampling(Uniform) \n Average Degrees")) 

%Create Plot for Average Coefficients
subplot(2,3,6)
for i = 0:11
   plot(frequencies, masterStatBank(6+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Average Coefficients') 
title(sprintf("Algorithm 2 - Edge Sampling(Uniform) \n Average Coefficients")) 

disp("Algorithm 2 Done")

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%ALGORITHM 2.2%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Go Through Every Data Set and Populate Master Stat Bank With It's
%Statistics After Being Reduced By Algorithm 1.2
masterStatBank = [];
for adjMatrix=1:size(adjData, 2)
    %Get all combination of Stats for a Matrix
    curr_adjMatrix = load(adjData(adjMatrix));
    data_adjMatrix = curr_adjMatrix.data;
    dataStatBank = zeros(6, 5);
    i = 1;
    
    %Get total num of edges
    [t1, totalEdges, t3, t4, t5, t6, t7, t8] = computeStatistics(data_adjMatrix, false);
    %Find Statistics For Every Kind of Frequencies
    for f=1:size(frequencies, 2)
        reduced = edgeSampling(data_adjMatrix, "weighted", totalEdges*frequencies(f));
        [s1, s2, s3, s4, s5, s6, s7, s8] = computeStatistics(reduced, true);
        dataStatBank(1, i) = s1;
        dataStatBank(2, i) = s2;
        dataStatBank(3, i) = s3;
        dataStatBank(4, i) = s4;
        dataStatBank(5, i) = s7;
        dataStatBank(6, i) = s8;
        i = i + 1;
    end
    %Append the Master Stat Bank With All Freq-Stats For the Dataset Ran Through Algorithm 1
    masterStatBank = [masterStatBank; dataStatBank];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate Line Graph     %
%for Every Stat vs. Freq %
%for Every Data Set      %
%for Algorithm 2.2       %
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
%Create Plot for Graph Size
subplot(2,3,1)
for i = 0:11
   plot(frequencies, masterStatBank(1+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
lgd = legend;
lgd.NumColumns = 4;
xlabel('Frequencies') 
ylabel('Graph Size') 
title(sprintf("Algorithm 2 - Edge Sampling(weighted) \n Graph Sizes")) 

%Create Plot for Num Edges
subplot(2,3,2)
for i = 0:11
   plot(frequencies, masterStatBank(2+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Number of Edges') 
title(sprintf("Algorithm 2 - Edge Sampling(weighted) \n Number of Edges"))  

%Create Plot for Volumes
subplot(2,3,3)
for i = 0:11
   plot(frequencies, masterStatBank(3+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Volume') 
title(sprintf("Algorithm 2 - Edge Sampling(weighted) \n Volume")) 

%Create Plot for Densities
subplot(2,3,4)
for i = 0:11
   plot(frequencies, masterStatBank(4+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Densities') 
title(sprintf("Algorithm 2 - Edge Sampling(weighted) \n Densities")) 

%Create Plot for Average Degrees
subplot(2,3,5)
for i = 0:11
   plot(frequencies, masterStatBank(5+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Average Degree') 
title(sprintf("Algorithm 2 - Edge Sampling(weighted) \n Average Degrees")) 

%Create Plot for Average Coefficients
subplot(2,3,6)
for i = 0:11
   plot(frequencies, masterStatBank(6+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Average Coefficients') 
title(sprintf("Algorithm 2 - Edge Sampling(weighted) \n Average Coefficients")) 

disp("Algorithm 2.2 Done")

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%ALGORITHM 3%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Go Through Every Data Set and Populate Master Stat Bank With It's
%Statistics After Being Reduced By Algorithm 3
masterStatBank = [];
for adjMatrix=1:size(adjData, 2)
    %Get all combination of Stats for a Matrix
    curr_adjMatrix = load(adjData(adjMatrix));
    data_adjMatrix = curr_adjMatrix.data;
    dataStatBank = zeros(6, 5);
    i = 1;
    
    %Get total num of edges
    [t1, totalEdges, t3, t4, t5, t6, t7, t8] = computeStatistics(data_adjMatrix, false);
    %Find Statistics For Every Kind of Frequencies
    for f=1:size(frequencies, 2)
        reduced = metropolisHastingsRW(data_adjMatrix, totalEdges*frequencies(f));
        [s1, s2, s3, s4, s5, s6, s7, s8] = computeStatistics(reduced, true);
        dataStatBank(1, i) = s1;
        dataStatBank(2, i) = s2;
        dataStatBank(3, i) = s3;
        dataStatBank(4, i) = s4;
        dataStatBank(5, i) = s7;
        dataStatBank(6, i) = s8;
        i = i + 1;
    end
    %Append the Master Stat Bank With All Freq-Stats For the Dataset Ran Through Algorithm 1
    masterStatBank = [masterStatBank; dataStatBank];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate Line Graph     %
%for Every Stat vs. Freq %
%for Every Data Set      %
%for Algorithm 3         %
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
%Create Plot for Graph Size
subplot(2,3,1)
for i = 0:11
   plot(frequencies, masterStatBank(1+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
lgd = legend;
lgd.NumColumns = 4;
xlabel('Frequencies') 
ylabel('Graph Size') 
title(sprintf("Algorithm 3 - Metropolis-Hastings Random Walk \n Graph Sizes")) 

%Create Plot for Num Edges
subplot(2,3,2)
for i = 0:11
   plot(frequencies, masterStatBank(2+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Number of Edges') 
title(sprintf("Algorithm 3 - Metropolis-Hastings Random Walk \n Number of Edges"))  

%Create Plot for Volumes
subplot(2,3,3)
for i = 0:11
   plot(frequencies, masterStatBank(3+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Volume') 
title(sprintf("Algorithm 3 - Metropolis-Hastings Random Walk \n Volume")) 

%Create Plot for Densities
subplot(2,3,4)
for i = 0:11
   plot(frequencies, masterStatBank(4+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Densities') 
title(sprintf("Algorithm 3 - Metropolis-Hastings Random Walk \n Densities")) 

%Create Plot for Average Degrees
subplot(2,3,5)
for i = 0:11
   plot(frequencies, masterStatBank(5+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Average Degree') 
title(sprintf("Algorithm 3 - Metropolis-Hastings Random Walk \n Average Degrees")) 

%Create Plot for Average Coefficients
subplot(2,3,6)
for i = 0:11
   plot(frequencies, masterStatBank(6+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Average Coefficients') 
title(sprintf("Algorithm 3 - Metropolis-Hastings Random Walk \n Average Coefficients")) 

disp("Algorithm 3 Done")

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%ALGORITHM 4%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Go Through Every Data Set and Populate Master Stat Bank With It's
%Statistics After Being Reduced By Algorithm 4
masterStatBank = [];
for adjMatrix=1:size(adjData, 2)
    %Get all combination of Stats for a Matrix
    curr_adjMatrix = load(adjData(adjMatrix));
    data_adjMatrix = curr_adjMatrix.data;
    dataStatBank = zeros(6, 5);
    i = 1;
    
    %Get total num of edges
    [t1, totalEdges, t3, t4, t5, t6, t7, t8] = computeStatistics(data_adjMatrix, false);
    %Find Statistics For Every Kind of Frequencies
    for f=1:size(frequencies, 2)
        reduced = frontierSampling(data_adjMatrix, totalEdges*frequencies(f));
        [s1, s2, s3, s4, s5, s6, s7, s8] = computeStatistics(reduced, true);
        dataStatBank(1, i) = s1;
        dataStatBank(2, i) = s2;
        dataStatBank(3, i) = s3;
        dataStatBank(4, i) = s4;
        dataStatBank(5, i) = s7;
        dataStatBank(6, i) = s8;
        i = i + 1;
    end
    %Append the Master Stat Bank With All Freq-Stats For the Dataset Ran Through Algorithm 1
    masterStatBank = [masterStatBank; dataStatBank];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate Line Graph     %
%for Every Stat vs. Freq %
%for Every Data Set      %
%for Algorithm 4         %
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
%Create Plot for Graph Size
subplot(2,3,1)
for i = 0:11
   plot(frequencies, masterStatBank(1+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
lgd = legend;
lgd.NumColumns = 4;
xlabel('Frequencies') 
ylabel('Graph Size') 
title(sprintf("Algorithm 4 - Frontier Sampling \n Graph Sizes")) 

%Create Plot for Num Edges
subplot(2,3,2)
for i = 0:11
   plot(frequencies, masterStatBank(2+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Number of Edges') 
title(sprintf("Algorithm 4 - Frontier Sampling \n Number of Edges"))  

%Create Plot for Volumes
subplot(2,3,3)
for i = 0:11
   plot(frequencies, masterStatBank(3+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Volume') 
title(sprintf("Algorithm 4 - Frontier Sampling \n Volume")) 

%Create Plot for Densities
subplot(2,3,4)
for i = 0:11
   plot(frequencies, masterStatBank(4+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Densities') 
title(sprintf("Algorithm 4 - Frontier Sampling \n Densities")) 

%Create Plot for Average Degrees
subplot(2,3,5)
for i = 0:11
   plot(frequencies, masterStatBank(5+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Average Degree') 
title(sprintf("Algorithm 4 - Frontier Sampling \n Average Degrees")) 

%Create Plot for Average Coefficients
subplot(2,3,6)
for i = 0:11
   plot(frequencies, masterStatBank(6+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Average Coefficients') 
title(sprintf("Algorithm 4 - Frontier Sampling \n Average Coefficients")) 

disp("Algorithm 4 Done")

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%ALGORITHM 4%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Go Through Every Data Set and Populate Master Stat Bank With It's
%Statistics After Being Reduced By Algorithm 4
masterStatBank = [];
for adjMatrix=1:size(adjData, 2)
    %Get all combination of Stats for a Matrix
    curr_adjMatrix = load(adjData(adjMatrix));
    data_adjMatrix = curr_adjMatrix.data;
    dataStatBank = zeros(6, 5);
    i = 1;
    
    %Get total num of edges
    [t1, totalEdges, t3, t4, t5, t6, t7, t8] = computeStatistics(data_adjMatrix, false);
    %Find Statistics For Every Kind of Frequencies
    for f=1:size(frequencies, 2)
        reduced = frontierSampling(data_adjMatrix, totalEdges*frequencies(f));
        [s1, s2, s3, s4, s5, s6, s7, s8] = computeStatistics(reduced, true);
        dataStatBank(1, i) = s1;
        dataStatBank(2, i) = s2;
        dataStatBank(3, i) = s3;
        dataStatBank(4, i) = s4;
        dataStatBank(5, i) = s7;
        dataStatBank(6, i) = s8;
        i = i + 1;
    end
    %Append the Master Stat Bank With All Freq-Stats For the Dataset Ran Through Algorithm 1
    masterStatBank = [masterStatBank; dataStatBank];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate Line Graph     %
%for Every Stat vs. Freq %
%for Every Data Set      %
%for Algorithm 4         %
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
%Create Plot for Graph Size
subplot(2,3,1)
for i = 0:11
   plot(frequencies, masterStatBank(1+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
lgd = legend;
lgd.NumColumns = 4;
xlabel('Frequencies') 
ylabel('Graph Size') 
title(sprintf("Algorithm 4 - Frontier Sampling \n Graph Sizes")) 

%Create Plot for Num Edges
subplot(2,3,2)
for i = 0:11
   plot(frequencies, masterStatBank(2+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Number of Edges') 
title(sprintf("Algorithm 4 - Frontier Sampling \n Number of Edges"))  

%Create Plot for Volumes
subplot(2,3,3)
for i = 0:11
   plot(frequencies, masterStatBank(3+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Volume') 
title(sprintf("Algorithm 4 - Frontier Sampling \n Volume")) 

%Create Plot for Densities
subplot(2,3,4)
for i = 0:11
   plot(frequencies, masterStatBank(4+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Densities') 
title(sprintf("Algorithm 4 - Frontier Sampling \n Densities")) 

%Create Plot for Average Degrees
subplot(2,3,5)
for i = 0:11
   plot(frequencies, masterStatBank(5+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Average Degree') 
title(sprintf("Algorithm 4 - Frontier Sampling \n Average Degrees")) 

%Create Plot for Average Coefficients
subplot(2,3,6)
for i = 0:11
   plot(frequencies, masterStatBank(6+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Average Coefficients') 
title(sprintf("Algorithm 4 - Frontier Sampling \n Average Coefficients")) 

disp("Algorithm 4 Done")

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%ALGORITHM 5%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Go Through Every Data Set and Populate Master Stat Bank With It's
%Statistics After Being Reduced By Algorithm 5
masterStatBank = [];
for adjMatrix=1:size(adjData, 2)
    %Get all combination of Stats for a Matrix
    curr_adjMatrix = load(adjData(adjMatrix));
    data_adjMatrix = curr_adjMatrix.data;
    dataStatBank = zeros(6, 5);
    i = 1;
    
    %Get total num of edges
    [totalVertices, t2, t3, t4, t5, t6, t7, t8] = computeStatistics(data_adjMatrix, false);
    %Find Statistics For Every Kind of Frequencies
    for f=1:size(frequencies, 2)
        reduced = snowBallExpansion(data_adjMatrix, totalVertices*frequencies(f));
        [s1, s2, s3, s4, s5, s6, s7, s8] = computeStatistics(reduced, true);
        dataStatBank(1, i) = s1;
        dataStatBank(2, i) = s2;
        dataStatBank(3, i) = s3;
        dataStatBank(4, i) = s4;
        dataStatBank(5, i) = s7;
        dataStatBank(6, i) = s8;
        i = i + 1;
    end
    %Append the Master Stat Bank With All Freq-Stats For the Dataset Ran Through Algorithm 5
    masterStatBank = [masterStatBank; dataStatBank];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate Line Graph     %
%for Every Stat vs. Freq %
%for Every Data Set      %
%for Algorithm 5         %
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
%Create Plot for Graph Size
subplot(2,3,1)
for i = 0:11
   plot(frequencies, masterStatBank(1+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
lgd = legend;
lgd.NumColumns = 4;
xlabel('Frequencies') 
ylabel('Graph Size') 
title(sprintf("Algorithm 5 - Snowball Expansion \n Graph Sizes")) 

%Create Plot for Num Edges
subplot(2,3,2)
for i = 0:11
   plot(frequencies, masterStatBank(2+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Number of Edges') 
title(sprintf("Algorithm 5 - Snowball Expansion \n Number of Edges"))  

%Create Plot for Volumes
subplot(2,3,3)
for i = 0:11
   plot(frequencies, masterStatBank(3+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Volume') 
title(sprintf("Algorithm 5 - Snowball Expansion \n Volume")) 

%Create Plot for Densities
subplot(2,3,4)
for i = 0:11
   plot(frequencies, masterStatBank(4+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Densities') 
title(sprintf("Algorithm 5 - Snowball Expansion \n Densities")) 

%Create Plot for Average Degrees
subplot(2,3,5)
for i = 0:11
   plot(frequencies, masterStatBank(5+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Average Degree') 
title(sprintf("Algorithm 5 - Snowball Expansion \n Average Degrees")) 

%Create Plot for Average Coefficients
subplot(2,3,6)
for i = 0:11
   plot(frequencies, masterStatBank(6+(i*6), :), 'DisplayName', adjData_names(i+1), 'LineWidth', 2)
   hold on
end
hold off
xlabel('Frequencies') 
ylabel('Average Coefficients') 
title(sprintf("Algorithm 5 - Snowball Expansion \n Average Coefficients")) 

disp("Algorithm 5 Done")