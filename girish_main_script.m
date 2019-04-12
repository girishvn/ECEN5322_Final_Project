% girish main script

f2f_A_files = ["f2f_adj_InVS13.mat", "f2f_adj_InVS15.mat", "f2f_adj_LH10.mat", "f2f_adj_LyonSchool.mat", "f2f_adj_SFHH.mat", "f2f_adj_Thiers13.mat"];
cp_A_files = ["cp_adj_InVS13.mat", "cp_adj_InVS15.mat", "cp_adj_LH10.mat", "cp_adj_LyonSchool.mat", "cp_adj_SFHH.mat", "cp_adj_Thiers13.mat"];

for i = 1:length(f2f_A_files) % iterate through data sets
    
    %Original, before sampling
    
    % load in f2f A
    f2fA =  load(f2f_A_files(i));
    f2fA = f2fA.data;
    
    % load in cp A
    cpA = load(cp_A_files(i));
    cpA = cpA.data;
    
    [f2fgraphSize, f2fnumEdges, f2fvolume, f2fdensity, f2fdegreeDistribution, f2fclusteringCoefficient, f2faverageDegree, f2favgCC] = computeStatistics(f2fA,false);
    [cpgraphSize, cpnumEdges, cpvolume, cpdensity, cpdegreeDistribution, cpclusteringCoefficient, cpaverageDegree, cpavgCC] = computeStatistics(cpA,false);

    figure()
    subplot(2,2,1)
    histogram(log(f2fclusteringCoefficient))
    set(gca, 'xscale','log')
    
    subplot(2,2,2)
    histogram(log(cpclusteringCoefficient))
    set(gca, 'xscale','log')
    
    subplot(2,2,3)
    %histogram(log(f2fdegreeDistribution))
    %set(gca, 'xscale','log')
    histogram(f2fdegreeDistribution)
    
    subplot(2,2,4)
    %histogram(log(cpdegreeDistribution))
    %set(gca, 'xscale','log')
    histogram(cpdegreeDistribution)
    
end







% TEEEEEEESSSSSSTTTTTTT
% computeStatistics(A);
% As1 = metropolisHastingsRW(A,1381*0.5);
% computeStatistics(As1);
% As2 = frontierSampling(A,1381*0.5);
% computeStatistics(As2);
% As3 = snowBallExpansion(A, 81*0.5);
% computeStatistics(As3);
% As4 = InducedGraphSampling(A, "weighted", 81*0.5);
% computeStatistics(As4);
% As5 = edgeSampling(A, "weighted", 1381*0.5);
% computeStatistics(As5);
% 
% %1381 vs 1156
% 
% subplot(2,5,1);
% image(A);
% title('Original Adjacency Matrix Lyon Hospital')
% 
% subplot(2,5,2);
% image(A);
% title('Original Adjacency Matrix Lyon Hospital')
% 
% subplot(2,5,3);
% image(A);
% title('Original Adjacency Matrix Lyon Hospital')
% 
% subplot(2,5,4);
% image(A);
% title('Original Adjacency Matrix Lyon Hospital')
% 
% subplot(2,5,5);
% image(A);
% title('Original Adjacency Matrix Lyon Hospital')
% 
% subplot(2,5,6);
% image(As4);
% title('Induced Sampling 50% Vertex Retension')
% 
% subplot(2,5,7);
% image(As5);
% title('Edge Sampling 50% Edge Retension')
% 
% subplot(2,5,8);
% image(As1);
% title('Metropolis Hastings RW 50% Edge Retension')
% 
% subplot(2,5,9);
% image(As2);
% title('Frontier Sampling 50% Edge Retension')
% 
% subplot(2,5,10);
% image(As3);
% title('Snow Ball Expansion 50% Vertex Retension')
