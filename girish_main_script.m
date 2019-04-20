%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Girish Narayanswamy
% ECEN 5322: Higher-Dimensional Datasets
% Final Project: Assignment 5 and 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the plotting script for the histograms
% for degree density and clustering coefficient 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f2f_A_files = ["f2f_adj_InVS13.mat", "f2f_adj_InVS15.mat", "f2f_adj_LH10.mat", "f2f_adj_LyonSchool.mat", "f2f_adj_SFHH.mat", "f2f_adj_Thiers13.mat"];
cp_A_files = ["cp_adj_InVS13.mat", "cp_adj_InVS15.mat", "cp_adj_LH10.mat", "cp_adj_LyonSchool.mat", "cp_adj_SFHH.mat", "cp_adj_Thiers13.mat"];
dataSets = ["InVS13", "InVS15", "LH10", "LyonSchool", "SFHH", "Thiers13"];
f = [0.5, 0.6, 0.7, 0.8, 0.9];

for i = 1:length(f2f_A_files) % iterate through data sets
    
    %Original, before sampling
    
    % load in f2f A
    f2fA =  load(f2f_A_files(i));
    f2fA = f2fA.data;
    
    % load in cp A
    cpA = load(cp_A_files(i));
    cpA = cpA.data;
    
    [f2fSize, f2fnumEdges, f2fvol, f2fdensity, f2fDD, f2fCC, f2favgDeg, f2favgCC] = computeStatistics(f2fA,false);
    [cpSize, cpnumEdges, cpvol, cpdensity, cpDD, cpCC, cpavgDeg, cpavgCC] = computeStatistics(cpA,false);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Figure titles and labels and stuff
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(30 + i) 
    sgtitle(strcat("Clustering Coeff. and Degree Distrib. for ", dataSets(i)))
    subplot(2,2,1)
    histogram(log(f2fCC))
    set(gca, 'xscale','log')
    title("Clustering Coeff. of F2F A Matrix")
    xlabel("log(Clustering Coeff)")
    ylabel("log(Occurrence)")
    
    subplot(2,2,2)
    histogram(log(cpCC))
    set(gca, 'xscale','log')
    title("Clustering Coeff. of CP A Matrix")
    xlabel("log(Clustering Coeff)")
    ylabel("log(Occurrence)")
    
    subplot(2,2,3)
    histogram(log(f2fDD))
    set(gca, 'xscale','log')
    title("Degree Distrib. of F2F A Matrix")
    xlabel("log(Degree)")
    ylabel("log(Occurrence)")

    subplot(2,2,4)
    histogram(log(cpDD))
    set(gca, 'xscale','log')
    title("Degree Distrib. of CP A Matrix")
    xlabel("log(Degree)")
    ylabel("log(Occurrence)")

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up figures
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(1 + i*5 - 5)
    set(gcf,'Position',[10 10 1200 1000])
    sgtitle(strcat("Induced Sampled Clustering Coeff. and Degree Distrib. for ", dataSets(i)))
    
    figure(2 + i*5 - 5)
    set(gcf,'Position',[10 10 1200 1000])
    sgtitle(strcat("Edge Sampled Sampled Clustering Coeff. and Degree Distrib. for ", dataSets(i)))
 
    figure(3 + i*5 - 5)
    set(gcf,'Position',[10 10 1200 1000])
    sgtitle(strcat("Metropolis Hastings RW Sampled Sampled Clustering Coeff. and Degree Distrib. for ", dataSets(i)))
   
    figure(4 + i*5 - 5)
    set(gcf,'Position',[10 10 1200 1000])
    sgtitle(strcat("Frontier Sampled Sampled Clustering Coeff. and Degree Distrib. for ", dataSets(i)))

    figure(5 + i*5 - 5)
    set(gcf,'Position',[10 10 1200 1000])
    sgtitle(strcat("Snow Ball Expansion Sampled Clustering Coeff. and Degree Distrib. for ", dataSets(i)))
    
    % Induced 
    for j = 1:length(f)
       fj = f(j);
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Run Sampling Algorithms on Adjacency Mtxs
       % Get Statistics on Sampled Adjacency Mtxs
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       % Sampled f2f graph info
       A1f2f = InducedGraphSampling(f2fA,"weighted",f2fSize*fj);
       [f2fSize1, f2fnumEdges1, f2fvol1, f2fdensity1, f2fDD1, f2fCC1, f2favgDeg1, f2favgCC1] = computeStatistics(A1f2f,true);
       
       A2f2f = edgeSampling(f2fA, "weighted", f2fnumEdges*fj);
       [f2fSize2, f2fnumEdges2, f2fvol2, f2fdensity2, f2fDD2, f2fCC2, f2favgDeg2, f2favgCC2] = computeStatistics(A2f2f,true);
       
       A3f2f = metropolisHastingsRW(f2fA,f2fnumEdges*fj);
       [f2fSize3, f2fnumEdges3, f2fvol3, f2fdensity3, f2fDD3, f2fCC3, f2favgDeg3, f2favgCC3] = computeStatistics(A3f2f,true);
       
       A4f2f = frontierSampling(f2fA,f2fnumEdges*fj);
       [f2fSize4, f2fnumEdges4, f2fvol4, f2fdensity4, f2fDD4, f2fCC4, f2favgDeg4, f2favgCC4] = computeStatistics(A4f2f,true);
       
       A5f2f = snowBallExpansion(f2fA,f2fSize*fj);
       [f2fSize5, f2fnumEdges5, f2fvol5, f2fdensity5, f2fDD5, f2fCC5, f2favgDeg5, f2favgCC5] = computeStatistics(A5f2f,true);
       
       % Sampled cp graph info
       A1cp = InducedGraphSampling(cpA,"weighted",cpSize*fj);
       [cpSize1, cpnumEdges1, cpvol1, cpdensity1, cpDD1, cpCC1, cpavgDeg1, cpavgCC1] = computeStatistics(A1cp,true);
       
       A2cp = edgeSampling(cpA, "weighted", cpnumEdges*fj);
       [cpSize2, cpnumEdges2, cpvol2, cpdensity2, cpDD2, cpCC2, cpavgDeg2, cpavgCC2] = computeStatistics(A2cp,true);
       
       A3cp = metropolisHastingsRW(cpA,cpnumEdges*fj);
       [cpSize3, cpnumEdges3, cpvol3, cpdensity3, cpDD3, cpCC3, cpavgDeg3, cpavgCC3] = computeStatistics(A3cp,true);
       
       A4cp = frontierSampling(cpA,cpnumEdges*fj);
       [cpSize4, cpnumEdges4, cpvol4, cpdensity4, cpDD4, cpCC4, cpavgDeg4, cpavgCC4] = computeStatistics(A4cp,true);
       
       A5cp = snowBallExpansion(cpA,cpSize*fj);
       [cpSize5, cpnumEdges5, cpvol5, cpdensity5, cpDD5, cpCC5, cpavgDeg5, cpavgCC5] = computeStatistics(A5cp,true);
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Figure formatting for DD and CC Plots
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       figure(1 + i*5 - 5) % algo 1 induced graph
       subplot(5,4,1 + j*4 - 4)
       histogram(log(f2fCC1))
       set(gca, 'xscale','log')
       title(strcat("CC of F2F f=",sprintf('%1.1f',fj)))
       xlabel("log(CC)")
       ylabel("log(Occur.)")
       
       subplot(5,4,2 + j*4 - 4)
       histogram(log(cpCC1))
       set(gca, 'xscale','log')
       title(strcat("CC of CP f=",sprintf('%1.1f',fj)))
       xlabel("log(CC)")
       ylabel("log(Occur.)")
       
       subplot(5,4,3 + j*4 - 4)
       histogram(log(f2fDD1))
       set(gca, 'xscale','log')
       title(strcat("DD of F2F f=",sprintf('%1.1f',fj)))
       xlabel("log(DD)")
       ylabel("log(Occur.)")
       
       subplot(5,4,4 + j*4 - 4)
       histogram(log(cpDD1))
       set(gca, 'xscale','log')
       title(strcat("DD of CP f=",sprintf('%1.1f',fj)))
       xlabel("log(DD)")
       ylabel("log(Occur.)")
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       figure(2 + i*5 - 5) % algo 2 edge sample graph
       subplot(5,4,1 + j*4 - 4)
       histogram(log(f2fCC2))
       set(gca, 'xscale','log')
       title(strcat("CC of F2F f=",sprintf('%1.1f',fj)))
       xlabel("log(CC)")
       ylabel("log(Occur.)")
       
       subplot(5,4,2 + j*4 - 4)
       histogram(log(cpCC2))
       set(gca, 'xscale','log')
       title(strcat("CC of CP f=",sprintf('%1.1f',fj)))
       xlabel("log(CC)")
       ylabel("log(Occur.)")
       
       subplot(5,4,3 + j*4 - 4)
       histogram(log(f2fDD2))
       set(gca, 'xscale','log')
       title(strcat("DD of F2F f=",sprintf('%1.1f',fj)))
       xlabel("log(DD)")
       ylabel("log(Occur.)")
       
       subplot(5,4,4 + j*4 - 4)
       histogram(log(cpDD2))
       set(gca, 'xscale','log')
       title(strcat("DD of CP f=",sprintf('%1.1f',fj)))
       xlabel("log(DD)")
       ylabel("log(Occur.)")
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       figure(3 + i*5 - 5) % algo 3 metropolis hastings rw graph
       subplot(5,4,1 + j*4 - 4)
       histogram(log(f2fCC3))
       set(gca, 'xscale','log')
       title(strcat("CC of F2F f=",sprintf('%1.1f',fj)))
       xlabel("log(CC)")
       ylabel("log(Occur.)")
       
       subplot(5,4,2 + j*4 - 4)
       histogram(log(cpCC3))
       set(gca, 'xscale','log')
       title(strcat("CC of CP f=",sprintf('%1.1f',fj)))
       xlabel("log(CC)")
       ylabel("log(Occur.)")
       
       subplot(5,4,3 + j*4 - 4)
       histogram(log(f2fDD3))
       set(gca, 'xscale','log')
       title(strcat("DD of F2F f=",sprintf('%1.1f',fj)))
       xlabel("log(DD)")
       ylabel("log(Occur.)")
       
       subplot(5,4,4 + j*4 - 4)
       histogram(log(cpDD3))
       set(gca, 'xscale','log')
       title(strcat("DD of CP f=",sprintf('%1.1f',fj)))
       xlabel("log(DD)")
       ylabel("log(Occur.)")
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       figure(4 + i*5 - 5) % algo 4 frontier sampling graph
       subplot(5,4,1 + j*4 - 4)
       histogram(log(f2fCC4))
       set(gca, 'xscale','log')
       title(strcat("CC of F2F f=",sprintf('%1.1f',fj)))
       xlabel("log(CC)")
       ylabel("log(Occur.)")
       
       subplot(5,4,2 + j*4 - 4)
       histogram(log(cpCC4))
       set(gca, 'xscale','log')
       title(strcat("CC of CP f=",sprintf('%1.1f',fj)))
       xlabel("log(CC)")
       ylabel("log(Occur.)")
       
       subplot(5,4,3 + j*4 - 4)
       histogram(log(f2fDD4))
       set(gca, 'xscale','log')
       title(strcat("DD of F2F f=",sprintf('%1.1f',fj)))
       xlabel("log(DD)")
       ylabel("log(Occur.)")
       
       subplot(5,4,4 + j*4 - 4)
       histogram(log(cpDD4))
       set(gca, 'xscale','log')
       title(strcat("DD of CP f=",sprintf('%1.1f',fj)))
       xlabel("log(DD)")
       ylabel("log(Occur.)")
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       figure(5 + i*5 - 5) % algo 5 snow ball expansion sampling graph
       subplot(5,4,1 + j*4 - 4)
       histogram(log(f2fCC5))
       set(gca, 'xscale','log')
       title(strcat("CC of F2F f=",sprintf('%1.1f',fj)))
       xlabel("log(CC)")
       ylabel("log(Occur.)")
       
       subplot(5,4,2 + j*4 - 4)
       histogram(log(cpCC5))
       set(gca, 'xscale','log')
       title(strcat("CC of CP f=",sprintf('%1.1f',fj)))
       xlabel("log(CC)")
       ylabel("log(Occur.)")
       
       subplot(5,4,3 + j*4 - 4)
       histogram(log(f2fDD5))
       set(gca, 'xscale','log')
       title(strcat("DD of F2F f=",sprintf('%1.1f',fj)))
       xlabel("log(DD)")
       ylabel("log(Occur.)")
       
       subplot(5,4,4 + j*4 - 4)
       histogram(log(cpDD5))
       set(gca, 'xscale','log')
       title(strcat("DD of CP f=",sprintf('%1.1f',fj)))
       xlabel("log(DD)")
       ylabel("log(Occur.)") 
    end    
end
