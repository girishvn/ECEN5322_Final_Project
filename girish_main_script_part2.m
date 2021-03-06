%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Girish Narayanswamy
% ECEN 5322: Higher-Dimensional Datasets
% Final Project: Assignment 7 - 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the plotting script for the histograms
% for degree density and clustering coefficient 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f2f_A_files = ["A_InVS13.mat", "A_InVS15.mat", "A_LH10.mat", "A_LyonSchool.mat", "A_SFHH.mat", "A_Thiers13.mat"];
dataSets = ["InVS13", "InVS15", "LH10", "LyonSchool", "SFHH", "Thiers13"];

for i = 1:length(f2f_A_files) % iterate through data sets
    
    % init 2 figures per data set
    figure(2*i - 1);
    sgtitle(strcat("Distribution of Recovered Nodes (100 Sims.) for f2f ", dataSets(i)));
    
    figure(2*i);
    sgtitle(strcat("Stats. for Epidemics with nr > 20% (100 Sims.) for f2f ", dataSets(i)));
    
    %Original, before sampling
    % load in f2f A
    A =  load(f2f_A_files(i));
    A = A.Z;
    
    [m,n] = size(A);
    
    B = (4e-4); % beta values
    k = [2, 4, 6, 8, 10]; % values of k / p0
    avgDeg = round(sum(sum(A))/n); % calculate average degree
    p0 = k*avgDeg; % value of p0
    u = B./k.*100; % mu array
    
    fracNr20Arr = zeros(length(u),1); % build arrays to hold plot for Assignment 8
    meanNr20Arr = zeros(length(u),1); % build arrays to hold plot for Assignment 9
    
    for j = 1:length(u)
        uj = u(j);
        numRcvrd100 = zeros(1,100);
        parfor sim = 1:100
            [numRcvrd, I, S, R] = Clint_SIR(A, B, uj);
            numRcvrd100(sim) = numRcvrd;
        end
        
        nrArr = numRcvrd100/n; % converet recovered node distirbution to nr distribution 
        numEpNr20 = sum(nrArr > 0.2); % number of epidemics with nr > 0.2
        fracNr20Arr(j) = numEpNr20/n; % fraction of epidemics of nr > 0.2
        meanNr20Arr(j) = sum(numRcvrd100(nrArr > 0.2))/numEpNr20; % mean recovered nodes for Epidemics with nr > 0.2
        
        figure(2*i - 1);
        subplot(2,3,j);
        histogram(numRcvrd100,n+20); % adjust bin size of histogram to number of nodes 
        title(strcat("p0 = ", num2str(p0(j), '%d')))
        xlabel("Recovered Nodes");
        ylabel("Occur.")
    end
    
    figure(2*i);
    subplot(2,1,1);
    plot(p0,fracNr20Arr);
    title("Fraction of Epidemics with nr > 20%")
    xlabel("p0 = B<n>/u");
    ylabel("Frac. of Epidemics")
    
    subplot(2,1,2);
    plot(p0,meanNr20Arr);
    title("Mean Recovered Nodes for Epidemics with nr > 20%")
    xlabel("p0 = B<n>/u");
    ylabel("Recovered Nodes")
    
end