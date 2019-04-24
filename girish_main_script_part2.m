%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Girish Narayanswamy
% ECEN 5322: Higher-Dimensional Datasets
% Final Project: Assignment 7 - 14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the plotting script for the histograms
% for degree density and clustering coefficient 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f2f_A_files = ["f2f_adj_InVS13.mat", "f2f_adj_InVS15.mat", "f2f_adj_LH10.mat", "f2f_adj_LyonSchool.mat", "f2f_adj_SFHH.mat", "f2f_adj_Thiers13.mat"];
cp_A_files = ["cp_adj_InVS13.mat", "cp_adj_InVS15.mat", "cp_adj_LH10.mat", "cp_adj_LyonSchool.mat", "cp_adj_SFHH.mat", "cp_adj_Thiers13.mat"];
dataSets = ["InVS13", "InVS15", "LH10", "LyonSchool", "SFHH", "Thiers13"];

for i = 1:length(f2f_A_files) % iterate through data sets
    
    %Original, before sampling
    
    % load in f2f A
    f2fA =  load(f2f_A_files(i));
    f2fA = f2fA.data;
    
    % load in cp A
    cpA = load(cp_A_files(i));
    cpA = cpA.data;
    
    B = 4*10^(-4);
    k = [1, 4, 16, 64, 256];
    u = B/k;
    
    for j = 1:length(u)
        uj = u(j);
        r100 = zeros(100,1);
        
        for k = 1:100
            
            
        end
        
    end
    
    
    
end