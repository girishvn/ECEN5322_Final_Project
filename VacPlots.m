%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clint Olsen
% ECEN 5322: Higher-Dimensional Datasets
% Final Project
% Plots for SIR Model with Vaccination (11-13)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - VacPlots.m: Creates plots for problems 11,12,13
%
% - Inputs: cp: "true" - run on subsampled co-presence
%               "false" - run on face-to-face matrix
%
% - Outputs: N/A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function VacPlots(cp)
    f2f_A_files = ["A_LH10", "A_InVS13.mat", "A_InVS15.mat", ...
                    "A_LyonSchool.mat", "A_SFHH.mat", "A_Thiers13.mat"];
    cp_A_files = ["A_pres_LH10", "A_pres_InVS13.mat", "A_pres_InVS15.mat", ...
                    "A_pres_LyonSchool.mat", "A_pres_SFHH.mat", "A_pres_Thiers13.mat"];
    f2f_file_names = ["LH10", "InVS13", "InVS15", "LyonSchool", "SFHH", "Thiers13"];


    % Simulation Parameters
    Beta = 4e-4;
    mu = 100*Beta;
    numSims = 100;

    % Index 1 is without vac and 2 is with random, and 3 is with high degree
    % vac
    recovered = zeros(1,3);
    recoveredWithVac = 0;
    recoveredWithVacNoRand = 0;
    recoveredWithoutVac = 0;
    
    for i = 1:size(f2f_A_files,2) - 5

        if(cp == "false") % Run on Face2Face
            A = load(f2f_A_files(i));
            A = A.Z;
            N = size(A,1);
        elseif(cp == "true") %Run on Co-presence subsampled with Met. Hastings
            a = load(cp_A_files(i));
            a = a.Z;
            N = size(a,1);
            numEdges = size(find(triu(a) > 0), 1);
            A = metropolisHastingsRW(a, numEdges * .8);
        end
        recovered(1) = 0;
        recovered(2) = 0;
        recovered(3) = 0;
        ratioIndex = 1;
        ratioIndex2 = 1;


        numInfNoVacIndex = 1;
        numInfVacRandIndex = 1;
        numInfVacDegIndex = 1;
        
        numInfNoVac = [];
        numInfVacRand = [];
        numInfVacDeg = [];


        for j = 1:numSims  
            recoveredWithoutVac = (Clint_SIR(A,Beta,mu));
            recoveredWithVac = (Clint_SIR_Vac(A,Beta,mu,"true"));
            recoveredWithVacNotRand = (Clint_SIR_Vac(A,Beta,mu,"false"));


            if(recoveredWithoutVac / N >= .2)
                recovered(1) = recovered(1) + 1;
                numInfNoVac(numInfNoVacIndex) = recoveredWithoutVac;
                numInfNoVacIndex = numInfNoVacIndex + 1;
            end

            if(recoveredWithVac / N >= .2)
                recovered(2) = recovered(2) + 1;
                numInfVacRand(numInfVacRandIndex) = recoveredWithVac;
                numInfVacRandIndex = numInfVacRandIndex + 1;
            end
            if(recoveredWithVacNotRand / N >= .2)
                recovered(3) = recovered(3) + 1;
                numInfVacDeg(numInfVacDegIndex) = recoveredWithVacNotRand;
                numInfVacDegIndex = numInfVacDegIndex + 1;
            end

            if(recoveredWithoutVac / N >= .2 && recoveredWithVac / N >= .2)
                ratio(ratioIndex) = recoveredWithVac / recoveredWithoutVac;
                ratioIndex = ratioIndex + 1;
            end

            if(recoveredWithoutVac / N >= .2 && recoveredWithVacNotRand / N >= .2)
                ratio2(ratioIndex2) = recoveredWithVacNotRand / recoveredWithoutVac;
                ratioIndex2 = ratioIndex2 + 1;
            end  
        end

        % Compute Median absolute distribution
        medRand(i) = (median(numInfVacRand) / median(numInfNoVac));
        medDeg(i) = (median(numInfVacDeg) / median(numInfNoVac));
        subplot(2,2,[1 2])
        bar([recovered]);
        title(sprintf("%s Total Number Infected", f2f_file_names(i)));
        set(gca, 'xticklabel', {sprintf("w/0 Vac"); sprintf("w/ Rand Vac"); sprintf("w/ Vac of High Degrees")})
        subplot(2,2,3)
        h = histogram(ratio);
        morebins(h);
        morebins(h);
        morebins(h);
        %errorbar(1:length(ratio),ratio, var(ratio));
        title(sprintf("%s Ratio between with and without random vaccination\nMean = %f Var = %f Std = %f", f2f_file_names(i),mean(ratio), var(ratio), std(ratio)));
        subplot(2,2,4)
        h = histogram(ratio2);
        morebins(h);
        morebins(h);
        morebins(h);
        %errorbar(1:length(ratio2),ratio2, var(ratio2));
        title(sprintf("%s Ratio between with and without high degree vaccination\nMean = %f Var = %f Std = %f", f2f_file_names(i), mean(ratio2), var(ratio2), std(ratio2)));
        figure;
    end
    
    % Assignment 12 Medians
    subplot(1,2,1);
    bar(medRand);
    set(gca, 'xticklabel', {sprintf("LH10"); sprintf("InVS13"); sprintf("InVS15"); sprintf("LyonSchool"); sprintf("SFHH"); sprintf("Thiers13")})
    title("Median Outbreak Ratio with Random Vaccination");
    subplot(1,2,2);
    bar(medDeg);
    set(gca, 'xticklabel', {sprintf("LH10"); sprintf("InVS13"); sprintf("InVS15"); sprintf("LyonSchool"); sprintf("SFHH"); sprintf("Thiers13")})
    title("Median Outbreak Ratio with High Degree Vaccination");
end
