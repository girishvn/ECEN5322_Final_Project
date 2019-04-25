%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clint Olsen
% ECEN 5322: Higher-Dimensional Datasets
% Final Project
% SIR Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Clint_SIR.m: Runs the SIR model on a given adjacency matrix
%
% - Inputs: Adjacency matrix (A), Infection rate (Beta), Recovery rate (mu)
%
% - Outputs: Number of recovered nodes and count distributions for
% infected, susceptible and recovered nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [numRecovered, infectedDistribution, susceptibleDistribution, recoveredDistribution] = Clint_SIR(A, Beta, mu)

    % Time step
    % Determine so probability is normalized between 0 and 1 for
    med = median(A(find(A~=0)));
    dt = (1/Beta)*(1/med)*0.001;

    %Generate random seedm for 1st node to be infected
    numInfected = 1;
    totalInfected = 1;
    %infected(numInfected) = randi([1 size(A,1)],1,1);
    infected = zeros(1,size(A,1));
    infected(1,randi([1 size(A,1)],1,1)) = 1;
    
    numRecovered = 0;
    plotIndex = 1;
    
    numSusceptible = size(A,1) - 1;
    susceptible = ones(1,size(A,1));
    susceptible(1,find(infected == 1)) = 0;
    while(numInfected > 0) % This is when the infection dies out
        
        infectedList = find(infected == 1); % Pick an infected node
        infectedIndex = infectedList(1);
           for j=2:size(A,1)    % Loop through neighbors of that infected node
               if(A(infectedIndex,j) > 0 && susceptible(j) == 1) % If index is a neighbor and is susceptible
                   if(binornd(1,Beta*A(infectedIndex,j)*dt) >= 1) % Node will become infected
                        numInfected = numInfected + 1;
                        infected(j) = 1; % Infect the Node
                        susceptible(j) = 0; % Infected implies no longer susceptible
                        numSusceptible = numSusceptible - 1;
                        
                        % Collect counts for returned data sets (for plots)
                        susceptibleDistribution(plotIndex) = numSusceptible;
                        infectedDistribution(plotIndex) = numInfected;
                        recoveredDistribution(plotIndex) = numRecovered;
                        plotIndex = plotIndex + 1;

                    end
               end
           end

            %Probability current infected node will transition to recovered
            % The current calculation for dt makes this value very small,
            % so the likelihood of all nodes getting infected is much
            % higher due to mu decreasing with increased k (increasing the
            % strength of the spread) rho becomes larger and larger than
            % 1
            if(binornd(1,mu*dt) >= 1) % Node will become recovered
                infected(infectedIndex) = 0;
                numRecovered = numRecovered + 1; % Node Recovers
                numInfected = numInfected - 1;
                susceptible(infectedIndex) = 0;

                % Collect counts for returned data sets (for plots)
                infectedDistribution(plotIndex) = numInfected;
                recoveredDistribution(plotIndex) = numRecovered;
                susceptibleDistribution(plotIndex) = numSusceptible;
                plotIndex = plotIndex + 1;
                
            end

    end
    
end