function [recoveredNodes] = Karros_SIR(A, B, u, dt, vo, V)

    %Initialize Necessary Constants
%     nTimeSteps = T/dt; %Create Timestep
    nTimeSteps = 100;
    q = u.*dt; %Calculate probability infected person recovers
    noV = length(V); %Get Total Number of Vertices
    
    %Initialize array structure to hold characteristics of each person
    x.infected = 0;
    x.recovered = 0;
    x.susceptible = 1;
    person = repmat(x, 1, noV);
    
    %Go through all initial infected nodes and update the array structure
    for k = 1:length(vo)
        person(vo(k)).infected = 1;
        person(vo(k)).susceptible = 0;
    end
    
    %Begin Simulation
    for t = 1:nTimeSteps %run for nTimeSteps
        for l=1:noV %loop through each person (l)
            if(person(l).infected == 1) %check if person is infected
                for i=1:noV %loop through the infected nodes neighbors (j)
                    if(A(i,l) >= 1) %node (j) is a neighbor of (i)
                        if(person(i).susceptible == 1) %check if the neighbor is susceptible to the disease
                            if(rand() < B.*A(i,l).*dt) %infect the person with probability p
                                person(i).infected = 1;
                                person(i).susceptible = 0;
                                vo = [vo, i]; %add infected person to list of infected people (vo)
                            end
                        end
                    end
                end %Looped through all the infected person's neighbors
                %Now check if the infected person will recover
                if(rand() < q) %infected person will recover w/ probability q
                    person(l).infected = 0;
                    person(l).recovered = 1;
                    person(l).susceptible = 0;
                end
            end
        end %looped through every infected person for this time step
    end %Finished Simulation
    
    %Generate and return a list of all the recovered nodes
    recoveredNodes = [];
    for l=1:noV
        if(person(l).recovered == 1)
            recoveredNodes = [recoveredNodes, l];
        end
    end
end