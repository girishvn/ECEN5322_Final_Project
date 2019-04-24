function [recoveredNodes, nr] = SIRModelSim(A, B, u, dt)

[m,n] = size(A); % get dimensions of A

% init vars
q = u*dt; % prob to recover with
v0 = randi(m); % random start node

% Init node info arrays
S = ones(m,1); % susceptible nodes
I = zeros(m,1); % infected nodes
R = zeros(m,1); % recovered nodes

S(v0) = 0; % initial infected node not susceptible
I(v0) = 1; % initial infected node is infected

while sum(I) ~= 0 % continue while infected nodes exist
    St = find(S); % find susceptible nodes 
    It = find(I); % find infected nodes
    NIt = find(I == 0); % find uninfected nodes 
    NRt = find(R == 0); % find unrecovered nodes
    
    for i = 1:length(It) % iterate through infected nodes
        
        v = I(i); % iterate through infected nodes v
        Nv = find(A(v,:) > 0); % Nv neighbor hood of v
        Nv = intersect(Nv, St); % Nv limited to susceptible nodes
        Nv = intersect(Nv, NRt); % Nv limited to non-recovered nodes
        Nv = intersect(Nv, NIt); % Nv limited to non-infected nodes 
        
        for j = 1:length(Nv) % iterate through infectable neighbors 
           u = Nv(j); % u in Nv
           p = B*A(v,u)*dt; % probability of getting infected 
           
           infected = binornd(1,p); % u infected with prob p
           I(u) = infected;
           S(u) = ~infected;
        end
        
        recovered = binornd(1,q); % v recovered with prob q
        R(v) = recovered;
        I(v) = ~recovered;
    end    
end

recoveredNodes = find(R);
nr = length(recoveredNodes)/n;
end

