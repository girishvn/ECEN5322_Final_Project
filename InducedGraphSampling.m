function [Vs, As] = InducedGraphSampling(G, p, k)
%INDUCED Summary of this function goes here
%   Detailed explanation goes here
    Vs = zeros(1, k);
    
    %Select K vertices based on P from G(E,V)
    [x, y] = size(G);
    
    if p == 'uniform'
        %uniform sampling without replacement
        Vs = randperm(x, k);
        
    elseif p == 'proportional'
        %calculate weights based on vertex degree's
        Vweights = sum(G);
        Volume = sum(Vweights);
        w = Vweights./Volume;
        
        %weighted sampling without replacement
        n = 1:x;
        Vs = zeros(1,k);
        for i = 1:k
            Vs(i) = randsample(n,1,true,w);
            w(n == Vs(i)) = 0;
        end        
    end
    
    %Select Edges 
    As = zeros(x);
    for m = 1:k
        for n = 1:k
            v = Vs(m);
            u = Vs(n);
            As(v,u) = G(v, u);
        end
    end            
    
end

