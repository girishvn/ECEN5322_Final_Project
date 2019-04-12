function [As] = InducedGraphSampling(G, p, k)
%INDUCED Summary of this function goes here
%   Detailed explanation goes here

    k = round(k);
    Vs = zeros(1, k);
    
    %Select K vertices based on P from G(E,V)
    [x, y] = size(G);
    
    if p == "uniform"
        %uniform sampling without replacement
        Vs = randperm(x, k);
        
    elseif p == "weighted"
        %calculate weights based on vertex degree's
        Vweights = sum(G);
        Volume = sum(Vweights);
        
        %weighted sampling without replacement
        w = 2.*Vweights./Volume;
        n = 1:x;
        
        isolatedNodes = sum(w == 0);
        if k > x - isolatedNodes
           k = x - isolatedNodes;
        end
        
        for i = 1:k
            Vs(i) = randsample(n,1,true,w);
            w = w(n ~= Vs(i));
            n = setdiff(n,Vs(i));
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
