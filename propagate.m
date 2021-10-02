function Clusters = propagate(U, thresh)


if max(U) < thresh
    return
else
    N = size(U);
    nx = N(1); ny = N(2);
    
    M = BuildConnect(nx,ny);
    
    Clusters = [];
    U = U(:);
    while max(U) > thresh
        CutOff = (U > thresh);        
        OldCluster = 0*CutOff;
        [~, ind] = max(U);
        OldCluster(ind) = 1;
        
        BreakCondition = 0;
        while BreakCondition == 0
            NewCluster = CutOff.*(M*OldCluster);
            NewCluster(NewCluster~=0) = 1;
            if sum(OldCluster)==sum(NewCluster)
                BreakCondition = 1;
            else
                OldCluster = NewCluster;
            end
        end        
        U(NewCluster == 1) = 0;
        Clusters = [Clusters NewCluster];
    end   
end