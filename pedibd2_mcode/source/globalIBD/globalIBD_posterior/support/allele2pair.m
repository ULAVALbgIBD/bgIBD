function pair = allele2pair(alleles)

    pair = [];
    [nind, c] = size(alleles);
    if( nind <= 0 || c ~= 2 )
        return;
    end    
    
    a1 = repmat(alleles(1:nind,1), [1,nind]);
    a2 = repmat(alleles(1:nind,2), [1,nind]);
    b1 = a1';
    b2 = a2';
    
    p1 = (a1==b1)+(a2==b2);
    p2 = (a1==b2)+(a2==b1);
    pair = max(p1,p2);
    
end