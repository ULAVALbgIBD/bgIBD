function consistency_HALFgenotype(inheritance, offset, haplotype, genotype)


[nloc, nind, ~] = size(inheritance);

% check half determined sites, unsolved problem
half1 = any(genotype(1:nloc,1:nind,1:2) == 0,3) & any(genotype(1:nloc,1:nind,1:2) == 1,3);
half2 = any(genotype(1:nloc,1:nind,1:2) == 0,3) & any(genotype(1:nloc,1:nind,1:2) == 2,3);    

% loci, ind1, allele1, ind2, allele2, relationship
constraintOR = false(nloc, nind, nind, 4);
conflict = false(nloc,1);

code = false(nloc, nind, 4);

code(half1(1:nloc,1:nind) & offset(1:nloc,1:nind,1) & offset(1:nloc,1:nind,2)) = 1;
code(half1(1:nloc,1:nind) & offest(1:nloc,1:nind,1) & ~offset(1:nloc,1:nind,2)) = 2;
code(half1(1:nloc,1:nind) & ~offset(1:nloc,1:nind,1) & offset(1:nloc,1:nind,2)) = 3;
code(half1(1:nloc,1:nind) & ~offset(1:nloc,1:nind,1) & ~offset(1:nloc,1:nind,2)) = 4;
code(half2(1:nloc,1:nind) & offset(1:nloc,1:nind,1) & offset(1:nloc,1:nind,2)) = 4;
code(half2(1:nloc,1:nind) & offest(1:nloc,1:nind,1) & ~offset(1:nloc,1:nind,2)) = 3;
code(half2(1:nloc,1:nind) & ~offset(1:nloc,1:nind,1) & offset(1:nloc,1:nind,2)) = 2;
code(half2(1:nloc,1:nind) & ~offset(1:nloc,1:nind,1) & ~offset(1:nloc,1:nind,2)) = 1;

% satuate all pairwise relationships
for i = 1:nind
    

    
    pmat = repmat(inh(1:nloc,i,1), [1,nind,2]);
    mmat = repmat(inh(1:nloc,i,2), [1,nind,2]);
    mbit = inh == mmat & mmat ~= 0;    
    
    relationshipI = [1:nloc, inh(1:nloc,i,1), inh(1:nloc,i,2), code(1:nloc)];
    constraintOR(nloc*(nind-1)+1:nloc*nind,1:4) = relationshipI;
    
end




end




