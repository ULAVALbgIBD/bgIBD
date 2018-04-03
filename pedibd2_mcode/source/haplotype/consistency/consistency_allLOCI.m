function [hap, consistency, error] = consistency_allLOCI(inheritance, genotype, family)

    error = 0;
    
    [nLOC, nIND, ~] = size(inheritance);
    
    hap = int8(zeros(nLOC,nIND,2));
    if( nIND < 2^7 )
        consistency = int8(zeros(nLOC, 1));
    else
        consistency = zeros(nLOC, 1);
    end

    if( nnz(inheritance(1:nLOC,family.family_range,1:2) == 0) > 0 )
        error = 1;
        disp('not all inheritance generated');
        return;
    end
    
    inh = inheritance;
    
    
    offset = false(nLOC, nIND, 2);
    
    allele1 = repmat(all(genotype(1:nLOC,1:nIND,1:2) == 1,3), [1,1,2]);
    allele2 = repmat(all(genotype(1:nLOC,1:nIND,1:2) == 2,3), [1,1,2]);
    hete = all(genotype(1:nLOC,1:nIND,1:2) ~= 0,3) & genotype(1:nLOC,1:nIND,1) ~= genotype(1:nLOC,1:nIND,2);
    homo1 = all(genotype(1:nLOC,1:nIND,1:2) == 1,3);
    homo2 = all(genotype(1:nLOC,1:nIND,1:2) == 2,3);
    homo = homo1 | homo2;
    
    alleleCOUNT = int8(zeros(nLOC,nIND,2));

    % ungenotyped individual is not assigned inheritance
    % if not imputed
    % left as "0"
    conflict = false(nLOC, 1);
    for i = 1:nIND

        disp(['          haplotyping individual: ', num2str(family.structure(i,8)), ' ...']);
        
        pmat = repmat(inh(1:nLOC,i,1), [1,nIND,2]);
        mmat = repmat(inh(1:nLOC,i,2), [1,nIND,2]);
        mbit = inh == mmat & mmat ~= 0;
        
        % homozygous i, heterzygous i, two sets or one set
        joined = inh(1:nLOC,i,1) == inh(1:nLOC,i,2);
        hetei_joined = hete(1:nLOC,i) & joined;
        hetei_indep = repmat(hete(1:nLOC,i) & ~joined, [1,nIND,2]);
        homoi_joined = homo(1:nLOC,i) & joined;
        homoi_indep = repmat(homo(1:nLOC,i) & ~joined, [1,nIND,2]);
        
        % already same set, check consistency
        cbit = hetei_joined & offset(1:nLOC,i,1) == offset(1:nLOC,i,2);
        conflict(cbit) = true;
        cbit = homoi_joined & offset(1:nLOC,i,1) ~= offset(1:nLOC,i,2);
        conflict(cbit) = true;
        
        % join two sets at hete individual
        repbit = mbit & hetei_indep;
        inh(repbit) = pmat(repbit);
        % change offset when merging
        repbitINV = repbit & repmat(offset(1:nLOC,i,1) == offset(1:nLOC,i,2), [1,nIND,2]);
        offset(repbitINV) = ~offset(repbitINV);
        % join two sets at homo individual
        repbit = mbit & homoi_indep;
        inh(repbit) = pmat(repbit);
        % change offset when merging
        repbitINV = repbit & repmat(offset(1:nLOC,i,1) ~= offset(1:nLOC,i,2), [1,nIND,2]);
        offset(repbitINV) = ~offset(repbitINV);
        
        % check fixed homozygous genotype consistency in merged set
        pbit = inh == pmat & pmat ~= 0;
        pbitp1 = any( any(pbit & offset & allele1, 3), 2);
        pbitp2 = any( any(pbit & offset & allele2, 3), 2);
        pbitn1 = any( any(pbit & ~offset & allele1, 3), 2);
        pbitn2 = any( any(pbit & ~offset & allele2, 3), 2);
        
        cbit = (pbitp1 & pbitp2) | (pbitp1 & pbitn1) | (pbitp2 & pbitn2) | (pbitn1 & pbitn2);
        conflict(cbit) = true;
        
        % assign haplotype
        hapbit = repmat(pbitp1 | pbitn2, [1, nIND, 2]);
        hap(pbit & offset & hapbit) = 1;
        hap(pbit & ~offset & hapbit) = 2;
        hapbit = repmat(pbitp2 | pbitn1, [1, nIND, 2]);
        hap(pbit & offset & hapbit) = 2;
        hap(pbit & ~offset & hapbit) = 1;
        
        bitpm = inheritance == repmat(inheritance(1:nLOC,i,1), [1,nIND,2]) | inheritance == repmat(inheritance(1:nLOC,i,2), [1,nIND,2]);
        bitALLELE = repmat(homo1(1:nLOC,i), [1,nIND,2]) & bitpm;
        alleleCOUNT(bitALLELE) = alleleCOUNT(bitALLELE) + 1;
        bitALLELE = repmat(homo2(1:nLOC,i), [1,nIND,2]) & bitpm;
        alleleCOUNT(bitALLELE) = alleleCOUNT(bitALLELE) - 1;
        
    end 
    

    
    consistency(conflict) = 0;
    consistency(~conflict) = nIND;
    
    % some missing genotypes can be assigned during above process
    % set "0" before voting
    bitVOTE = repmat(conflict, [1,nIND,2]);
    hap(bitVOTE) = 0;
    % set all homozygous loci
    hap(bitVOTE & allele1) = 1;
    hap(bitVOTE & allele2) = 2;
    % set heterozygous loci
    % 1 is on left, 2 is on right
    bit1p = alleleCOUNT(1:nLOC,1:nIND,1) >= alleleCOUNT(1:nLOC,1:nIND,2);
    % 2 is on left, 1 is on right
    bit2p = ~bit1p;
    bit1 = cat(3,bit1p,bit2p);
    bit2 = cat(3,bit2p,bit1p);
    bitHETE = bitVOTE & repmat(hete, [1,1,2]);
    hap(bitHETE & bit1) = 1;
    hap(bitHETE & bit2) = 2;
    
    % verified the same as haplotype_allLOCI, on July 24, 2012
    % on ceufam data

end

    
    
    
    
    