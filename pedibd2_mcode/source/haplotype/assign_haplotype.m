function [inheritance, haplotype, ...
    consistency, inconsistency, error] ...
    ...
    = assign_haplotype ...
    (assignment, family, markers1, ...
    ...
    family_genotype, listIND2, markers2)

    error = 0;
    inheritance = [];
    haplotype = [];
    consistency = [];
    inconsistency = [];
        
    [ancestralSharing, boundaries, error] = check_allsegments(assignment, family);
    if( error ~= 0 )
        error = 1;
        return;
    end
    [geno, map, error] ...
        = generate_markermap(boundaries, family, family_genotype, ...
        markers1, markers2, ...
        listIND2);
    if( error ~= 0 )
        error = 1;
        return;
    end
    
    nIND = size(family.structure, 1);
    nLOC = length(map);
    [nSEG, c] = size(boundaries);
    if( c ~= 2 )
        error = 1;
        disp('error in generating chromosome boundaries');
        return;
    end
    if( any(map > nSEG) )
        error = 1;
        disp('error in generating chromsome boundaries');
        return;
    end
    
    display(' ');
    display(['haplotyping ...']);
    
    time = cputime;
        
    global complexity;
    global debug_mode;
    complexity = 0;    
    
    inheritance = ancestralSharing(map(map > 0), 1:nIND, 1:2);
    geno = geno(map > 0, :, :);    
    
%     [haplotype, consistency, typing_error, error] = haplotype_allLOCI(inheritance, geno, map);
    [haplotype, consistency, error] = consistency_allLOCI(inheritance, geno, family);
    
    if( ~all(map > 0) )
        temp = int8(zeros(nLOC, nIND, 2));
        temp(map > 0, :, :) = haplotype;
        haplotype = temp;
        temp = zeros(nLOC, 1);
        temp(map > 0) = consistency;
        consistency = temp;        
    end
    
    inconsistency = zeros(nSEG,1);
    for i = 1:nSEG
        inconsistency(i) = nnz(consistency(map == i) < nIND)./nnz(map == i);
    end
    
    if( debug_mode == 1 )
        ratio = nnz(consistency==nIND & map > 0)/nnz(map > 0);
        disp(['haplotype consistency: ', num2str(ratio*100, '%.2f'), ' % of loci']);
    end
    
    display(['haplotype reconstruction time: ', num2str(cputime - time), ' seconds']);
    
    clear global complexity;
    
end






