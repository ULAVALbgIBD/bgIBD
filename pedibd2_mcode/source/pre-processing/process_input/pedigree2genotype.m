function [expanded error] = pedigree2genotype(pedigree, all_data)

    error = 0;
    [nind, cols] = size(pedigree);
    if( nind <= 0 || cols < 12 )
        error = 1;
        disp('error in pedigree data');
        return;
    end
    
    genotype = all_data.genotype;
    [nloc, fields, c] = size(genotype);
    if( nloc <= 0 || fields < 0 || c ~= 2 )
        error = 1;
        disp('error in genotype data');
        return;
    end
    header = all_data.header;
    [d1, d2] = size(header);
    if( d1 ~= fields || d2 ~= 6 )
        error = 1;
        disp('error in genotype data');
        return;
    end
    
    expanded = int8(zeros(nloc, nind, 2));
    map(1:nind) = 0;
    for i = 1:nind
        family_id = pedigree(i,12);
        original_id = pedigree(i,8);
        genotyped = pedigree(i,7);
        index = find( header(:,1) == family_id & header(:,2) == original_id );
        if( length(index) == 1 )
            map(i) = index;
            temp_genotype = genotype(1:nloc, index, 1:2);
            typed_ratio = nnz(temp_genotype)/(nloc * 2);
            if( typed_ratio < 0.5 && genotyped == 1 )
                disp(['warning: family ', num2str(family_id), ': individual ', num2str(original_id), ', marked genotyped but missing ', num2str((1-typed_ratio)*100, '%.2f'), '% genotypes']);
            end
            if( typed_ratio > 0 && genotyped == 0 )
                disp(['warning: family ', num2str(family_id), ': individual ', num2str(original_id), ', marked ungenotyped, genotypes ignored']);
            end
            
            if( genotyped == 1 )
                expanded(1:nloc, i, 1:2) = genotype(1:nloc, index, 1:2);
            else
                expanded(1:nloc, i, 1:2) = 0;
            end            
        else
            if( length(index) > 1 )
                disp(['family ', num2str(family_id), ': multiple genotype lines for individual ', num2str(original_id)]);
                error = 1;
                return;
            end
            if( length(index) == 0 && genotyped == 1 )
                disp(['family ', num2str(family_id), ': cannot find genotype for indivudal ', num2str(original_id)]);
                error = 1;
                return;
            end
            map(i) = 0;
            expanded(1:nloc, i, 1:2) = 0;
        end
    end    
end









