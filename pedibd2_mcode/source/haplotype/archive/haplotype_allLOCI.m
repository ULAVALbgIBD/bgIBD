function [haplotype, consistency, typing_error, error] = haplotype_allLOCI(inheritance, geno, map)

    error = 0;
    
    [num_loci, num_ind, ~] = size(inheritance);
    
    haplotype = int8(zeros(num_loci, num_ind, 2));
    if( num_ind < 2^7 )
        consistency = int8(zeros(num_loci, 1));
    else
        consistency = zeros(num_loci, 1);
    end
    typing_error = zeros(num_loci, 1);
    
    seg = 0;
    for i = 1:num_loci

        if( map(i) ~= seg )
            seg = map(i);
            disp(['          haplotyping segment: ', num2str(seg), ' ...']);
        end
      
        [haplotype(i, 1:num_ind, 1:2), consistency(i), typing_error(i), ~, error] = assign_locus_flat( reshape(inheritance(i, 1:num_ind, 1:2), [num_ind,2]), reshape(geno(i, 1:num_ind, 1:2), num_ind, 2) );
        if( error ~= 0 )
            disp('error in generating haploid genotypes');
            return;
        end
        
        
        
    end

end
