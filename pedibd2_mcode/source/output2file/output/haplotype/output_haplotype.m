function [error] = output_haplotype(fid, all_families, haplotype, markers)
    
    error = 0;

    nfam = length(all_families);
    if( nfam <= 0 )
        error = 1;
        return;
    end
    
    if( nfam ~= length(haplotype) )
        error = 1;
        return;
    end
    
    if( ndims(markers) ~= 2 )
        error = 1;
        return;
    end
    
    [nloc c] = size(markers);
    if( nloc <= 0 || c ~= 2 )
        error = 1;
        return;
    end
    
    fprintf(fid, '*****        output format:              *****\n');
    fprintf(fid, '*****        individual information      *****\n');
    fprintf(fid, '*****        paternal haplotype          *****\n');
    fprintf(fid, '*****        maternal haplotype          *****\n\n\n\n');
    
    for i = 1:length(all_families)
        error = check_haplotype1family(all_families{i}, haplotype{i}, markers);
        if( error ~= 0 )
            continue;
        end
        
        family = all_families{i}.structure;
        [nind, cols] = size(family);
        if( nind <= 0 || cols < 12 )
            error = 1;
            return;
        end
        
        for j = 1:nind
            self = j; % local index
            id = family(self,8);
            f = family(self,12);
            father = family(self,3);
            if( father ~= 0 )
                father_id = family(father,8);
            else
                father_id = 0;
            end
            mother = family(self,4);
            if( mother ~= 0 )
                mother_id = family(mother,8);  
            else
                mother_id = 0;
            end
            sex = family(self,5);
            status = family(self,6);
            result = reshape(haplotype{i}(1:nloc,self,1:2), [nloc,2]);
            fprintf(fid, '%d\t%d\t%d\t%d\t%d\t%d\t', f, id, father_id, mother_id, sex, status);
            % let the proceeding line be the paternal line
            fprintf(fid, '\n');
            fprintf(fid, '%d ', result(1:nloc,1));
            fprintf(fid, '\n');
            fprintf(fid, '%d ', result(1:nloc,2));
            fprintf(fid, '\n');
        end
    end
end











