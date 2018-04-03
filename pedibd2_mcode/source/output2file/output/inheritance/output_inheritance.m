function [error] = output_inheritance(fid, all_families, all_intervals, assignment, markers)
    
    error = 0;

    nfam = length(all_families);
    if( nfam <= 0 )
        error = 1;
        disp('error in family structures');
        return;
    end
    
    if( nfam ~= length(assignment) )
        error = 1;
        disp('error in allele assignment');
        return;
    end    
    
    [nseg c] = size(all_intervals);
    if( nseg <= 0 || c ~= 2 )
        error = 1;
        disp('error in intervals');
        return;
    end
    
    [nmarker c] = size(markers);
    if( nmarker <= 0 || c ~= 2 )
        error = 1;
        disp('error in marker positions');
        return;
    end
    
    if( any(all_intervals(1:nseg,2) < all_intervals(1:nseg,1)) )        
        error = 1;
        disp('error in intervals');
        return;
    end
    
    if( nseg > 1 && any(all_intervals(2:nseg,1) - all_intervals(1:nseg-1,2) ~= 1) )
        error = 1;
        disp('error in intervals');
        return;
    end
    
    if( all_intervals(1,1) ~= 1 || all_intervals(nseg,2) ~= nmarker )
        error = 1;
        disp('error in intervals');
        return;
    end
    
    % map global intervals to local intervals in each family
    index = zeros(nfam, nseg);
    for i = 1:nfam
        intervals = assignment{i}.intervals(:,1:2);
        for j= 1:nseg
            a = all_intervals(j,1);
            b = all_intervals(j,2);
            bit = (intervals(:,1) <= a) & (intervals(:,2) >= b);
            ind = find(bit);
            if( length(ind) ~= 1 )
                error = 1;
                return;
            end
            index(i,j) = ind;
        end
    end
    
    fprintf(fid, '*****        output format:              *****\n');
    fprintf(fid, '*****        individual information      *****\n');
    fprintf(fid, '*****        paternal haplotype          *****\n');
    fprintf(fid, '*****        maternal haplotype          *****\n\n\n\n');
    
    for i = 1:nfam
        
        if( error ~= 0 )
            continue;
        end
        
        family = all_families{i}.structure;
        [nind, cols] = size(family);
        if( nind <= 0 || cols < 12 )
            error = 1;
            return;
        end
        
        alleles_all = assignment{i}.alleles_all;
              
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
            
            fprintf(fid, '%d\t%d\t%d\t%d\t%d\t%d\n', f, id, father_id, mother_id, sex, status);
            % let the proceeding line be the paternal line
            for k = 1:nseg
                fprintf(fid, '%3d ', alleles_all{index(i,k)}(self,1));                
            end
            fprintf(fid, '\n');
            for k = 1:nseg
                fprintf(fid, '%3d ', alleles_all{index(i,k)}(self,2));              
            end
            fprintf(fid, '\n');
        end
        
    end
    
end











