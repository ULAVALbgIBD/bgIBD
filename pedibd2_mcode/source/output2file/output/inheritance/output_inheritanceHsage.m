function [error] = output_inheritanceHsage(fid, ...
    family_range, ...
    pedigree_info, ...
    assignment, parameters)

error = 0;

if( isempty(family_range) || isempty(assignment) )
    error = 1;
    return;
end

nfam = length(family_range);

if( nfam <= 0 )
    error = 1;
    return;
end

len = length(assignment);
if( len <= 0 || len ~= nfam )
    error = 1;
    disp('not all families assigned');
    return;
end

sampled_markerlist = parameters.sampled_markerlist;
[nmarkers, cols] = size(sampled_markerlist);
if( nmarkers <= 0 || cols ~= 2 )
    error = 1;
    disp('error in marker list');
    return;
end

chr = parameters.chr;
map = sampled_markerlist(1:nmarkers,2);


fprintf(fid, '*****       inheritance      *****\n\n\n\n');


for i = 1:nfam

    family = family_range{i}.structure;
    family_id = family_range{i}.family_id;
    intervals = assignment{i}.intervals;
    alleles = assignment{i}.alleles_all;
    [nseg, nIND, ~] = size(alleles);
    
    fprintf(fid, 'family %6d\n', family_id);
    fprintf(fid, 'chr\t  position\t          \t\t');
    for k = 1:nIND
        if( family(k,7) == 1 )
            fprintf(fid, 'ker%d\t\t', pedigree_info(family(k,8),8));
        end
    end
    fprintf(fid, '\n');
        
    for j = 1:nseg
        s = intervals(j,1);
        t = intervals(j,2);
        fprintf(fid, '%3d\t', chr);
        fprintf(fid, '%10d\t%10d\t\t', map(s), map(t)); 
        for k = 1:nIND
            if( family(k,7) == 1 )
                fprintf(fid, '%3d/%3d\t', 77 + alleles(j,k,1:2));
            end
        end
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n');
    
end

end









































