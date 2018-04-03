function [error] = output_recombination(fid, family_range, assignment, parameters)

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

chr = parameters.chr;
sampled_markerlist = parameters.sampled_markerlist;
if( chr < 0 )
    error = 1;
    disp('error in marker list, negative chromosome code');
    return;
end
[nmarkers, cols] = size(sampled_markerlist);
if( nmarkers <= 0 || cols ~= 2 )
    error = 1;
    disp('error in marker list');
    return;
end
map(1:nmarkers) = sampled_markerlist(1:nmarkers,2);

for i = 1:nfam
    if( isempty(family_range{i}) || isempty(assignment{i}) )
        error = 1;
        return;
    end
    family = family_range{i}.structure;
    [nind, fields] = size(family);
    if( nind <= 0 || fields < 12 )
        error = 1;
        return;
    end
    recombination = assignment{i}.recombination.list;
    [r, c] = size( recombination );
    if( r > 0 )
        if( c ~= 4 )
            error = 1;
            return;
        end
        for j = 1:r
            if( any(recombination(j,1:2) > nmarkers) )
                error = 1;
                disp('error in recombination positions, out of bound')
                return;
            end
            if( any(recombination(j,3:4) > nind) )
                error = 1;
                disp('error in recombination gamete genitor, out of bound');
                return;
            end
        end
    end
end

fprintf(fid, '*****       recombination positions      *****\n\n\n\n');

fprintf(fid, '                                     gamete from \t\t  family\n');
fprintf(fid, 'chr        crossover position        parent      \t\t');
fprintf(fid, '\n');

for i = 1:nfam

    family = family_range{i}.structure;
    family_id = family_range{i}.family_id;
    recombination = assignment{i}.recombination.list;
    [r, c] = size( recombination );
    
    for j = 1:r
        fprintf(fid, '%3d   ', chr);
        s = recombination(j,1);
        t = recombination(j,2);
        p1 = recombination(j,3);
        p2 = recombination(j,4);
        fprintf(fid, '%10d ~ %10d        ', map(s), map(t)); 
        if( p1 == 0 && p2 == 0 )
            error = 1;
            return;
        end
        if( p1 ~= 0 && p2 == 0 )
            fprintf(fid, '%6d          ', family(p1,8));
        end
        if( p1 == 0 && p2 ~= 0 )
            fprintf(fid, '%6d          ', family(p2,8));
        end
        if( p1 ~= 0 && p2 ~= 0 )
            fprintf(fid, '%6d or %6d', family(p1,8), family(p2,8));
        end
        fprintf(fid, '\t\t  %6d', family_id);
        fprintf(fid, '\n');
    end
    
end

end









































