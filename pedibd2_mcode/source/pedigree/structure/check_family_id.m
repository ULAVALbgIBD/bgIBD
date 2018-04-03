function error = check_family_id(pedigree)


% check to see whether all pedigree ids are forming consecutive indexing

error = 0;

[rows cols] = size(pedigree);

if( rows <=0 || cols < 7 )
    error = 1;
    disp('error in pedigree structure');
    return;
end

fids = unique(pedigree(:,1));

for f = fids'
    f_index = find(pedigree(:,1) == f);
    if( isempty(f_index) )
        error = 1;
        disp('error in pedigree structure');
        return;
    end
    if( length(unique(f_index)) ~= length(f_index) )
        error = 1;
        disp('error in pedigree indexing');
        return;
    end
    if( max(f_index) - min(f_index) + 1 ~= length(f_index) )
        error = 1;
        disp('error in pedigree indexing');
        return;
    end
    for i = 1:length(f_index)
        if( pedigree(f_index(i),2) ~= i )
            error = 1;
            disp('error in pedigree indexing');
            return;
        end
    end
end


end





