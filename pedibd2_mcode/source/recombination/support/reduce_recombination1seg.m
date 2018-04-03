
function [new_alleles num_rec error] = reduce_recombination1seg(ref_alleles, alleles, family)

error = 0;
num_rec = 0;
new_alleles = [];

if( isempty(ref_alleles) || isempty(alleles) )
    error = 1;
    disp('invalid global IBD');
    return;
end


[nind, cols] = size(family);
if( nind <= 0 || cols < 12 )
    error = 1;
    disp('error in family structures');
    return;
end

founder(1:nind) = 0;
for i = 1:nind
    if( family(i,3) == 0 && family(i,4) == 0 )
        founder(i) = 1;
    end
end


[r, c] = size(ref_alleles);
if( r ~= nind || c ~= 2 )
    error = 1;
    return;
end


[r, c] = size(alleles);
if( r ~= nind || c ~= 2 )
    error = 1;
    return;
end


% swap founder alleles to reduce recombination
% swap individuals with no genotyped ancestors?
new_alleles = alleles;
min_diff = nnz(ref_alleles-alleles);
for i = 1:nind
    p = new_alleles(i,1);
    m = new_alleles(i,2);     
    if( p ~= 0 && m ~= 0 )
        if( p == -i && m == i )
            % only swap founder alleles
            temp_alleles = new_alleles;
            for j = 1:nind
                if( new_alleles(j,1) == p && new_alleles(j,2) == m )
                    continue;
                end
                if( new_alleles(j,1) == m && new_alleles(j,2) == p )
                    continue;
                end
                temp_alleles(j, new_alleles(j,1:2) == p) = m;
                temp_alleles(j, new_alleles(j,1:2) == m) = p;
            end
            % swap_alleles is not exported, just for measuring purposes
%             [swap_alleles, diff] = realign(ref_alleles, temp_alleles);
%             if( isempty(swap_alleles) )
%                 error = 1;
%                 return;
%             end
            diff = nnz(ref_alleles - temp_alleles);
            if( diff < min_diff )
                min_diff = diff;
                new_alleles = temp_alleles;
            end
        end
    end
end


num_rec = min_diff;

end

















