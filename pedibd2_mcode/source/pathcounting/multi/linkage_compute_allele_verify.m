

function compare = linkage_compute_allele_verify(inheritance, family, ii1, ii2)

global debug_mode;
% this linkage computing is not correct, modify this to correct errors
lk_result = linkage_compute_allele(family, ii1, ii2);

fast_index = inheritance.fast_index;
content = inheritance.content;


compare = lk_result;
compare(:,:) = 0;

slots = 2;
alleles = [ii1*2-1, ii1*2, ii2*2-1, ii2*2];


% extract directly from inheritance is correct.
count = 0;
for p = 1:4
    for q = p+1:4
        i = alleles(p);
        j = alleles(q);
        count = count + 1;
        if( fast_index(i,j) ~= 0 )
            temp = content{fast_index(i,j)};
            for k = 1:length(temp(:,1))
                len = sum(temp(k,slots+1+1:end));
                if( temp(k,2) == temp(k,3) && len ~= 0 )
                    compare(count, len) = temp(k,1);
                end
            end
        end
    end
end

% not all families can be processed using this function
if( debug_mode == 1 )
    [d1 d2] = size(compare);
    [d3 d4] = size(lk_result);
    if( d1 ~= d3 || d2 ~= d4 || nnz(xor(compare, lk_result)) ~= 0 )
        disp(ii1);
        disp(ii2);
        disp(lk_result);
        disp(compare);
        disp('discrepancy in calculating inheritance path');
    end
end

end



