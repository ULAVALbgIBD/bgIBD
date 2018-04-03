
function [new_alleles, diff] = realign(ref_alleles, alleles)

%swap paternal and maternal alleles to reduce difference

new_alleles = [];
diff = 0;

[r1, c1] = size(ref_alleles);
if( r1 <=0 || c1 ~= 2 )
    return;
end

[r2, c2] = size( alleles );
if( r2 ~= r1 || c2 ~= 2 )
    return;
end

diff = 0;
new_alleles = alleles;
for i = 1:r1
    d1 = nnz(ref_alleles(i,1:2) - alleles(i,1:2));
    d2 = nnz(ref_alleles(i,1:2) - alleles(i,2:-1:1));
    if( d1 <= d2 )
        diff = diff + d1;
        new_alleles(i,1:2) = alleles(i,1:2);
    else
        diff = diff + d2;
        new_alleles(i,1:2) = alleles(i,2:-1:1);
    end
end


end

