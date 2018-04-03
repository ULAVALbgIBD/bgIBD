function accuracy = verify_imputation(haplotype, data, family, id)

accuracy = 0;

[nind, fields] = size(family);

if( nind <= 0 || fields < 12 )
    return;
end

if( id <= 0 )
    return;
end

outerid = family(1:nind, 8);
innerid = find(outerid == id);
if( length(innerid) ~= 1 )
    return;
end

family_id = family(1,12);

[r c] = size(data);
if( r <= 0 || c < 6 || mod(c - 6, 2) ~= 0 )
    return;
end

index = find(data(1:r,2) == id & data(1:r,1) == family_id);
if( length(index) ~= 1 )
    return;
end
genotype = data(index, 7:c);

nmarkers = (c-6)/2;

[r c] = size(haplotype);
if( r ~= nind || c ~= 2 * nmarkers )
    return;
end

count = 0;
hap = 0;
bench = 0;
valid(1:nmarkers) = 0;
for i = 1:nmarkers
    if( any(haplotype(innerid,i*2-1:i*2) == 0) )
        continue;
    end
    hap = hap + 1;
    if( any(genotype(i*2-1:i*2) == 0) )
        continue;
    end
    bench = bench + 1;
    valid(i) = 1;
    if( ibs(haplotype(innerid,i*2-1:i*2), genotype(i*2-1:i*2)) ~= 2 )
        count = count + 1;
        valid(i) = 2;
    end
end

accuracy = valid;

end

