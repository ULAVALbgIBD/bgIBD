function error = plot_genotype(family_genotype, family, list, coor, limit)

error = 0;

nloc = length(coor);
if( nloc <= 0 )
    error = 1;
    return;
end

[nind, c] = size(family);
if( nind <= 0 || c < 12 )
    error = 1;
    return;
end

sbps = limit(1);
ebps = limit(2);

mid = mean(coor(limit));
[num c] = size(list);

bit(1:nloc) = false;
bit(sbps:ebps) = true;
bit_ind = family_genotype(1:num, 1:nloc, 1) ~= family_genotype(1:num, 1:nloc, 2);



temp = smooth(double(bit),10);
plot(coor(sbps:ebps,1), 10 * temp(sbps:ebps));

count = sum(bit_ind(1:num,1:nloc), 1);

bit = bit & (count >= num);
nnz(bit)

hold on;

for i = 1:num
    
    text(mid, i*4 - 2 + 1, num2str(family(list(i),8)), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized', 'FontSize', 0.05);
    id = list(i);
    
    oblist = family_genotype(id, 1:nloc, 1) ~= family_genotype(id, 1:nloc, 2);

    
    if( ebps - sbps > 500 )

        plot(coor(bit,1), i*4 + random_ibs((oblist(bit))), '.');

    else
        
        for j = sbps:ebps
            text(coor(j), i-0.1, num2str(family_genotype(id,j,1)));
            text(coor(j), i+0.1, num2str(family_genotype(id,j,2)));
        end    
        ylim([0,5]);
        xlim([sbps,ebps]);
    end

end


ylim([0,num*4+4]);


figure;
hist(count, max(count)+1);

















