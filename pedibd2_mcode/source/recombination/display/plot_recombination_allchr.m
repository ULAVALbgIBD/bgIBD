function error = plot_recombination_allchr(output, input, marker_list_allchr)

error = 0;

fid = 1;

[nloc c] = size( marker_list_allchr );
if( nloc <= 0 || c ~= 2 )
    error = 1;
    disp('error in marker list');
    return;
end

figure;
hold on;

pcount = 0;
mcount = 0;

for chr = 1:22
    
    rec_list = output{chr}.merged_assignment{fid}.recombination.list;
    [nr c] = size(rec_list);
    if( nr <= 0 )
        continue;
    end
    if( c ~= 4 )
        error = 1;
        disp('error in recombination format');
        return;
    end
    
    family = input{chr}.family_range{fid}.structure;
    
    [nind c] = size(family);
    if( nind <= 0 || c ~= 12 )
        error = 1;
        disp('error in family structure');
        return;
    end
    sex_code = family(1:nind, 5);
    if( any( sex_code ~= 1 & sex_code ~= 2 ) )
        error = 1;
        disp('error in sex code');
        return;
    end
    
    bps = marker_list_allchr(marker_list_allchr(:,1) == chr, 2);
    y = zeros(length(bps),1);
    sexbit = zeros(length(bps),1);
    
    bit = rec_list(1:nr,3) ~= 0 &  rec_list(1:nr,4) == 0;
    sexbit(bit) = sex_code(rec_list(bit,3));
    pbit = sexbit == 1;
    mbit = sexbit == 2;
    pcount = pcount + nnz(pbit);
    mcount = mcount + nnz(mbit);
    
    loc = round(mean(rec_list(1:nr,1:2),2));
    
    if( chr <= 11 )
    
        scatter(bps(loc(pbit)), y(pbit) + chr - 0.1, 40, sexbit(pbit), 'x', 'LineWidth', 2);
        scatter(bps(loc(mbit)), y(mbit) + chr + 0.1, 40, sexbit(mbit), 'x', 'LineWidth', 2);
        
        text(-10^7, chr, ['chr ', num2str(chr)], 'HorizontalAlignment', 'right');
    
    else
        
        scatter(3.2 * 10^8 - bps(loc(pbit)), y(pbit) + 23 - chr - 0.1, 40, sexbit(pbit), 'x', 'LineWidth', 2);
        scatter(3.2 * 10^8 - bps(loc(mbit)), y(mbit) + 23 - chr + 0.1, 40, sexbit(mbit), 'x', 'LineWidth', 2);        

        text(3.2 * 10^8 + 10^7,  23 - chr, ['chr ', num2str(chr)], 'HorizontalAlignment', 'left');
    
    end
    
end

set(gca,'YDir','reverse');
axis off;

disp(['paternal recombinations ', num2str(pcount)]); 
disp(['maternal recombinations ', num2str(mcount)]);

end











