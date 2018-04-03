function plot_table_vertical(family_display, markersize)


n_ind = family_display.num_ind;
if( n_ind <= 50 )
    height = 10;
else
    if( n_ind <= 100 )
        height = 20;
    else
        n_ind = 100;
        % remaining alleles skipped
    end
end

cols = ceil(n_ind/height);

allele_string = family_display.allele_string;
id_string = family_display.id_string;

width = 4;
dist = 2;

fontsize = min(0.8/(height),0.8/(4*cols));
if( height < 15 )
    fontsize = min(0.8/(15),0.8/(4*cols));
end

for i = 1:n_ind

    order = floor((i-1)/height);
    y = dist * mod(i-1,height)+1;
    text(order*width+1, -y, id_string{i}, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized', 'FontSize', fontsize);
    text(order*width+2, -y, allele_string{i,1}, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized', 'FontSize', fontsize);
    text(order*width+3, -y, allele_string{i,2}, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized', 'FontSize', fontsize);
    text(order*width+2.5, -y, 'x', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized', 'FontSize', fontsize);
    x_pos(i) = order*width + 0.5;
    y_pos(i) = -y;
end
hold on;

%%%% mark each individual

plot_individuals(x_pos, y_pos, markersize, family_display, 0); 

% bound box
leftmost = 0;
if( cols < 4 )
    rightmost = width * 4 + 0.5;
else
    rightmost = width * cols + 0.5;
end
uppermost = 0.5;
lowermost = -dist * height - 1;


line([leftmost,rightmost], [uppermost, uppermost], 'linewidth', 1);
line([leftmost,rightmost], [lowermost, lowermost], 'linewidth', 1);


xlim([leftmost, rightmost]);
ylim([lowermost,uppermost]);

axis off;

end











