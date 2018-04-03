function plot_table(family_display)


n_ind = family_display.num_ind;
len = 20;
lines = ceil(n_ind/20);

allele_string = family_display.allele_string;
id_string = family_display.id_string;

for i = 1:n_ind

    order = floor((i-1)/len);
    x = mod(i-1,len)+1;
    text(x,-(1+order*4),id_string{i}, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontSize', 8);
    text(x,-(2+order*4),allele_string{i,1}, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontSize', 8);
    text(x,-(3+order*4),allele_string{i,2}, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontSize', 8);
    
end

for i = 1:lines+1
    order = i - 1;
    line([0.5,len+0.5], [-1.5-order*4,-1.5-order*4]);
end
line([0.5,len+0.5], [0,0], 'linewidth', 2);
line([0.5,len+0.5], [-lines*4,-lines*4], 'linewidth', 1);

xlim([0.5,len+0.5]);
ylim([-lines*4-0.5,0.5]);

axis off;

end