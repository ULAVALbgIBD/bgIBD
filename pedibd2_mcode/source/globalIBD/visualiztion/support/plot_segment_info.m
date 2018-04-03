

function plot_segment_info(chr, basepair, p_value)



row = -1;
text(1, row, ['CHR ', num2str(chr)], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized');
row = row - 1;
% reserve one row for region
region_row = row;
row = row - 1;
text(1, row, [insertCommas(basepair(1)), ' bp'], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized');
row = row - 0.8;
fromto_row = row;
row = row - 0.8;
text(1, row, [insertCommas(basepair(2)), ' bp'], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized');
row = row - 1;
p_value_row = row;
row = row - 1;
if( ~isempty(p_value) && p_value >= 0 && p_value <= 1 )    
    degree = -log10(p_value);
    if( degree > 2 )
        precision = '%.5f';
    else
        precision = '%.2f';
    end
    text(1, row, [num2str(p_value, precision)], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized');
else
    text(1, row, ['N/A'], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized');
end


% all text in current axes
h = findobj(gca,'Type','text');
set(h, 'FontSize', 1/(-row+2));

h = text(1, region_row, 'region', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized');
set(h, 'FontSize', 0.8/(-row+2));

h = text(1, fromto_row, '~', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized');
set(h, 'FontSize', 0.8/(-row+2));

h = text(1, p_value_row, 'p-value', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized');
set(h, 'FontSize', 0.8/(-row+2));

% bound box
leftmost = 0;
rightmost = 2;
uppermost = 0;
lowermost = row - 1;

line([leftmost,rightmost], [uppermost, uppermost]);
line([leftmost,rightmost], [lowermost, lowermost]);
line([leftmost,leftmost], [lowermost, uppermost]);
line([rightmost,rightmost], [lowermost, uppermost]);


xlim([leftmost,rightmost]);
ylim([lowermost,uppermost]);
axis off;

end