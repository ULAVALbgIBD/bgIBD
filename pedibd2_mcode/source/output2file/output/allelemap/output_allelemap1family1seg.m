
function error = output_allelemap1family1seg(filename, alleles_all, family, chr, basepair, p_value)

error = 0;

[family_display error] = generate_family_display(alleles_all, family);

if( error ~= 0 )
    error = 1;
    disp('error in global IBD');
    return;
end

if( isempty(family_display) )
    return;
end

% draw allele map
scale = 1;
hFig = figure('Position',[100 100 scale * 1000  scale * 700], 'Visible', 'off');
markersize = view_mendelian_map(family_display, 0, scale * 800);

% draw allele assignment
hax = axes('Position', [0.16, 0.15, 0.72, 0.3]);
plot_table_vertical(family_display, markersize);


% plot allele information
hax = axes('Position', [0.75, 0.78, 0.15, 0.14]);
plot_segment_info(chr, basepair, p_value);

hax = axes('Position', [0.6, 0.07, 0.3, 0.03]);
plot_title();


I = getframe(gcf);
imwrite(I.cdata, filename);

set(hFig, 'Visible', 'off');
delete(findobj(0, 'type', 'figure'));


end


