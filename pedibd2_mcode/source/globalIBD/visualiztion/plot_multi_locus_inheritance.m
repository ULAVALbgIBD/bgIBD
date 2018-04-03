
% assigned_list specifies which individual to be included in the global
% inheritance

function plot_multi_locus_inheritance( file_head, pm, input_family_range, assignment, pedigree_all_missing )


chr = pm.chr;
stat_interval = assignment.intervals;
stat_assignment_all = assignment.alleles;
f = pedigree_all_missing(input_family_range(1),1);

sampled_markerlist = pm.sampled_markerlist;

assigned_list = 1:length(input_family_range);

ix = length(stat_interval(1,:));

newdir = sprintf('%s/family%d/chr%d', file_head, f, chr);

mkdir(newdir);

for l = 1:length(stat_interval(:,ix))

    
    scrsz = get(0,'ScreenSize');
    figure('OuterPosition',[1 0 scrsz(3) scrsz(4)]);
    clear scrsz;

    plot_one_locus_inheritance( assigned_list, input_family_range, pedigree_all_missing, stat_assignment_all{l} );
    
    title(['Descent graph: chr', num2str(chr), ': ', num2str(sampled_markerlist(stat_interval(l,1),2)), '-', num2str(sampled_markerlist(stat_interval(l,2),2)), 'bps score: ', num2str((stat_interval(l,ix)))], 'fontsize', 6, 'color', 'g');
    
    fOut = sprintf('%s/family%d/chr%d/chr%02d_%d_%d_%d', file_head, f, chr, chr, l, sampled_markerlist(stat_interval(l,1),2), sampled_markerlist(stat_interval(l,2),2));
    fOut = sprintf('%s/family%d/chr%d/chr%02d_%d', file_head, f, chr, chr, l);

    
    saveas(gcf, fOut, 'jpg')
    close(gcf);
    
    
end



end

