function plot_linkage( allOUTPUT, marker_list )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
figure;
base = 0;
for chr = 1:length(allOUTPUT)
   
    
    sampled_markerlist = [];
    sampled_markerlist(:,2) = marker_list((marker_list(:,1) == chr), 2);
    sampled_markerlist(:,1) = 1:length(sampled_markerlist(:,1));
    
    %subplot(1,22,chr);
    all_score = allOUTPUT{chr}.combined_score;
    intervals = all_score.all_intervals;
    pvalue = all_score.all_p_value;
    inconsistency = allOUTPUT{chr}.consistency{1} < max(allOUTPUT{chr}.consistency{1});
    pscore = -1;
    for j = 1:size(intervals,1)
        p = pvalue(j);
        score = -log10(p + 1e-6); % number of permutation is 1e5
        bit = inconsistency(intervals(j,1):intervals(j,2));
        if( nnz(bit)/length(bit) >= 0.2 )
            continue;
        end
        line( base + sampled_markerlist(intervals(j,1:2),2), [score, score], 'LineWidth', 2, 'Color', [mod(chr-1,2),0,mod(chr,2)]);
        if( pscore > 0 )
            line( base + sampled_markerlist([intervals(j-1,2),intervals(j,1)],2), [pscore, score], 'LineWidth', 2, 'Color', [mod(chr-1,2),0,mod(chr,2)]);
        end
        pscore = score; % previous score

    end
    text(base + mean(sampled_markerlist(:,2)), -0.2,  [num2str(chr)], 'Color', [mod(chr-1,2),0,mod(chr,2)], 'FontSize', 18);
    base = base + max(sampled_markerlist(end,2));
    
end

xlabel({'';'Chromosome position'}, 'FontSize', 18);
ylabel('-log10(p value)', 'FontSize', 18);
title('Non-parametric linkage score', 'FontSize', 18);
set(gca, 'xtick',[]);
ylim([0,7]);

end

