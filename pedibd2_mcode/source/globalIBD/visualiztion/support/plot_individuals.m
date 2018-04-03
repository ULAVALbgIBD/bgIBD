function plot_individuals(x_pos, y_pos, markersize, family_display, show_legend)


% pixels in area
markerarea = markersize^2;

legend_list = [];

group1 = family_display.group1;
group2 = family_display.group2;
group3 = family_display.group3;
group4 = family_display.group4;
group5 = family_display.group5;
group6 = family_display.group6;
scatter(x_pos(group1),y_pos(group1), markerarea*0.9, 'black', 's', 'MarkerFaceColor', 'k');
scatter(x_pos(group2),y_pos(group2), markerarea, 'black', 'o', 'MarkerFaceColor', 'k');
scatter(x_pos(group3),y_pos(group3), markerarea*0.9, 'black', 's', 'MarkerFaceColor', 'w');
scatter(x_pos(group4),y_pos(group4), markerarea, 'black', 'o', 'MarkerFaceColor', 'w');
light_grey = [0.8,0.8,0.8];
scatter(x_pos(group5),y_pos(group5), markerarea*0.9, 'black', 's', 'MarkerFaceColor', light_grey);
scatter(x_pos(group6),y_pos(group6), markerarea, 'black', 'o', 'MarkerFaceColor', light_grey); 

count = 0;
% group is a logical index length of all genotyped individuals
if( any(group1) )
    count = count + 1;
    legend_list{count} = 'affected male';
end

if( any(group2) )
    count = count + 1;
    legend_list{count} = 'affected female';
end

if( any(group3) )
    count = count + 1;
    legend_list{count} = 'non-affected male';
end

if( any(group4) )
    count = count + 1;
    legend_list{count} = 'non-affected female';
end

if( any(group5) )
    count = count + 1;
    legend_list{count} = 'non-status male';
end

if( any(group6) )
    count = count + 1;
    legend_list{count} = 'non-status female';
end



if( show_legend && ~isempty(legend_list) )
    h = legend(legend_list);
    set(h, 'Location', 'Northwest');
    set(h, 'EdgeColor', 'blue');
end




end