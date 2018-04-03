function plot_title()


text(-0.1, 0, 'PedIBD', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right', 'FontUnits', 'normalized', 'FontSize', 0.5, 'FontName', 'Comic Sans MS', 'FontWeight', 'bold');
text(0, 0, datestr(now), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left', 'FontUnits', 'normalized', 'FontSize', 0.5, 'FontName', 'Arial Rounded MT');

xlim([-1,1]);
ylim([-1,1]);

axis off;

end