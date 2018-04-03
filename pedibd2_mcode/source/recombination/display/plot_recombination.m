
% takes input of 1 column only

function plot_recombination(rec)

scrsz = get(0,'ScreenSize');
figure('OuterPosition',[1 scrsz(4)/3 scrsz(3)*3/4 scrsz(4)/3]);


r = sort(rec);

hold on;
i = 0;
y = [];
y(1:length(r)) = i;
scatter(r, y, 'x');

ylim([i - 0.2, i + 0.2]);


end





