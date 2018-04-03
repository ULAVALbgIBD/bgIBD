
function plot_ibd(haplotype1, haplotype2, coor, id1, id2, limit)


count = 0;
ks = limit(1);
js = limit(2);


hold on;
text((coor(ks)+coor(js))/2,2.5,[num2str(id1),'-',num2str(id2)]);

mark(1:length(coor)) = -0.1;
mark = random_ibs(mark);


for i = 1:2
    subplot(2,1,i);
    temp1 = haplotype1(:,1) == haplotype2(:,i);
    temp2 = haplotype1(:,2) == haplotype2(:,i);
    hold on;
    plot(coor(temp1(ks:js)), mark(temp1(ks:js)), '.', 'color', 'r', 'MarkerSize', 2);
    plot(coor(temp2(ks:js)), mark(temp2(ks:js)), '.', 'color', 'b', 'MarkerSize', 2);
    ylim([-1,1]);
end




end

