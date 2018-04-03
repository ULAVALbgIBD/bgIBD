
function plot_haplotype(haplotype1, haplotype2, coor, id1, id2, limit)


count = 0;
ks = limit(1);
js = limit(2);


hold on;
text((coor(ks)+coor(js))/2,2.5,[num2str(id1),'-',num2str(id2)]);


for i = 1:2
    for j = 1:2
        temp = haplotype1(:,i) == haplotype2(:,j);
        count = count + 1;
        plot(coor(ks:js), count*2 + random_ibs(temp(ks:js)), '.');
%         plot(coor(ks:js), count*2 + smooth(temp(ks:js), 100), '.');
%         plot(coor(ks:js), count*2 + (temp(ks:js)), '.');
    end
end






end