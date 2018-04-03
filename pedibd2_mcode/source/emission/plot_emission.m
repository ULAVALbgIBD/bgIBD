function plot_emission()

figure;

eo = generate_emission_option(1);

count = 0;
for p = 0:0.01:1
    count = count + 1;
    q = 1-p;
    freq(count,1:2) = [p,q];
    ms = 0;
    te = 0;
    ep_merge_genotype(count,1:3,1:10) = dc_generate_emission_genotype(freq(count,:), ms, te, eo);
end

g11 = 1; d{1} = '11';
g12 = 2; d{2} = '12';
g22 = 3; d{3} = '22';
g0 = 10;

for i = 1:3
    for j = 1:3
        index = (i-1)*3 + j;
        subplot(4,3,index);
        plot(freq(:,1), ep_merge_genotype(:,1:3,index));
        text(freq(20,1), ep_merge_genotype(20, 1, index), 'ibd0');
        text(freq(20,1), ep_merge_genotype(20, 2, index), 'ibd1');
        text(freq(20,1), ep_merge_genotype(20, 3, index), 'ibd2');

        
        title([d{i}, '/', d{j}]);
    end
end


subplot(4,3,10);
temp = sum(ep_merge_genotype(:,1:3,[3,7]),3);
plot(freq(:,1), temp);
title('ibs0');
text(freq(20,1), temp(20, 1), 'ibd0');
text(freq(20,1), temp(20, 2), 'ibd1');
text(freq(20,1), temp(20, 3), 'ibd2');

subplot(4,3,11);
temp = sum(ep_merge_genotype(:,1:3,[2,4,6,8]),3);
plot(freq(:,1), temp);
title('ibs1');
text(freq(20,1), temp(20, 1), 'ibd0');
text(freq(20,1), temp(20, 2), 'ibd1');
text(freq(20,1), temp(20, 3), 'ibd2');

subplot(4,3,12);
temp = sum(ep_merge_genotype(:,1:3,[1,5,9]),3);
plot(freq(:,1), temp);
title('ibs2');
text(freq(20,1), temp(20, 1), 'ibd0');
text(freq(20,1), temp(20, 2), 'ibd1');
text(freq(20,1), temp(20, 3), 'ibd2');

end


