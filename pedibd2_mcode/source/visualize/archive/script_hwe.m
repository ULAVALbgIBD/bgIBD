%%
%%Hardy-Weinberg Equilibrium

odd(1:32550) = 0;

for i = 1:32550
    %all_hetero{i} = get_hetero_locus(all_data, i, member_list);
    odd(i) = (all_af{i}(1)*all_af{i}(1)+all_af{i}(2)*all_af{i}(2));
end

%%
a = 10980;
b = 11020;
cdfplot(odd(a:b));
hold on;
cdfplot(all_hr(a:b));


%%

temp = [];
temp(1:32550, 1:3) = 0;

for i = 1:32550
    temp1 = get_hetero_locus(data, i, 2);
    temp(i,1:3) = temp1;
    %all_af{i} = get_afreq(all_freq{i}, de_freq);
end
