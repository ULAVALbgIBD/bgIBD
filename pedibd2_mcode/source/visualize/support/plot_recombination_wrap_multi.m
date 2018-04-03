function plot_recombination_wrap_multi(input, output, parameters, f, list1, list2, limit)


count = 0;
for i = 1:length(list1)
    for j = 1:length(list2)
        
        id1 = list1(i);
        id2 = list2(j);
        temp = input.family_range{f}.structure(input.family_range{f}.family_range,8);
        in1 = find(temp == id1);
        in2 = find(temp == id2);
        if( length(in1) ~= 1 || length(in2) ~= 1 || id1 == id2)
            continue;
        end
        count = count + 1;
        list(count,1:2) = [id1, id2];
    end
end

scrsz = get(0,'ScreenSize');
figure('OuterPosition',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2]);

for k = 1:count;
    
    subplot(count,1,k);
    plot_recombination_wrap(input, output, parameters, f, list(k,1), list(k,2), limit)

end

end