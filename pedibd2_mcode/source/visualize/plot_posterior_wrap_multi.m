function plot_posterior_wrap_multi(input, output, parameters, f, list1, list2, limit)


count1 = 0;
count2 = 0;
temp = input.family_range{f}.structure(input.family_range{f}.family_range,8);   
for i = 1:length(list1)
    id1 = list1(i);
    in1 = find(temp == id1);
    if( length(in1) ~= 1 )
        continue;
    end
    count1 = count1 + 1;
    count2 = 0;
    for j = 1:length(list2)       
        id2 = list2(j);             
        in2 = find(temp == id2);
        if( length(in2) ~= 1 || in1 == in2 )
            continue;
        end
        count2 = count2 + 1;
        list(count1, count2, 1:2) = [id1, id2];
    end
end

figure;
count = 0;
for i = 1:count1
    for j = 1:count2
        count = count + 1;
        if( count1 == 1 )
            subplot(count2, count1, count);
        else
            subplot(count1, count2, count);
        end
        plot_posterior_wrap(input, output, parameters, f, list(i,j,1), list(i,j,2), limit);   
    end
end

end

