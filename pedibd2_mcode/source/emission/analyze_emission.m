function analyze_emission(data, ref, code, cate)





temp = find(code == cate);


for i = 1:length(temp)
    temp1(i,1:4) = data(temp(i),1:4)./sum(data(temp(i),1:4));
    temp1(i,3) = temp1(i,3)+temp1(i,4);
end

for i = 1:length(temp)
    temp2(i,1:4) = ref(temp(i),1:4)./sum(ref(temp(i),1:4));
    temp2(i,3) = temp2(i,3)+temp2(i,4);
end

figure;
plot(temp1(:,1:3));
hold on;
plot(temp2(:,1:3));

m(1:3) = mean(temp1(:,1:3),1);
for i = 1:3
    text(length(temp)/2, m(i), ['average: ',num2str(m(i))]);
end


count = 0;
for p = 0.25:0.01:0.3
    count = count + 1;
    q = 1-p;
    freq(count,1:2) = [p,q];
    ms = 0;
    te = 0.02;
    ep_merge_ibs(count,1:3,1:4) = dc_generate_emission_ibs(freq(count,:), ms, te);
    for i = 1:length(temp)
        temp3(i,1:3) = ep_merge_ibs(count,cate,1:3);
    end
    plot(temp3, '--', 'linewidth', (p-0.2)*30);
end




end