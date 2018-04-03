
function analyse_result(family, all_inheritance, parameters, family_genotype, input_pair, limit)



parameters.typing_error = 0.001;

[viterbi, posterior, ~, oblist, emission, ~, ~] = posterior_decode_1pair(family, all_inheritance, parameters, family_genotype, input_pair);


sampled_markerlist = parameters.sampled_markerlist;
sbps = limit(1);
ebps = limit(2);
coor = sampled_markerlist(:,1);

iv = 0;
[d1, d2, d3] = size( emission );
nibs = d3;

[nloc c] = size(sampled_markerlist);
if( c ~= 2 || nloc <= 0 )
    return;
end

%1 ibd0
%2 ibd1, paternal
%8 ibd1, maternal
%10 ibd2

scrsz = get(0,'ScreenSize');
figure('OuterPosition',[1 scrsz(4)/3 scrsz(3) scrsz(4)/3]);
hold on;

nloc = length(coor);

likelihood = zeros(nloc, 4); % hidden states 4: [1, 2, 8, 10]
for i = 1:4 % emission states 4
    likelihood(oblist == i, 1:4) = emission(oblist==i, [1, 2, 8, 10], i);
end

% temp0(1:nloc) = emission(1:nloc, 1, 1:nibs); %ibd0
% %for siblings, paternal and maternal sharing emission is different
% temp1(1:nloc) = emission(1:nloc, 2, 1:nibs); %ibd1, paternal
% temp2(1:nloc) = emission(1:nloc, 8, 1:nibs); %ibd1, maternal
% temp3(1:nloc) = emission(1:nloc, 10, 1:nibs);%ibd2


likelihood = log10(likelihood);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot observation
subplot(4,1,1);
temp0 = oblist == 1;
temp1 = oblist == 2;
temp2 = oblist == 3 | oblist == 4;

if( ebps - sbps > 500 )

    plot(sampled_markerlist(sbps:ebps,1), random_ibs(oblist(sbps:ebps), []), '.');

else

    id1 = input_pair(1);
    id2 = input_pair(2);
    fa1 = family(id1,3);
    ma1 = family(id1,4);
    fa2 = family(id2,3);
    ma2 = family(id2,4);
    nmarkers = length(coor);
    geno(1:4, 1:nmarkers, 1:2) = 0;
    geno(1, 1:nmarkers, 1:2) = family_genotype(1:nmarkers, id1, 1:2);
    geno(2, 1:nmarkers, 1:2) = family_genotype(1:nmarkers, id2, 1:2);
    if( fa1 == fa2 && fa1 ~= 0 && ma1 == ma2 && ma1 ~= 0 )
        geno(3, 1:nmarkers, 1:2) = family_genotype(1:nmarkers, fa1, 1:2);
        geno(4, 1:nmarkers, 1:2) = family_genotype(1:nmarkers, ma1, 1:2);
    end
    for i = sbps:ebps
        for j = 1:4
            text(coor(i), j-0.1, num2str(geno(j,i,1)));
            text(coor(i), j+0.1, num2str(geno(j,i,2)));
        end
    end    
    ylim([0,5]);
    xlim([sbps,ebps]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot emission

subplot(4,1,2);
hold on;

plot(sampled_markerlist(sbps:ebps,1), smooth(likelihood(sbps:ebps, 1), iv), 'r');%ibd0
plot(sampled_markerlist(sbps:ebps,1), smooth(likelihood(sbps:ebps, 2), iv), 'g');%ibdL
% plot(sampled_markerlist(sbps:ebps,1), smooth(likelihood(sbps:ebps, 3), iv), 'y');%ibdR
plot(sampled_markerlist(sbps:ebps,1), smooth(likelihood(sbps:ebps, 4), iv), 'b');%ibd2

% ylim([-2,-0.5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot posterior

subplot(4,1,3);
hold on;
dec = posterior;
plot(sampled_markerlist(sbps:ebps,1), dec(sbps:ebps,1), 'r');
plot(sampled_markerlist(sbps:ebps,1), dec(sbps:ebps,2), 'g');
plot(sampled_markerlist(sbps:ebps,1), dec(sbps:ebps,3), 'b');


 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot viterbi

subplot(4,1,4);
hold on;
plot(sampled_markerlist(sbps:ebps,1), viterbi(sbps:ebps,1), 'b', 'linewidth', 3);
plot(sampled_markerlist(sbps:ebps,1), viterbi(sbps:ebps,2), 'g', 'linewidth', 3);
ylim([0.9,3.1]);


end





