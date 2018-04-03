%%%%%%%%%%allele frequency%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%

parameters.allele_freq = allele_frq(all_data, parameters.member_list, parameters.markerlist);


%%

x = parameters.markerlist(:,1);
y = parameters.allele_freq.all_af(:,1:3);
% y = simulation.allele_freq.all_af(:,1:3);

for i = 1:length(y(:,1))
    if( y(i,1) < y(i,2) )
        temp = y(i,1);
        y(i,1) = y(i,2);
        y(i,2) = temp;
    end
end

y(:,1) = smooth(y(:,1), 20);
y(:,2) = smooth(y(:,2), 20);
y(:,3) = smooth(y(:,3), 20);

plot(x,y);

clear i x y;

%%





