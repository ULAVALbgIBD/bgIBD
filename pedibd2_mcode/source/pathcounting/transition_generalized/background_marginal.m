function [ prob ] = background_marginal( is, m )


% ancestral group method k = 1/m
% m is the probability that given two alleles they are ibd


temp(1:length(is(:,1))) = 0;
prob(1:length(is(:,1))) = 1;

for i = 1:length(is(:,1))
    for j = 1:length(is(i,:))
        if( is(i,j) ~= j )
            temp(i) = temp(i) + 1;
        end
    end
end



for i = 1:length(temp)
    prob(i) = prob(i) * m^(temp(i));
    for j = 1:length(is(1,:))-temp(i)-1
        prob(i) = prob(i) * (1-j*m);
    end
end

end

