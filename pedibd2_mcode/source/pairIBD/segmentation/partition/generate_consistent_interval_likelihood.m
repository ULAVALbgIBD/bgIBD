function [ stat_interval ] = generate_consistent_interval_likelihood( likelihood )



interval = 0;
count = 0;

stat_interval = [];

former_r(1:length(likelihood)) = 0;

threshold = 1;

for locus = 1:1:length(likelihood{1}(:,1))
    
    same = 1;

    for i = 1:length(likelihood) % i index number of relative pairs
        [temp1, temp2] = max(likelihood{i}(locus,:));
        if( former_r(i) ~= temp2 )         
            same = 0;
        end
        former_r(i) = temp2;
    end

    if( same == 1 )
        interval = interval + 1;
    else
        if( interval > threshold )
            count = count + 1;
            stat_interval(count,1:3) = [locus-interval-1, locus-1, round(locus-interval/2-1)];
        end
        interval = 0;
    end
end

if( interval > threshold )
    count = count + 1;
    stat_interval(count,1:3) = [locus-interval, locus, round(locus-interval/2-1)];
end

end

