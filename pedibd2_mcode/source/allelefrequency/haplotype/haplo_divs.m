function [ result ] = haplo_divs( data, list )

d = data(:,7:end);

inv = 5;

count(1:floor(length(d(1,:))/(2*inv)),1:2^inv) = 0;

for i = 1:length(list)
    for j = 1:floor(length(d(list(i),:))/(inv*2))
        code = d(list(i),inv*2*(j-1)+1:inv*2*(j-1)+inv*2);
        de = 1;
        for k = 1:inv
            if( code(2*k-1) == 2 )
                de = de + 2^(k-1);
            end
        end
        count(j,de) = count(j,de) + 1;
        de = 1;
        for k = 1:inv
             if( code(2*k) == 2 )
                de = de + 2^(k-1);
            end           
        end
        count(j,de) = count(j,de) + 1;
    end
end

result.count = count;
[result.sort, temp] = sort(count, 2, 'descend');

diff(1:floor(length(d(1,:))/(2*inv))) = 0;

for i = 1:length(temp(:,1))

    code1 = temp(i,1);
    code2 = temp(i,2);
    for k = 1:inv
        if( bitget(code1,k) ~= bitget(code2,k) )
            diff(i) = diff(i) + 1;
        end
    end

end

result.diff = diff;


count = 0;

for i = 1:length(list)
    count = count + 1;
    code(count,1:length(d(1,:))/2) = d(list(i),1:2:end);
    count = count + 1;
    code(count,1:length(d(1,:))/2) = d(list(i),2:2:end);
end

inv2 = 200;
diversity(1:count*(count-1)/2, 1:floor(1+length(code(1,:))/inv2)) = 0; 

for i = 1:count
    for j = 1:i-1
        loc = (i-1)*(i-2)/2 + j;
        for k = 1:length(code(1,:))
            if( code(i,k) ~= code(j,k) )
                diversity(loc, 1+floor(k/inv2)) = diversity(loc, 1+floor(k/inv2)) + 1;
            end
        end
    end
end

result.div = max(diversity);

end

