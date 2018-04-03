function [align, false] = recombination_align(res, ref, sampled_markerlist)

map = sampled_markerlist(:,2);

a = length(res)+1;
[b, b1] = size(ref);
b = b + 1;
mismatch = 1000;
score(1:a, 1:b) = 100000;
score(1,1) = 0;
pre(1:a,1:b,1:2) = 0;
for i = 1:a
    for j = 1:b
        if( i > 1 )
            if( score(i,j) > score(i-1,j) + mismatch )
                score(i,j) = score(i-1,j) + mismatch;
                pre(i,j,1:2) = [i-1,j];
            end
        end
        if( j > 1 )
            if( score(i,j) > score(i,j-1) + mismatch )
                score(i,j) = score(i,j-1) + mismatch;
                pre(i,j,1:2) = [i,j-1];
            end
        end
        if( i > 1 && j > 1 )
            temp_i = i - 1;
            temp_j = j - 1;
            temp_lower = ref(temp_j,1);
            temp_higher = ref(temp_j,2);
            if( temp_lower > temp_higher )
                temp_lower = ref(temp_j,2);
                temp_higher = ref(temp_j,1);
                ref(temp_j,2) = temp_higher;
                ref(temp_j,1) = temp_lower;
            end
            if( res(temp_i) <= ref(temp_j,2) && res(temp_i) >= ref(temp_j,1) )
                difference = 0;
            else
                difference = min(abs(res(temp_i)-ref(temp_j,1)), abs(res(temp_i)-ref(temp_j,2)));
            end
            if( score(i,j) > score(i-1,j-1) + difference )
                score(i,j) = score(i-1,j-1) + difference;
                pre(i,j,1:2) = [i-1,j-1];
            end
        end
    end
end

align(1:a+b-2, 1:3) = 0;
align(1:a+b-2, 4) = -1;

i = a;
j = b;
count = 0;
fp = 0;
fn = 0;
while ( i ~= 1 || j ~= 1 )
    temp_i = pre(i,j, 1);
    temp_j = pre(i,j, 2);
    if( temp_i < i && temp_j >= j )
        count = count + 1;
        align(count, 3) = res(temp_i);
        fp = fp + 1;
    end
    if( temp_j < j && temp_i >= i )
        count = count + 1;
        align(count, 1:2) = ref(temp_j,1:2);
        fn = fn + 1;
    end
    if( temp_i < i && temp_j < j )
        count = count + 1;
        align(count, 1:2) = ref(temp_j,1:2);
        align(count, 3) = res(temp_i);
        if( res(temp_i) <= ref(temp_j,2) && res(temp_i) >= ref(temp_j,1) )
            difference = 0;
        else
            difference = min(abs(res(temp_i)-ref(temp_j,1)), abs(res(temp_i)-ref(temp_j,2)));
            difference = min(abs(map(res(temp_i))-map(ref(temp_j,1))), abs(map(res(temp_i))-map(ref(temp_j,2))));
        end        
        align(count, 4) = difference;
    end
    i = temp_i;
    j = temp_j;
end

align(count+1:end,:) = [];

false = [fp,a-1,fn,b-1];

end
