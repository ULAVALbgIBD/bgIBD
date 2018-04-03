function p_value = categorical_null_ibs(emission, observation, ibd)

global temp1;
global temp2;

len1 = length(emission);
len2 = length(observation);
len3 = length(ibd);
p_value = 1;

if( len1 ~= len2 || len2 ~= len3 )
    disp('error in check data consistency');
    return;
end

if( len1 <= 0 )
    disp('no genotype data for evaluating relationship');
    return;
end


len = len1;

sum1 = 0;
sum2 = 0;
sum3 = 0;
sum4 = 0;
count = 0;
for i = 1:len
    if( observation(i) == 4 )
        continue;
        % do not include missing observations
    else
        count = count + 1;
    end
    temp = sum(emission{i}(ibd(i),1:3));
    em = emission{i}(ibd(i),1:3)./temp;
    sum1 = sum1 + normalize_stat(em, observation(i));
    sum3 = sum3 + observation(i)-1;
    sum4 = sum4 + em*[1;2;3]-1;
end
stat = sum1/(sqrt(count));

% test is assuming all observations are independent
% therefore not considering background ibd, or LD

% stat is approximately normal distribution N(0,1)

[h,p_value] = ztest(stat, 0, 1);

temp2 = temp2 + 1;
temp1(temp2,1) = sum3/count;
temp1(temp2,2) = sum4/count;
temp1(temp2,3) = stat;
temp1(temp2,4) = p_value;

end