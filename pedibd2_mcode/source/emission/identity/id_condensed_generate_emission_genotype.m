% double chain

function ep = id_condensed_generate_emission_genotype(freq, is, ms, te)


e = id_condensed_generate_emission_genotype_nonmissing(freq, is);

num = length(is(1,:));
num_es = 2^num;
num_is = length(is(:,1));

% typing error emission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = [0.5,0.5];

er(1:num_es) = 0;
code(1:num) = 1;

while(1)
    shift = 1;
    prob = 1;
    es = 1;
    for i = 1:num
        if( shift == 1 )
            if( code(i) == 1 )
                code(i) = 2;
                shift = 0;
            else
                code(i) = 1;
            end
        end
        es = es + (code(i)-1) * 2^(num-i);
        prob = prob * p(code(i));
    end
    er(es) = prob;
    if( shift == 1 )
        break;
    end
end

e2(1:num_is,1:num_es) = 0;

for i = 1:num_is
    for j = 1:num_es
        e2(i,j) = (1-te) * e(i,j) + te * er(j);
    end
end

% add effect of missing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



e3(1:num_is,1:num_es) = 0;
% assume ms is value for genotype missing
nonmissing = ((1 - ms)^0.5)^num;


for i = 1:num_is
    for j = 1:num_es
        e3(i,j) = nonmissing * e2(i,j);
    end
end

e3(1:num_is,num_es+1) = 1 - nonmissing;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ep = e3;

end


