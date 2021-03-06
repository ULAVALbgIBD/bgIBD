% double chain

function e = id_condensed_generate_emission_genotype_nonmissing(freq, is)

num = length(is(1,:));
num_es = 2^(num-1);
num_is = length(is(:,1));

e(1:num_is,1:num_es) = 0;


p1 = freq(1);
p2 = freq(2);

for i = 1:num_is
    id = is(i,:);
    code(1:num) = 1;
    while(1)     
        shift = 1;
        prob = 1;
        es = 1;
        for j = 1:num
            if( id(j) ~= j )
                code(j) = code(id(j));
            else
                if( shift == 1 )
                    if( code(j) == 1 )
                        code(j) = 2;
                        shift = 0;
                    else
                        code(j) = 1;
                    end
                end
                prob = prob * freq(code(j));
            end
            es = es + (code(j) - 1) * 2^(num-j);
        end
        e(i,es) = prob;
        if( shift == 1 )
            break;
        end
    end
end


end


