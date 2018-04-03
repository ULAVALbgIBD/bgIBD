function result = generate_genotype(ind1, ind2, emission_option)

result = [];


[r1, c1] = size(ind1);
[r2, c2] = size(ind2);

if( r1 ~= r2 || r1 <= 0 )
    disp('error in processing genotype');
    return;
end
if( c1 ~= c2 || c1 ~= 2 ) 
    disp('error in processing genotype');
    return;
end

nloc = r1;

static_code(1:3,1:3,1:3,1:3) = 0;
for i1 = 0:2
    for i2 = 0:2
        for j1 = 0:2
            for j2 = 0:2
                pair_code = code_pair([i1,i2,j1,j2], emission_option);
                if( isempty(pair_code) )
                    disp('error in genotype coding');
                    result = [];
                    return;
                end                
                static_code(i1+1,i2+1,j1+1,j2+1) = pair_code;
            end
        end
    end
end

result(1:nloc) = 0;
for i = 1:nloc
    pair_code = static_code(ind1(i,1)+1, ind1(i,2)+1, ind2(i,1)+1, ind2(i,2)+1);
    result(i) = pair_code;
end



end







