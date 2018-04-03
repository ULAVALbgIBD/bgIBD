function [weight error] = founder_source1family(family)


error = 0;
[nind, ~] = size(family);

% nearest genotyped ancestor or trace until founder
% half-founder (one parent missing) not counted as founder
ancestor = zeros(nind, nind);
degree = zeros(nind, nind);


genotyped = family(1:nind,7) == 1;
founder = family(1:nind,3) == 0 & family(1:nind,4) == 0;


% marker the closest genotyped ancestor or founder
for i = 1:nind
    id = family(i,2);
    father = family(i,3);
    mother = family(i,4);
    if( father ~= 0 )       
        if( genotyped(father) || founder(father) )
            ancestor(id,father) = 1;
            degree(id,father) = 1 * 0.5;
        else
            ancestor(id,:) = ancestor(father,:) + ancestor(id,:);
            degree(id,:) = degree(father,:) * 0.5 + degree(id,:);
        end
    end
    if( mother ~= 0 )
        if( genotyped(mother) || founder(mother) )
            ancestor(id,mother) = 1;
            degree(id,mother) = 1 * 0.5;
        else
            ancestor(id,:) = ancestor(mother,:) + ancestor(id,:);
            degree(id,:) = degree(mother,:) * 0.5 + degree(id,:);
        end
    end
end

% marker where the founder frequency comes from
% trace all closest genotyped descendants
source = zeros(nind, nind);
acc_founder = 0;
for i = 1:nind
    if( founder(i) )
        if( genotyped(i) )
            % if founder is genotyped, count its own frequency
            source(i,i) = 1;
        else
            geno_desc = 0;
            for j = 1:nind
                if( degree(j,i) ~= 0 && genotyped(j) == 1 )
                    geno_desc = geno_desc + degree(j,i);
                end
            end
            if( geno_desc <= 0 )
                continue;
            end
            for j = 1:nind
               if( degree(j,i) ~= 0 && genotyped(j) == 1 )
                   source(j,i) = degree(j,i)/geno_desc;
               end
            end
        end
        acc_founder = acc_founder + 1;
    end   
    % num go founders accessible to genotyped individuals
end

weight = zeros(nind, 1);
geno_count = 0;
if( acc_founder > 0 )
    for i = 1:nind
        weight(i) = sum(source(i,:));
        if( weight(i) > 0 )
            geno_count = geno_count + 1;
        end
    end
end



end











