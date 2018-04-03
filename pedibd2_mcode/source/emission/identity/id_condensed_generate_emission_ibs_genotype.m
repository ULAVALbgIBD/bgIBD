function ep_merge_ibs = id_condensed_generate_emission_ibs_genotype(ep_merge_genotype)

% designed only for 4 alleles

num_is = length(ep_merge_genotype(:,1));

ep_merge_ibs(1:num_is, 1:4) = 0;

ep_merge_ibs(1:num_is,4) = ep_merge_genotype(:,17);

for i = 1:num_is
    for j = 0:15
        temp1 = 0;
        temp2 = 0;
        if( bitget(j,1) == bitget(j,3) )
            temp1 = temp1 + 1;
        end
        if( bitget(j,1) == bitget(j,4) )
            temp2 = temp2 + 1;
        end
        if( bitget(j,2) == bitget(j,4) )
            temp1 = temp1 + 1;
        end
        if( bitget(j,2) == bitget(j,3) )
            temp2 = temp2 + 1;
        end
        temp = max(temp1,temp2) + 1;
        ep_merge_ibs(i,temp) = ep_merge_ibs(i,temp) + ep_merge_genotype(i,j+1);
    end
end

end


