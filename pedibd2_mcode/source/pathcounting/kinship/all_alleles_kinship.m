function [kinship, kinship2, kinship2ex, error] = all_alleles_kinship(all_inheritance, range, full_range)
    
    % specially designed for two alleles only
    error = 0;
    kinship = [];
    kinship2 = [];
    list = all_inheritance.list;
    if( isempty(list) )
        error = 1;
        disp('error in calculating kinship');
        return;
    end
    slots = size(list, 2);
    if( slots ~= 2 )
        error = 1;
        disp('error in calculating kinship');
        return;
    end
    list_prob = zeros(size(list,1), 1);
    
    nIND = length(full_range);
    nGENO = length(range);
    
    content = all_inheritance.content;
    fast_index = all_inheritance.fast_index;
    
    % check correcness of the data
    for ii = 1:nGENO
        id1 = range(ii);
        for jj = 1:nGENO
            id2 = range(jj);
            for i = [2*id1-1, 2*id1]
                for j = [2*id2-1, 2*id2]
                    if( fast_index(i,j) == 0 )
                        error = 1;
                        disp('error in calculating kinship');
                        return;
                    else
                        if( fast_index(i,j) > size(list,1) )
                            error = 1;
                            disp('error in calculating kinship');
                            return;
                        else
                            if( ~isempty(setxor([i,j], list(fast_index(i,j),:))) )
                                error = 1;
                                disp('error in calculating kinship');
                                return;
                            end
                        end
                    end
                end
            end
        end
    end
    
    % content format
    % first column is coefficient
    % 2~slots+1, is the allele involved
    % 2+slots~, is the power of each term
    
    k = zeros(nIND * 2, nIND * 2);
    for i = 1:size(list,1)
        temp = content{i};
        for j = 1:size(temp,1)
            if( temp(j,2) == temp(j,3) )
                % ibd
                list_prob(i) = list_prob(i) + temp(j,1) * ( 0.5 ^ (sum(temp(j,slots+2:end))) );
            end
        end
        k(list(i,1),list(i,2)) = list_prob(i);
        k(list(i,2),list(i,1)) = list_prob(i);
        % kinship between 2 alleles
    end
    
    kinship2ex = zeros(nIND, nIND, 2, 2);
    for i = 1:nIND
        for j = 1:nIND
            kinship2ex(i,j,1:2,1:2) = k(i*2-1 : i*2, j*2-1 : j*2);
        end
    end
    
    kinship2 = kinship2ex(range, range, 1:2, 1:2);
    kinship = 0.25 * sum(sum(kinship2(1:nGENO,1:nGENO, 1:2, 1:2), 4),3);
    
    
end



