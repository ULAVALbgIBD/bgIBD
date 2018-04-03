

function lk_result = linkage_compute_allele(family, ind1, ind2)

% allele level kinship can be reduced to individual level kinship by
% considering their parents
% but need to treat ancestor-offspring relationship as special cases

max_path = 20;
f1 = family(ind1,3);
m1 = family(ind1,4);
f2 = family(ind2,3);
m2 = family(ind2,4);

parents = [f1,m1,f2,m2];


mapping(1,2) = 1;   
mapping(1,3) = 2;   % paternal sharing
mapping(1,4) = 3;
mapping(2,3) = 4;
mapping(2,4) = 5;   % maternal sharing
mapping(3,4) = 6;

mapping(2,1) = 1;
mapping(3,1) = 2;
mapping(4,1) = 3;
mapping(3,2) = 4;
mapping(4,2) = 5;
mapping(4,3) = 6;

lk_result(1:6,1:max_path) = 0;

% not correct for ancestor relationships

temp_family = family;
temp_family(ind1,3) = 0;
temp_family(ind1,4) = 0;
temp_family(ind2,3) = 0;
temp_family(ind2,4) = 0;

% make parents unconnected

for i = 1:2
    % common ancester            
    for j = 3:4
        p1 = parents(i);  
        p2 = parents(j);
        if( p1 ~= 0 && p2 ~= 0 )
            temp = linkage_compute(temp_family, p1, p2, max_path);
            lk_result(mapping(i,j), 1+2:max_path) = temp(1:max_path-2);
        end
    end
end
for i = 1:4
    % ancestor offspring
    p1 = parents(i);
    if( p1 ~= 0 )    
        if( i == 3 || i == 4 )
            temp = linkage_compute(temp_family, p1, ind1, max_path);
            lk_result(mapping(i,1),1+1:max_path) = lk_result(mapping(i,1),1+1:max_path) + 0.5 * temp(1:max_path-1);
            lk_result(mapping(i,2),1+1:max_path) = lk_result(mapping(i,2),1+1:max_path) + 0.5 * temp(1:max_path-1);
        end
        if( i == 1 || i == 2 )
            temp = linkage_compute(temp_family, p1, ind2, max_path);
            lk_result(mapping(i,3),1+1:max_path) = lk_result(mapping(i,3),1+1:max_path) + 0.5 * temp(1:max_path-1);
            lk_result(mapping(i,4),1+1:max_path) = lk_result(mapping(i,4),1+1:max_path) + 0.5 * temp(1:max_path-1);
        end
    end
end
for i = [1,3]
    % inbreeding
    % use original family structures
    % the connection between two parents may require the connection around
    % the other pair of parents to be present
    j = i + 1;
    p1 = parents(i);
    p2 = parents(j);
    if( p1 ~= 0 && p2 ~= 0 )
        temp = linkage_compute(family, p1, p2, max_path);
        lk_result(mapping(i,j), 1+2:max_path) = temp(1:max_path-2);
    end
end


lk_result(:,1:max_path-1) = lk_result(:,1+1:max_path);

%lk_result(:,i) number of path of length i


end



