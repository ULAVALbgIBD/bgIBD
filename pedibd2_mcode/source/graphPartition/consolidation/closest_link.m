function [T error] = closest_link(T, group1, group2, affinity, fconf)

error = 0;

num = length(T);
if( num <= 0 )
    error = 1;
    disp('no chromosomal identity array');
    return;
end

if( ~isa(group1, 'logical') || ~isa(group2, 'logical') || ~isa(fconf, 'logical') )
    error = 1;
    disp('cluster type error');
    return;
end

if( length(group1) ~= num || length(group2) ~= num )
    error = 1;
    disp('cluster size not matched');
    return;
end

if( any(T > num | T < 0) )
    error = 1;
    disp('cluster id invalid');
    return;
end


[d1 d2] = size(affinity);
if( d1 <= 0 || d1 ~= num || d2 ~= d1 )
    error = 1;
    disp('chromosomal affinity matrix error');
    return;
end

[d1 d2] = size(fconf);
if( d1 <= 0 || d1 ~= num || d2 ~= d1 )
    error = 1;
    disp('chromosomal exclusion matrix error');
    return;
end

if( any(group1 & group2) )
    return;
end

% mark positions of chromosomal identities with occurences in either group
id1(1:num) = false;
id2(1:num) = false;
id1(T(group1)) = true;
id2(T(group2)) = true;

if( any(id1 & id2) )
    return;
end

% enumerating on chromosomal identity, representative chromosome
coherence = 0;
maxc = 0;
mi = 0;
mj = 0;
% mark the chromosomal identities with closest affinity
for i = 1:num
    if( ~id1(i) )
        continue;
    end
    % actual chromosome list;
    list1 = (T == i);
    for j = 1:num
        if( ~id2(j) )
            continue;
        end
        list2 = (T == j);
        if( any(any(fconf(list1,list2))) )
            continue;
        end
        coherence = mean(mean(affinity(list1,list2)));
        if( coherence > maxc )
            maxc = coherence;
            mi = i;
            mj = j;
        end
    end
end


% merge two identity groups
if( maxc > 0 && mi > 0 && mj > 0 )
    if( mi < mj )
        T( T == mj ) = mi;
    else
        T( T == mi ) = mj;
    end
end


end















