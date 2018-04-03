function [p_clique, m_clique, error] = kinship_partition(index, neighbor1, kinship2, posterior)

% partition neighbors of index
% 1 according to predetermined parental origin, all non-inbreeding
% ancestors going through determined parent
% 2 according to very confident siblings with typed parents
% 3 some siblings and all children will remain undetermined


error = 0;
p_clique = [];
m_clique = [];

[d1, d2, d3, d4] = size(kinship2);
if( d1 <= 0 || d1 ~= d2 || d3 ~= 2 || d4 ~= 2 )
    error = 1;
    disp('error in kinship');
    return;
end
num = d1;

[d1, d2, d3, d4] = size(posterior);
if( d1 ~= num || d2 ~= d1 || d3 ~= 2 || d4 ~= 2 )
    error = 1;
    disp('error in posterior sharing status');
    return;
end

if( index > num || index <= 0 )
    error = 1;
    disp('wrong individual id');
    return;
end

p(1:num,1:2) = kinship2(index,1:num,1,1:2);
m(1:num,1:2) = kinship2(index,1:num,2,1:2);

temp1 = sum( p(1:num,1:2) , 2);
temp2 = sum( m(1:num,1:2) , 2);
temp3 = temp1+temp2;

p_only = find((temp1~=0) & (temp2==0));
m_only = find((temp1==0) & (temp2~=0));
pm_both = find((temp1~=0) & (temp2~=0));

if( ~isempty(p_only) )
    [X,I] = sort(temp1(p_only), 'descend');
    p_only = p_only(I);
end
if( ~isempty(m_only) )
    [X,I] = sort(temp2(m_only), 'descend');
    m_only = m_only(I);
end
if( ~isempty(pm_both) )
    [X,I] = sort(temp3(pm_both), 'descend');
    pm_both = pm_both(I);
end

% two cliques are generated out of neighbor1

for i = 1:length(p_only)
    if( ismember(p_only(i), neighbor1) )
        p_clique = union(p_clique, p_only(i));
    end
end


for i = 1:length(m_only)
    if( ismember(m_only(i), neighbor1) )
        m_clique = union(m_clique, m_only(i));
    end
end


% process predisposed siblings, because of both typed parents

p_post(1:num,1:2) = posterior(index,1:num,1,1:2);
m_post(1:num,1:2) = posterior(index,1:num,2,1:2);


% single lineal relationships are already resolved
threshold = 10;
for i = 1:length(pm_both)
    if( ismember(pm_both(i), neighbor1) )
        p_kinship_temp = sum(p(pm_both(i),1:2));
        m_kinship_temp = sum(m(pm_both(i),1:2));
        if( p_kinship_temp <= 0 || m_kinship_temp <= 0 )
            % pm_both must have both prios non-zeros
            disp('error in sharing status');
        end
        p_post_temp = sum(p_post(pm_both(i),1:2));
        m_post_temp = sum(m_post(pm_both(i),1:2));        
        ptemp = p_post_temp;
        mtemp = m_post_temp;
        if( ptemp / mtemp > threshold )
            p_clique = union(p_clique, pm_both(i));
        end
        if( mtemp / ptemp > threshold )
            m_clique = union(m_clique, pm_both(i));
        end
    end
end


end

