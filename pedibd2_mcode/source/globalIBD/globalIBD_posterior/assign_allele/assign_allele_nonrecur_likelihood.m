function [ max_assignment, state_code, error, overflow ] = assign_allele_nonrecur_likelihood( match, match_prob, prob_IBD1, triples, map, pre_all )

error = 0;
overflow = 0;
equal_assignment = [];
max_assignment = [];
state_code = [];

output_complexity = 0;
num_solutions = 0;

num_geno = length(map.id_genotyped);
num_all = length(map.id);
if( num_all <= 0 || num_geno > num_all || any(map.id_genotyped > num_all) )
    error = 1;
    disp('error in family structures');
    return;
end
[r, c] = size(match);
if( r ~= num_geno || c ~= num_geno )
    error = 1;
    disp('error in viterbi decoding');
    return;
end
[d1, d2, d3] = size( match_prob );
if( d1 ~= num_geno || d2 ~= num_geno || d3 ~= 3 )
    error = 1;
    disp('error in posterior probability');
    return;
end
[d1, d2, d3, d4] = size( prob_IBD1 );
if( d1 ~= num_geno || d2 ~= num_geno || d3 ~= 2 || d4 ~= 2 )
    error = 1;
    disp('error in posterior probability');
    return;
end

if( isempty(triples) )
    error = 1;
    disp('error in triplets');
    return;
end
num_triples = triples.num_triples;
triple_viewpair = triples.triple_viewpair;
valid_triples = triples.ambig_triples;
[d1, d2, d3] = size(valid_triples);
if( d1 <= 0 || d1 ~= d2 || d2 ~= d3 || d1 ~= num_geno )
    error = 1;
    disp('error in triplets');
    return;
end
if( num_triples ~= max(max(max(valid_triples))) )
    error = 1;
    disp('error in triplets');
    return;
end
[d1, d2, d3] = size(triple_viewpair);
if( d1 ~= d2 || d2 ~= d3 || d1 ~= num_geno )
    error = 1;
    disp('error in triplets');
    return;
end

[r, c] = size(pre_all);
if( r ~= num_all || c ~= 2 )
    error = 1;
    disp('error in former IBD');
    return;
end

output_diff = 0;

pre = pre_all(map.id_genotyped, 1:2);
max_likelihood = 1;

for i = 1:num_geno
    for j = 1:i-1
        ibd_num = count_ibd(pre(i,1:2), pre(j,1:2));
        if( ibd_num == 1 )
            % if there is directional information
            [a,b] = count_ibd1(pre(i,1:2), pre(j,1:2));
            max_likelihood = max_likelihood * prob_IBD1(i,j,a,b);
        else
            max_likelihood = max_likelihood * match_prob(i,j,ibd_num+1);
        end
        if( match(i,j) ~= ibd_num + 1 )
            output_diff = output_diff + 1;
        end
    end
end
output_likelihood = max_likelihood;

global assignment;
global all_assignment;
global mismatch;
global likelihood;
global paternal;
global maternal;
global paternal_num;
global maternal_num;

assignment(1:num_geno,1:2) = 0;
mismatch(1:num_geno) = 0;
likelihood(1:num_geno) = 1;
change = 2*num_geno;

current = 1;
all_assignment(1:num_all,1:2) = 0;
max_assignment = pre_all;
p(1:num_all) = 0;
m(1:num_all) = 0;
reset(1:num_all) = 1;

paternal(1:num_all,1:2) = 0;
maternal(1:num_all,1:2) = 0;
paternal_num(1:num_all) = 0;
maternal_num(1:num_all) = 0;

paternal_alleles = map.paternal_alleles;
maternal_alleles = map.maternal_alleles;
[r, c] = size(paternal_alleles);
if( r ~= num_all || c ~= 2 )
    error = 1;
    disp('error in family structures');
    return;
end
[r, c] = size(maternal_alleles);
if( r ~= num_all || c ~= 2 )
    error = 1;
    disp('error in family structures');
    return;
end
for i = 1:num_all
    if( any(paternal_alleles(i,1:2) == 0 | paternal_alleles(i,1:2) > i | paternal_alleles(i,1:2) < -i ) )
        error = 1;
        disp('error in family structures');
        return;
    end
    if( any(maternal_alleles(i,1:2) == 0 | maternal_alleles(i,1:2) > i | maternal_alleles(i,1:2) < -i ) )
        error = 1;
        disp('error in family structures');
        return;
    end
end

while(1)

    if( current < 1 )
        break;
    end
    
    if( map.relevance(current) ~= 1 )
        if( reset(current) == 1 )
            reset(current) = 0;
            current = current + 1;            
        else
            reset(current) = 1;
            current = current - 1;            
        end
        continue;
    end
    
    if( reset(current) == 1 )
        new_assignment(current, paternal_alleles, maternal_alleles);
        p(current) = 1;
        m(current) = 0;
        reset(current) = 0;
    end

    i = p(current);
    j = m(current);

    backward = false;
    
    while(1)

        correct = 1;
        if( j < maternal_num(current) )
            j = j + 1;
        else
            if( i < paternal_num(current) )
                i = i + 1;
                j = 1;
            else
                backward = true;
                break;
            end
        end
        
        output_complexity = output_complexity + 1;

        if( output_complexity > 10^7 )
            overflow = 1;
            return;
        end
        
        correct = compute_likelihood(current, i, j, match, match_prob, prob_IBD1, valid_triples, triple_viewpair, max_likelihood, map);
        
        order = map.reverse_list(current);
        if( correct == 1 && order == num_geno )
            
            t_change = nnz(pre_all(1:num_all,1:2) - all_assignment(1:num_all,1:2));          
            if( likelihood(num_geno) == max_likelihood )
                if( t_change >= change )
                    correct = 0;
                end
            end
            if( likelihood(num_geno) == max_likelihood )
                num_solutions = num_solutions + 1;                
            end
            if( likelihood(num_geno) > max_likelihood )
                num_solutions = 1;
            end
            if( correct == 1 )
                max_likelihood = likelihood(num_geno);
                max_assignment = all_assignment;
                output_likelihood = max_likelihood;
                output_diff = mismatch(num_geno);
                change = t_change;
            end
            
        end
        
        if( correct == 1 && order < num_geno )
            backward = false;
            break;
        end
        
    end
    
    p(current) = i;
    m(current) = j;
    
    if( backward )
        reset(current) = 1;
        current = current - 1;
    else
        current = current + 1;
    end
    
end

state_code = [change, output_diff, output_likelihood, output_complexity, num_solutions];

clear global assignment;
clear global all_assignment;
clear global mismatch;
clear global likelihood;
clear global paternal;
clear global maternal;
clear global paternal_num;
clear global maternal_num;

end






