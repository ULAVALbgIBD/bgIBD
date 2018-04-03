function [ assignment, consistency, typing_error, error ] = assign_locus( constraint, genotype )

error = 0;
assignment = [];
consistency = [];
typing_error = 0;

global debug_mode;

[row1, col1] = size(constraint);
[row2, col2] = size(genotype);
if( col1 ~= 2 || row1 <= 0 || row1 ~= row2 )
    error = 1;
    disp('genotype and IBD do not match');
    return;
end


num = row1; % number of all individuals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global ancestral;
% paternal side on the left
% maternal side on the right
% founder f := [-f,f]

% ancestral covers the whole family, indexed by local id
% constraint may only cover genotyped range

ancestral(1:num,1) = -1:-1:-num;
ancestral(1:num,2) = 1:num;

ancestral(1:num,3:4) = 0;
ancestral(1:num,5:6) = 0;

depth = 0;
consistency = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:num
    % not missing, homozygous    
    if( genotype(i,1) ~= 0 && genotype(i,2) ~= 0 )
        if( genotype(i,1) == genotype(i,2) )
            ancestral(i, 5:6) = genotype(i,1:2);
        end
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%update IBD constraints to ancestral
%constraint must be labeled by ancestors

for i = 1:num
    if( constraint(i,1) ~= 0 )
        if( constraint(i,1) ~= -i )
            consistency = consistency * join_set(-i, constraint(i,1), 0, true);
        end
    end
    if( constraint(i,2) ~= 0 )
        if( constraint(i,2) ~= i )
            consistency = consistency * join_set(i, constraint(i,2), 0, true);
        end
    end
    if( consistency == 0 )
        depth = i - 1;
        break;
    end
    depth = i;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculating vote of each allele in each representative
for i = 1:num
    update_ancestor(-i);
    update_ancestor(i);
end

[vote_assignment, typing_error, error] = assign_locus_vote(ancestral(1:num,1:2), genotype);
if( error ~= 0 )
    disp('error in global IBD');
    return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%update IBS constraints to ancestral

if( consistency == 1 )
    for i = 1:num
        if( genotype(i,1) ~= 0 && genotype(i,2) ~= 0 )
            if( genotype(i,1) == genotype(i,2) )
                [consistency] = join_set(-i, i, 0, false);
                if( consistency == 1 )
                    [consistency] = update_ancestor(i);
                    [anc, offset, assign] = find_immediate(i);
                    [consistency] = set_assign(anc, genotype(i,1), offset);
                end
            else
                [consistency] = join_set(-i, i, 1, false);
            end
        end
        if( consistency == 0 )
            depth = i - 1;
            break;
        end
        depth = i;
    end
end

if( depth == 0 )
    depth;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%generate final haplotype solution
assignment(1:num,1:2) = 0;

if( consistency == 1 )
    for i = 1:num
        update_ancestor(-i);
        update_ancestor(i);
        assignment(i, 1:2) = ancestral(i, 5:6);
    end
    if( debug_mode == 1 )
        if( typing_error ~= 0 )
            disp(['          consistent assignment against majority vote: ', num2str(typing_error)]);
        end
    end
else
    assignment = vote_assignment;   
end


consistency = depth;

clear global ancestral;

end

