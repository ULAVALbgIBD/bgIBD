function [ output_assignment error ] ...
    = generate_consistent_assignment_2methods ...
    ( pairs, triples, range, ...
    kinship2, singleLINEAL, ...
    option )

error = 0;
output_assignment = [];


if( isempty(range) )
    error = 1;
    disp('error in family structures')
    return;
end

nIND = length(range.pedigree_range_full);
if( nIND < 1 )
    error = 1;
    disp('error in family structures');
    return;
end

if( isempty(pairs) )
    error = 1;
    disp('error in segmenting chromosomes');
    return;
end

nGENO = length(range.family_range);

% unGENOtyped individual is not assigned inheritance
% inheritance left as 0
% in enumeration method, inheritance will be randomly assigned in some
% is not constrained for unGENOtyped individuals
if( nGENO < 2 )
    if( ~isempty(kinship2) )
        disp('error in generating global IBD');
        error = 1;
        return;
    end
    intervals = pairs.intervals;
    [nseg, cols] = size( intervals );
    if( nseg ~= 1 || cols ~= 3 )
        error = 1;
        return;
    end
    config(1:nIND,1:2) = 0;
    assignment = [];    
    assignment.intervals(1,1:3) = intervals(1, 1:3);
    if( ~isempty(range.family_range) )
        config(range.family_range, 1) = (-1) * range.family_range;
        config(range.family_range, 2) = range.family_range;
    end
    assignment.alleles_all(1,1:nIND,1:2) = config;
    if( ~isempty(range.family_range) )
        assignment.alleles(1, 1:nGENO, 1:2) = config(range.family_range, 1:2);
    else
        assignment.alleles = [];
    end    
    output_assignment = assignment;
    return;
end


display(' ');
display(['family ',num2str(range.family_id), ': separating homologous chromosomes']);
time = cputime;


if( option(2) == 1 )

    if( nGENO <= 10 && nIND <= 18 )
        if( option(4) == 0 )
            disp('warning: low marker density: exact ibd mode enforced');
        end
        [ assignment error overflow ] ...
            = generate_consistent_assignment_likelihood( pairs, triples, range );
        if( error ~= 0 )
            disp('error in generating global IBD');
            return;
        end
    else
        overflow = 1;
    end

    if( overflow == 1 )
    % use graph partition for dense markers only
    % for bigger than 20 genotyped individuals, use graph partition
        [ assignment error ] ...
            = generate_consistent_assignment_partition ...
            ( pairs, triples, ...
            kinship2, ...
            singleLINEAL, ...
            [], 0, range );
        if( error ~= 0 )
            disp('error in generating global IBD');
            return;
        end
    end


end


% consider adding ibd re-mapping to reduce recombination changes

assignment.alleles = assignment.alleles_all(:,range.family_range,1:2);

difference = difference2viterbi(assignment.alleles, pairs.pair_vit);
assignment.diff = difference;

output_assignment = assignment;

display(['separating costs time: ', num2str(cputime - time), ' seconds']);

end






