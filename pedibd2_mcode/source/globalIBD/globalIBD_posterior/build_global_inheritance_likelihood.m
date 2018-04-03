% build global relationship from pairwise relationship

function [ stat_assignment_all, stat_assignment_code error overflow ] = build_global_inheritance_likelihood( pairs, triples, map )

stat_assignment_all = [];
stat_assignment_code = [];
error = 0;
overflow = 0;
global debug_mode;

markers = pairs.intervals;

[nseg, cols] = size(markers);
if( nseg <= 0 || cols ~= 3 )
    error = 1;
    return;
end

ngeno = length(map.id_genotyped);
pair_pos = pairs.pair_pos;
pair_max = pairs.pair_max;
pair_posIBD1 = pairs.pair_posIBD1;
pair_vit = pairs.pair_vit;

if( ndims(pair_pos) ~= 4 || any(size(pair_pos) ~= [nseg, ngeno, ngeno, 3]) )
    error = 1;
    return;
end

if( ndims(pair_max) ~= 3 || any(size(pair_max) ~= [nseg, ngeno, ngeno]) )
    error = 1;
    return;
end

if( ndims(pair_vit) ~= 3 || any(size(pair_vit) ~= [nseg, ngeno, ngeno]) )
    error = 1;
    return;
end

if( ndims(pair_posIBD1) ~= 5 || any(size(pair_posIBD1) ~= [nseg, ngeno, ngeno, 2, 2]) )
    error = 1;
    return;
end

if( ngeno > 0 && any(size(triples.triple_viewpair) ~= [nseg,ngeno,ngeno,ngeno]) )
    error = 1;
    return;
end


% initialize the first assignment
[stat_assignment error] = assign_allele_random(map);
if( error ~= 0 || isempty(stat_assignment) )
    error = 1;
    disp('error in generating initial IBD assignment');
    return;
end

[nind, ~] = size(stat_assignment);

% state code
recombination = 1;
diff2ml = 2;
likelihood = 3;
complexity = 4;
solutions = 5;


stat_assignment_code = zeros(nseg, 5);
all_assignment = zeros(nseg, nind, 2);

seg_triples.num_triples = triples.num_triples;
seg_triples.ambig_triples = triples.ambig_triples;

% squeeze may accidentally transpose the matrix
for seg = 1:nseg


    stat_relation_prob = reshape(pair_pos(seg,1:ngeno,1:ngeno,1:3), [ngeno,ngeno,3]);
    stat_relation = reshape(pair_vit(seg,1:ngeno,1:ngeno), [ngeno,ngeno]);
    stat_relation = reshape(pair_max(seg,1:ngeno,1:ngeno), [ngeno,ngeno]);
    prob_IBD1 = reshape(pair_posIBD1(seg,1:ngeno,1:ngeno,1:2,1:2), [ngeno,ngeno,2,2]);

    seg_triples.triple_viewpair = reshape(triples.triple_viewpair(seg,1:ngeno,1:ngeno,1:ngeno), [ngeno,ngeno,ngeno]);
    
    % time likelihood together to optimize, branch and bound search
    % use former configuration to quicken search
    
    disp(['          processing segment ', num2str(seg), ' ...']);
    [stat_assignment, state_code, error, overflow] = assign_allele_nonrecur_likelihood(stat_relation, stat_relation_prob, prob_IBD1, seg_triples, map, stat_assignment);
    
    if( error ~= 0 )
        disp('error in generating global IBD');
        return;
    end
    
    if( overflow ~= 0 )
        return;
    end
    
    if( debug_mode == 1 )
        disp(['          complexity: ', num2str(state_code(complexity))]);
    end
    if( error ~= 0 )
        % complexity too high;
        if( seg == 1 )
            disp(['extremely high IBD complexity estimated for low marker density']);
        end
    end
    
    stat_assignment_code(seg,1:5) = state_code;
    all_assignment(seg, 1:nind, 1:2) = stat_assignment;
      
end

stat_assignment_all.all = all_assignment; 


end

