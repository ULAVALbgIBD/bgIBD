

function [output error] = sibship_isomorphism(assignment, range)

error = 0;
sibship = [];
new_alleles = [];
global debug_mode;

[alleles_all, intervals, error] = check_allsegments(assignment, range);
if( error ~= 0 )
    disp('error in extracting recombination positions');
    return;
end

family = range.structure;

[nIND, fields] = size(family);
if( nIND <= 0 || fields < 12 )
    error = 1;
    disp('error in family structures');
    return;
end


if( ndims(alleles_all) ~= 3 )
    error = 1;
    disp('error in global IBD');
    return;
else
    [nSEG, d2, d3] = size(alleles_all);
    if( nSEG <= 0 || d2 ~= nIND || d3 ~= 2 )
        error = 1;
        disp('error in global IBD');
        return;
    end
end

[r, c] = size(intervals);
if( r <= 0 || r ~= nSEG || c < 2 )
    error = 1;
    disp('error in global IBD');
    return;
end




genotyped(1:nIND) = 0;
for i = 1:nIND
    if( family(i,7) == 1 )
        genotyped(i) = 1;
    end
end

ismother(1:nIND,1:nIND) = 0;
isfather(1:nIND,1:nIND) = 0;
isparent(1:nIND,1:nIND) = 0;
for i = 1:nIND
    id = i;
    father = family(id,3);
    mother = family(id,4);
    if( father ~= 0 && genotyped(id) == 1 ) 
        isfather(father,id) = 1;
    end
    if( mother ~= 0 && genotyped(id) == 1 )
        ismother(mother,id) = 1;
    end
end
for i = 1:nIND
    if( any(isfather(i,:) & ismother(i,:)) )
        error = 1;
        disp(['error in family structures: ', num2str(family(i,8)), 'being both father and mother']);
        return;
    end
    isparent(i,:) = isfather(i,:) | ismother(i,:);
end


sibship = zeros(nIND, nIND, nIND);
for i = 1:nIND
    for j = 1:nIND
        if( any( isfather(i,:) & ismother(j,:) ) )
            sibship(i,j, 1:nIND) = - (isfather(i,:) & ismother(j,:));
        end
        if( any( ismother(i,:) & isfather(j,:) ) )
            sibship(i,j, 1:nIND) = ismother(i,:) & isfather(j,:);
        end
    end
end


% if one of the parents is genotyped
% the phase should be determined, check for exact match

% check whether children's alleles appear in parents

all_config(1:nSEG, 1:nIND, 1:2) = 0;
sibship_valid(1:nIND,1:nIND,1:nSEG) = true;
for k = 1:nSEG
    config = reshape(alleles_all(k, 1:nIND, 1:2), [nIND,2]);
    for i = 1:nIND
        for j = 1:nIND
            children = [];
            children(1:nIND) = sibship(i,j,1:nIND);
            if( any( children ~= 0 ) )
                paternal_alleles = config( sibship(i,j,1:nIND) ~= 0, 1 );
                maternal_alleles = config( sibship(i,j,1:nIND) ~= 0, 2 );
                alleles = unique(union(paternal_alleles,maternal_alleles));
                num_alleles = length(alleles);
                if( num_alleles > 4 )
                    sibship_valid(i,j,k) = false;
                end
                all_config(k, children~=0, 1:2) = isomorphism_replace(config(children~=0, 1:2));
            end
            if( any( children < 0 ) )
                if( genotyped(i) == 1 )
                    if( ~all( paternal_alleles == config(i,1) | paternal_alleles == config(i,2) ) )
                        sibship_valid(i,j,k) = false;
                    end
                    if( length(unique(maternal_alleles)) > 2 )
                        sibship_valid(i,j,k) = false;
                    end
                end        
            end
            if( any( children > 0 ) )
                if( genotyped(i) == 1 )
                    if( ~all( maternal_alleles == config(i,1) | maternal_alleles == config(i,2) ) )
                        sibship_valid(i,j,k) = false;
                    end
                    if( length(unique(paternal_alleles)) > 2 )
                        sibship_valid(i,j,k) = false;
                    end
                end             
            end
        end
    end
end

valid(1:nSEG) = true;
for i = 1:nSEG
    valid(i) = all(all(sibship_valid(1:nIND,1:nIND,i)));
end
if( debug_mode == 1 )
    disp(['family ', num2str(range.family_id), ': total consistent loci: ', num2str(sum(intervals(valid>0,2)-intervals(valid>0,1)+1))]);
end


% gamete genitor

sibship_status(1:nIND,1:nIND,1:nSEG) = -1;
for i = 1:nIND
    for j = 1:nIND
        children = [];
        children(1:nIND) = sibship(i,j,1:nIND);
        nchi = nnz(children);
        if( nchi <= 0 )
            continue;
        end
        for k = 1:nSEG
            if( sibship_valid(i,j,k) )
                first_seg = k;
                break;
            end
        end
        if( k == nSEG )
            break;
        end
        pre_seg = first_seg;
        sibship_status(i,j,pre_seg) = 0;
        for k = first_seg + 1 : nSEG
            if( ~sibship_valid(i,j,k) )
                continue;
            end
            cur_seg = k;
            pre_config(1:nchi,1:2) = all_config(pre_seg,children~=0,1:2);
            cur_config(1:nchi,1:2) = all_config(cur_seg,children~=0,1:2); 
            if( genotyped(i) == 1 || genotyped(j) == 1 )
                sex = family(i,5);
                if( any(any(cur_config(1:nchi,sex)-pre_config(1:nchi,sex))) )
                    sibship_status(i,j,cur_seg) = mod(sibship_status(i,j,pre_seg)+1,2);
                else
                    sibship_status(i,j,cur_seg) = sibship_status(i,j,pre_seg);
                end
            end
            if( genotyped(i) == 0 && genotyped(j) == 0 )
                [swap_config, diff] = realign(pre_config(1:nchi,1:2), cur_config(1:nchi,1:2));
                if( diff > 0 )
                    sibship_status(i,j,cur_seg) = mod(sibship_status(i,j,pre_seg)+1,2);
                else
                    sibship_status(i,j,cur_seg) = sibship_status(i,j,pre_seg);
                end
            end
            pre_seg = cur_seg;
        end
    end
end


output.sibship_status = sibship_status;


if( debug_mode == 1 )
    figure;
    hold on;
    for i = 1:length(valid)
        line(intervals(i,1:2), [valid(i), valid(i)], 'linewidth', 2);
        text(mean(intervals(i,1:2)), 0.5, num2str(i));
    end
    ylim([-0.1,1.1]);
    hold off;
end




end




































