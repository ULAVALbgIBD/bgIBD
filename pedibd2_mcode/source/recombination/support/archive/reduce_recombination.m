
function [new_alleles, error] = reduce_recombination(family, valid, alleles_all)

error = 0;
new_alleles = [];


[nind, fields] = size(family);
if( nind <= 0 || fields < 12 )
    error = 1;
    return;
end


if( ndims(alleles_all) ~= 3 )
    error = 1;
    disp('error in global IBD');
    return;
else
    [nseg, d2, d3] = size(alleles_all);
    if( nseg <= 0 || d2 ~= nind || d3 ~= 2 )
        error = 1;
        disp('error in global IBD');
        return;
    end
end


len = length(valid);
if( len <= 0 || len ~= nseg )
    error = 1;
    return;
end

new_alleles = alleles_all;

if( nseg < 2 )
    return;
end


% let founder alleles be representatives
for i = 1:nseg
    config = reshape(alleles_all(i, 1:nind, 1:2), [nind,2]);
    [config, error] = reassign_founder1seg(config, family);
    if( error ~= 0 )
        disp('error in global IBD');
        return;
    end
    new_alleles(i, 1:nind, 1:2) = config;    
end

if( ~any(valid) )
    return;
end


% reduce recombination between nearby segments
for i = 1:nseg
    if( valid(i) == 1 )
        first_seg = i;
        break;
    end
end
pre_seg = first_seg;
pre_config = reshape(new_alleles(pre_seg, 1:nind, 1:2), [nind,2]);
for i = first_seg+1:nseg
    if( valid(i) ~= 1 )
        %continue;
    end
    config = reshape(new_alleles(i, 1:nind, 1:2), [nind,2]);
    [r1, c1] = size(pre_config);
    [r2, c2] = size(config);
    if( r1 <= 0 || r1 ~= r2 || r1 ~= nind || c1 ~= 2 || c2 ~= 2 )
        error = 1;
        disp('error in global IBD');
        return;
    end
    [config, nrec, error] = reduce_recombination1seg(pre_config, config, family);
%     disp([num2str(i), ' ', num2str(tdiff), ' ', num2str(nrec)]);
    if( error ~= 0 )
        disp('error in global IBD');
        return;
    end
    new_alleles(i, 1:nind, 1:2) = config;
    pre_config = config;
end


end










