function plot_genotype_wrap_multi(input, output, parameters, f_id, list1, limit)


nfam = length(input.family_range);
if( nfam <= 0 )
    return;
end

f = 0;
for i = 1:nfam
    if( input.family_range{i}.family_id == f_id )
        f = i;
        break;
    end
end

if( f == 0 )
    return;
end

family = input.family_range{f}.structure;
genotyped = family(:,7) == 1;

list = [];
count = 0;
for i = 1:length(list1)
    id1 = list1(i);
    temp = family(:,8);
    in1 = find(temp == id1);
    if( length(in1) ~= 1 ) 
        continue;
    end
    count = count + 1;        
    list(count,1) = in1;
end



if( isempty(list1) )
    list = family(genotyped, 2);
end

if( isempty(list) )
    return;
end

scrsz = get(0,'ScreenSize');
figure('OuterPosition',[1 scrsz(4)/2 scrsz(3) scrsz(4)/2]);

global disp_mode;
if( disp_mode == 0 ) 
    coor = parameters.sampled_markerlist(:,1);
else
    coor = parameters.sampled_markerlist(:,2);
end

if( length(limit) ~= 2 )
    limit = [1,length(parameters.sampled_markerlist(:,1))];
else
    if( limit(1) > limit(2) || limit(1) < 1 || limit(2) > length(parameters.sampled_markerlist(:,1)) )
        return;
    end
end

error = plot_genotype(input.family_genotype{f}, family, list, coor, limit);

end













