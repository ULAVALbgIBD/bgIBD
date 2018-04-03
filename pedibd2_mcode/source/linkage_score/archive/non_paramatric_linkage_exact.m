

function stat_score = non_paramatric_linkage_exact(stat_interval, stat_assignment_all, input_family_range, pedigree_all_missing)


stat_score(1:length(stat_interval(:,1))) = 0;


affected = [];
count = 0;
for i = 1:length(input_family_range)
    if( pedigree_all_missing(input_family_range(i),6) == 2 )
        count = count + 1;
        affected(count) = i;     
    end
end 

for i = 1:count
    fac(i) = factorial(i);
end

allele_configure = [];
configure_count = [];

time = cputime;

for l = 1:length(stat_interval(:,1))
    
    stat_assignment = stat_assignment_all{l}(affected,1:2);
    t1 = min(min(stat_assignment));
    t2 = max(max(stat_assignment));
    range = t2-t1+1;
    offset = t1;
    stat_assignment = stat_assignment - offset + 1;
    

    allele_configure(1,1:range) = 0;
    configure_count(1) = 1;
    index = 1;
    
    for i = 1:length(affected)
        p = stat_assignment(i,1);
        m = stat_assignment(i,2);
        if( p ~= m )
            allele_configure(index+1:2*index,:) = allele_configure(1:index,:);
            configure_count(index+1:2*index) = configure_count(1:index);
            allele_configure(1:index,p) = allele_configure(1:index,p) + 1;
            allele_configure(index+1:2*index,m) = allele_configure(index+1:2*index,m) + 1;
            index = 2*index;
        else
            allele_configure(1:index,p) = allele_configure(1:index,p) + 1;
            configure_count(1:index) = 2*configure_count(1:index);
        end
    end
    
    total = 0;
    for i = 1:index
        temp = 1;
        for k = 1:range
            if( allele_configure(i,k) ~= 0 )
                temp = temp * fac(allele_configure(i,k));
            end
        end
        temp = temp * configure_count(i);
        total = total + temp;
    end
    
    l;
    stat_score(l) = total/(2^length(affected));
    disp = stat_score(l);
    
end

display(['nonparametric linkage computing time: ', num2str(cputime - time), ' seconds']);


end




