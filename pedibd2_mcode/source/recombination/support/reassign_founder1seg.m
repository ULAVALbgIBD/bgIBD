
function [new_alleles error] = reassign_founder1seg(alleles, family)

error = 0;
new_alleles = [];

global debug_mode;

error = check_1segment(alleles, family);
if( error ~= 0 )
    disp('error in global IBD');
    return;
end

[nind, fields] = size(family);
if( nind <= 0 || fields < 12 )
    error = 1;
    return;
end

[nrows, ncols] = size(alleles);
if( nrows ~= nind || ncols ~= 2 )
    error = 1;
    return;
end

founder(1:nind) = 0;
for i = 1:nind
    if( family(i,3) == 0 && family(i,4) == 0 )
        founder(i) = 1;
    end
end

new_alleles = alleles;
for i = 1:nind
    if( founder(i) == 1 )
        if( new_alleles(i,1) ~= 0 && new_alleles(i,1) ~= -i )
            p = new_alleles(i,1);
            temp_alleles = new_alleles;
            temp_alleles(new_alleles(:,1) == p, 1) = -i;
            temp_alleles(new_alleles(:,2) == p, 2) = -i;
            new_alleles = temp_alleles;
        end
        if( new_alleles(i,2) ~= 0 && new_alleles(i,2) ~= i )
            m = new_alleles(i,2);
            temp_alleles = new_alleles;
            temp_alleles(new_alleles(:,1) == m, 1) = i;
            temp_alleles(new_alleles(:,2) == m, 2) = i;
            new_alleles = temp_alleles;
        end
    end
end


for i = 1:nind    
    if( founder(i) == 1 )
        p = new_alleles(i,1);
        m = new_alleles(i,2);        
        if( p ~= 0 && m ~= 0 )
            if( p ~= -i || m ~= i )
                if( debug_mode == 1 )
                    error = 1;
                    disp('conflict in assigning founder alleles');
                    disp('founders may be related');
                    return;
                end
            end
        end
    end
end

end
