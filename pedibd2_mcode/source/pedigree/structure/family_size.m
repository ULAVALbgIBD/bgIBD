

function [stat_family_size pedigree error] = family_size(pedigree)

% split families into connected components
% statistics about family size
% assign new family id

error = 0;


[rows, cols] = size(pedigree);

stat_family_size(1:rows) = 0;

if( rows == 0 || cols < 7 )
    error = 1;
    disp('error in pedigree structure');
    return;
end

if( cols >= 10 )
    % column 10 is reserved for relatedness
    if( nnz(pedigree(:,10)) > 0 )
        error = 1;
        disp('error in pedigree structure');
        return;
    end
end

if( cols < 12 )
    % column 12 is reserved for original family structure
    error = 1;
    disp('error in pedigree structure');
    return;
end

ncols = max(cols, 12);


related = 10;
new_f = 0;



pedigree(:,12) = pedigree(:,1);

list_f = unique(pedigree(:,12));

for f = list_f'
    
    f_index = find(pedigree(:,12) == f);
    if( isempty(f_index) )
        continue;
    end
    
    % find all disconnected components
    % in each original family
    family = pedigree( f_index, : );
    % initialize all to be unrelated
    family(:,related) = 1:length(family(:,1));
    for i = 1:length(family(:,1))
        id = family(i,2);
        if( id <= 0 )
            error = 1;
            disp('error in pedigre structures');
            return;
        end
        father = family(i,3);
        mother = family(i,4);
        if( father ~= 0 )
            rep0 = family(id,related);
            rep1 = family(father,related);
            if( rep0 < rep1 )
                temp = ( family(:,related) == rep1 );            
                family(temp,related) = rep0;
            end
            if( rep0 > rep1 )
                temp = ( family(:,related) == rep0 );            
                family(temp,related) = rep1;
            end            
        end
        if( mother ~= 0 )
            rep0 = family(id,related);
            rep2 = family(mother,related);
            % ordering is not necessary, just rank by founder order
            if( rep0 < rep2 )
                temp = ( family(:,related) == rep2 );            
                family(temp,related) = rep0;
            end
            if( rep0 > rep2 )
                temp = ( family(:,related) == rep0 );            
                family(temp,related) = rep2;
            end            
        end      
    end
    
    components = unique(family(:,related));
    if( isempty(components) )
        error = 1;
        disp('error in processing families');
        return;
    else
        nSUBfamily = length(unique(family(:,related)));
        if(  nSUBfamily > 1 )
            disp(['warning: family ', num2str(f), ': not all family members are related']);
            disp(['warning: family ', num2str(f), ' has ', num2str(nSUBfamily), ' disconnected portions']);                        
            disp(' ');
        end
    end
    
    % re-number family id in each original family
    pedigree( f_index, related ) = family(:, related);
    
    c1 = 0;
    tp = [];
    for j = components'
        t1 =  find(family(:,related) == j) ;
        if( isempty(t1) )
            error = 1;
            disp('error in processing families');
            return;
        end
        map(t1) = 1:length(t1);
        t2 = find( family(t1,7) == 1 );
        %partition families and reorder;
        if( ~isempty(t1) )
            new_f = new_f + 1;
            family( t1, 1 ) = new_f;
            tp(c1+1:c1+length(t1), 1:ncols) = family( t1, 1:ncols );
            
            % transform all family ids
            for k = 1:length(t1)
                for r = 2:4
                    if( tp(c1+k,r) ~= 0 )
                        tp(c1+k,r) = map(tp(c1+k,r));
                    end
                end
            end            
            c1 = c1+length(t1);          
        end
        
        if( ~isempty(t2) )
            % record number of typed individuals in a family
            stat_family_size(length(t2)) = stat_family_size(length(t2)) + 1;
        end
    end
    pedigree( f_index, 1:ncols ) = tp( 1:c1, 1:ncols );
end











