
% this reference check requires all allele assignment 
% is indexed by actual individual ids

% this is efficient for dis-joint set management
% but not necessary where complexity is no concern

function [error] = check_1segment(alleles, family)

error = 0;
global debug_mode;

if( isempty(alleles) )
    error = 1;
    disp('error in global IBD');
    return;
end

if( isempty(family) )
    error = 1 ;
    disp('error in family structure');
    return;
end

[nind, ncols] = size(family);
if( nind <= 0 || ncols < 12 )
    error = 1;
    disp('error in family structure');
    return;
end

[r, c] = size( alleles );
if( r ~= nind || c ~= 2 )
    error = 1;
    disp('error in global IBD');
    return;
end


ids = unique(alleles);
if( any(abs(ids) > nind) )
    error = 1;
    disp('error in global IBD');
    return;
end

for i = 1:nind
    if( alleles(i,1) == i || alleles(i,2) == -i )
        error = 1;
        disp('error in global IBD: circular relationship');
        return;
    end
    % circular relationship
end

% relationship not updated to immediate representative
for i = 1:nind
    p = alleles(i,1);
    m = alleles(i,2);
    if( p ~= 0 )
        if( p ~= -i )
            if( p > 0 )
                if( alleles(p,2) ~= p )
                    error = 1;
                    disp('error in global IBD');
                    return;
                end
            end
            if( p < 0 )
                if( alleles(-p,1) ~= p )
                    error = 1;
                    disp('error in global IBD');
                    return;
                end
            end
        end
    end
    if( m ~= 0 )
        if( m ~= i )
            if( m > 0 )
                if( alleles(m,2) ~= m )
                    error = 1;
                    disp('error in global IBD');
                    return;
                end
            end
            if( m < 0 )
                if( alleles(-m,1) ~= m )
                    error = 1;
                    disp('error in global IBD');
                    return;
                end
            end
        end
    end
end



founder(1:nind) = 0;

for i = 1:nind
    if( family(i,3) == 0 && family(i,4) == 0 )
        founder(i) = 1;
    end
end

assigned(1:nind) = 0;
for i = 1:nind
    if( alleles(i,1) ~= 0 && alleles(i,2) ~= 0 )
        assigned(i) = 1;
    end
end



founder_alleles = unique(alleles(founder == 1 & assigned == 1));
if( length(founder_alleles) ~= nnz(founder == 1 & assigned == 1) )
    if( debug_mode == 1 )
        %error = 1;
        disp('error in global IBD');
        disp('founders related');
        %return;
    end
end




end
























