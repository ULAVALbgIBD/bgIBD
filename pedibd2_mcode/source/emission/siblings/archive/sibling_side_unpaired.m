function e = sibling_side_unpaired(f_geno, m_geno, freq)

% ibd not because of sharing parental alleles
% farther sharing

% not yet considering sharing between one individual

% a1:a2, a2:b2
% a1:a2, b1:a1
ibd1 = 2;
% a1:a2, a2:a1
ibd2 = 3;


ibs0 = 1;
ibs1 = 2;
ibs2 = 3;
ibs3 = 4; %missing

% change this, 1/2 does not reflect actual alleles, 1 is major allele

p1 = freq(1);
p2 = freq(2);

e(1:4,1:3) = 0;

if( f_geno(1) ~= 0 && f_geno(2) ~= 0 && m_geno(1) ~= 0 && m_geno(2) ~= 0 )
    if( f_geno(1) == f_geno(2) )
        if( m_geno(1) == m_geno(2) )
            % aa, cc
            if( g_geno(1) == m_geno(1) )
                e(ibd1, ibs0) = 0;
                e(ibd1, ibs1) = 0;
                e(ibd1, ibs2) = 1;
                e(ibd2, ibs0) = 0;
                e(ibd2, ibs1) = 0;
                e(ibd2, ibs2) = 1;
            else
                e(ibd1, ibs0) = 0;
                e(ibd1, ibs1) = 0;
                e(ibd1, ibs2) = 0;
                e(ibd2, ibs0) = 0;
                e(ibd2, ibs1) = 0;
                e(ibd2, ibs2) = 0;
            end            
        else
            % aa, cd
            if( g_geno(1) == m_geno(1) )
                e(ibd1, ibs0) = 0;
                e(ibd1, ibs1) = 0;
                e(ibd1, ibs2) = 1;
            else
                e(ibd1, ibs0) = 0;
                e(ibd1, ibs1) = 0;
                e(ibd1, ibs2) = 0;
            end             
            e(ibd2, ibs0) = 0;
            e(ibd2, ibs1) = 0;
            e(ibd2, ibs2) = 0; 
        end
    else
        if( m_geno(1) == m_geno(2) )
            % ab, cc
            if( g_geno(1) == m_geno(1) )
                e(ibd1, ibs0) = 0;
                e(ibd1, ibs1) = 1;
                e(ibd1, ibs2) = 0;                 
            else
                e(ibd2, ibs0) = 0;
                e(ibd2, ibs1) = 0;
                e(ibd2, ibs2) = 0;                 
            end
            e(ibd2, ibs0) = 0;
            e(ibd2, ibs1) = 0;
            e(ibd2, ibs2) = 0;            
        else
            % ab, cd
            % besides this situation, all situations are absorbed by ibd0
            e(ibd1, ibs0) = 0;
            e(ibd1, ibs1) = 0;
            e(ibd1, ibs2) = 1;
            e(ibd2, ibs0) = 0;
            e(ibd2, ibs1) = 0;
            e(ibd2, ibs2) = 1;
        end    
    end
end

% father missing
if( (f_geno(1) == 0 || f_geno(2) == 0) && m_geno(1) ~= 0 && m_geno(2) ~= 0 )
    if( m_geno(1) == m_geno(2) )
        % ?? aa

        if( m_geno(1) == 1 )
            e(ibd1, ibs0) = 0;
            e(ibd1, ibs1) = p2;
            e(ibd1, ibs2) = p1;
        else
            e(ibd1, ibs0) = 0;
            e(ibd1, ibs1) = p1;
            e(ibd1, ibs2) = p2;            
        end
        
        e(ibd2, ibs0) = 0;
        e(ibd2, ibs1) = 0;
        e(ibd2, ibs2) = 1;
        
    else
        % ?? cd
        e(ibd1, ibs0) = 0;
        e(ibd1, ibs1) = 0.5;
        e(ibd1, ibs2) = 0.5;
        
        e(ibd2, ibs0) = 0;
        e(ibd2, ibs1) = 0;
        e(ibd2, ibs2) = 1;        
      
    end
end

% mother missing
if( f_geno(1) ~= 0 && f_geno(2) ~= 0 && ( m_geno(1) == 0 || m_geno(2) == 0 ) )
    if( f_geno(1) == f_geno(2) )
        %aa ??
        if( f_geno(1) == 1 )
            e(ibd1, ibs0) = 0;
            e(ibd1, ibs1) = p2;
            e(ibd1, ibs2) = p1;
        else
            e(ibd1, ibs0) = 0;
            e(ibd1, ibs1) = p2;
            e(ibd1, ibs2) = p1;            
        end
        e(ibd2, ibs0) = 0;
        e(ibd2, ibs1) = 0;
        e(ibd2, ibs2) = 1;      
    else
        % ab ??
        e(ibd1, ibs0) = 0;
        e(ibd1, ibs1) = 0.5;
        e(ibd1, ibs2) = 0.5;

        e(ibd2, ibs0) = 0;
        e(ibd2, ibs1) = 0;
        e(ibd2, ibs2) = 1;      
    end
end



% both missing %checked with summation to 1
if( (f_geno(1) == 0 || f_geno(2) == 0) && (m_geno(1) == 0 || m_geno(2) == 0) ) 

    e(ibd1, ibs0) = 0;
    e(ibd1, ibs1) = 2 * p1 * p2;
    e(ibd1, ibs2) = p1*p1 + p2*p2;

    e(ibd2, ibs0) = 0;
    e(ibd2, ibs1) = 0;
    e(ibd2, ibs2) = 1;
 
end



end