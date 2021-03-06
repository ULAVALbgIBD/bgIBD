function e = sibling_side_paired(f_geno, m_geno, freq)

% paired distant sharing

ibd_l = 1;
ibd_r = 2;
ibd_lr = 3;

ibs0 = 1;
ibs1 = 2;
ibs2 = 3;
ibs3 = 4; %missing

p1 = freq(1);
p2 = freq(2);

e(1:4,1:3) = 0;

if( f_geno(1) ~= 0 && f_geno(2) ~= 0 && m_geno(1) ~= 0 && m_geno(2) ~= 0 )
    if( f_geno(1) == f_geno(2) )
        if( m_geno(1) == m_geno(2) )
            % aa, cc
            e(ibd_l, ibs0) = 0;
            e(ibd_l, ibs1) = 0;
            e(ibd_l, ibs2) = 0;
            
            e(ibd_l, ibs0) = 0;
            e(ibd_l, ibs1) = 0;
            e(ibd_l, ibs2) = 0;
            
            e(ibd_lr, ibs0) = 0;
            e(ibd_lr, ibs1) = 0;
            e(ibd_lr, ibs2) = 1;            
        else
            % aa, cd
            e(ibd_l, ibs0) = 0;
            e(ibd_l, ibs1) = 1;
            e(ibd_l, ibs2) = 0; 
            
            e(ibd_r, ibs0) = 0;
            e(ibd_r, ibs1) = 0;
            e(ibd_r, ibs2) = 0;
            
            e(ibd_lr, ibs0) = 0;
            e(ibd_lr, ibs1) = 0;
            e(ibd_lr, ibs2) = 0;            
        end
    else
        if( m_geno(1) == m_geno(2) )
            % ab, cc

            e(ibd_l, ibs0) = 0;
            e(ibd_l, ibs1) = 0;
            e(ibd_l, ibs2) = 0;

            e(ibd_r, ibs0) = 0;
            e(ibd_r, ibs1) = 1;
            e(ibd_r, ibs2) = 0;

            e(ibd_lr, ibs0) = 0;
            e(ibd_lr, ibs1) = 0;
            e(ibd_lr, ibs2) = 0;
        else
            % ab, cd
            % critical condition

            e(ibd_l, ibs0) = 0;
            e(ibd_l, ibs1) = 0;
            e(ibd_l, ibs2) = 0;

            e(ibd_r, ibs0) = 0;
            e(ibd_r, ibs1) = 0;
            e(ibd_r, ibs2) = 0;

            e(ibd_lr, ibs0) = 0;
            e(ibd_lr, ibs1) = 0;
            e(ibd_lr, ibs2) = 0;        
        end    
    end
end

% father missing
if( (f_geno(1) == 0 || f_geno(2) == 0) && m_geno(1) ~= 0 && m_geno(2) ~= 0 )
    if( m_geno(1) == m_geno(2) )

        e(ibd_l, ibs0) = 0;
        e(ibd_l, ibs1) = 0;
        e(ibd_l, ibs2) = 0;

        e(ibd_r, ibs0) = 0;
        e(ibd_r, ibs1) = 2*p1*p2;
        e(ibd_r, ibs2) = p1*p1 + p2*p2;

        e(ibd_lr, ibs0) = 0;
        e(ibd_lr, ibs1) = 0;
        e(ibd_lr, ibs2) = 1;        
    else
        e(ibd_l, ibs0) = 0;
        e(ibd_l, ibs1) = 0;
        e(ibd_l, ibs2) = 1;

        e(ibd_r, ibs0) = 0;
        e(ibd_r, ibs1) = 0;
        e(ibd_r, ibs2) = 0;

        e(ibd_lr, ibs0) = 0;
        e(ibd_lr, ibs1) = 0;
        e(ibd_lr, ibs2) = 0;        
    end
end

% mother missing
if( f_geno(1) ~= 0 && f_geno(2) ~= 0 && ( m_geno(1) == 0 || m_geno(2) == 0 ) )
    if( f_geno(1) == f_geno(2) )
        e(ibd_l, ibs0) = 0;
        e(ibd_l, ibs1) = 2*p1*p2;
        e(ibd_l, ibs2) = p1*p1 + p2*p2;

        e(ibd_r, ibs0) = 0;
        e(ibd_r, ibs1) = 0;
        e(ibd_r, ibs2) = 0;

        e(ibd_lr, ibs0) = 0;
        e(ibd_lr, ibs1) = 0;
        e(ibd_lr, ibs2) = 1;          
    else
        e(ibd_l, ibs0) = 0;
        e(ibd_l, ibs1) = 0;
        e(ibd_l, ibs2) = 0;

        e(ibd_r, ibs0) = 0;
        e(ibd_r, ibs1) = 1;
        e(ibd_r, ibs2) = 0;

        e(ibd_lr, ibs0) = 0;
        e(ibd_lr, ibs1) = 0;
        e(ibd_lr, ibs2) = 0;          
    end
end



% both missing %checked with summation to 1
if( (f_geno(1) == 0 || f_geno(2) == 0) && (m_geno(1) == 0 || m_geno(2) == 0) ) 

    e(ibd_l, ibs0) = 0;
    e(ibd_l, ibs1) = 2 * p1 * p2;
    e(ibd_l, ibs2) = p1*p1 + p2*p2;

    e(ibd_r, ibs0) = 0;
    e(ibd_r, ibs1) = 2 * p1 * p2;
    e(ibd_r, ibs2) = p1*p1 + p2*p2;

    e(ibd_lr, ibs0) = 0;
    e(ibd_lr, ibs1) = 0;
    e(ibd_lr, ibs2) = 1;  
end



end