function e = sibling_core(f_geno, m_geno, freq)

ibd0 = 1;
ibd_l = 2;
ibd_r = 3;
ibd_lr = 4;

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
            e(ibd0, ibs0) = 0;
            e(ibd0, ibs1) = 0;
            e(ibd0, ibs2) = 1;

            e(ibd_l, ibs0) = 0;
            e(ibd_l, ibs1) = 0;
            e(ibd_l, ibs2) = 1;

            e(ibd_r, ibs0) = 0;
            e(ibd_r, ibs1) = 0;
            e(ibd_r, ibs2) = 1;

            e(ibd_lr, ibs0) = 0;
            e(ibd_lr, ibs1) = 0;
            e(ibd_lr, ibs2) = 1;
        else
            % aa, cd
            e(ibd0, ibs0) = 0;
            e(ibd0, ibs1) = 1;
            e(ibd0, ibs2) = 0;

            e(ibd_l, ibs0) = 0;
            e(ibd_l, ibs1) = 1;
            e(ibd_l, ibs2) = 0;

            e(ibd_r, ibs0) = 0;
            e(ibd_r, ibs1) = 0;
            e(ibd_r, ibs2) = 1;

            e(ibd_lr, ibs0) = 0;
            e(ibd_lr, ibs1) = 0;
            e(ibd_lr, ibs2) = 1;        
        end
    else
        if( m_geno(1) == m_geno(2) )
            % ab, cc
            e(ibd0, ibs0) = 0;
            e(ibd0, ibs1) = 1;
            e(ibd0, ibs2) = 0;

            e(ibd_l, ibs0) = 0;
            e(ibd_l, ibs1) = 0;
            e(ibd_l, ibs2) = 1;

            e(ibd_r, ibs0) = 0;
            e(ibd_r, ibs1) = 1;
            e(ibd_r, ibs2) = 0;

            e(ibd_lr, ibs0) = 0;
            e(ibd_lr, ibs1) = 0;
            e(ibd_lr, ibs2) = 1;
        else
            % ab, cd
            % critical condition
            e(ibd0, ibs0) = 0.5;
            e(ibd0, ibs1) = 0;
            e(ibd0, ibs2) = 0.5;

            e(ibd_l, ibs0) = 0;
            e(ibd_l, ibs1) = 1;
            e(ibd_l, ibs2) = 0;

            e(ibd_r, ibs0) = 0;
            e(ibd_r, ibs1) = 1;
            e(ibd_r, ibs2) = 0;

            e(ibd_lr, ibs0) = 0;
            e(ibd_lr, ibs1) = 0;
            e(ibd_lr, ibs2) = 1;        
        end    
    end
end

% father missing
if( (f_geno(1) == 0 || f_geno(2) == 0) && m_geno(1) ~= 0 && m_geno(2) ~= 0 )
    if( m_geno(1) == m_geno(2) )
        e(ibd0, ibs0) = 0;
        e(ibd0, ibs1) = 2*p1*p2;
        e(ibd0, ibs2) = p1*p1 + p2*p2;

        e(ibd_l, ibs0) = 0;
        e(ibd_l, ibs1) = 0;
        e(ibd_l, ibs2) = 1;

        e(ibd_r, ibs0) = 0;
        e(ibd_r, ibs1) = 2*p1*p2;
        e(ibd_r, ibs2) = p1*p1 + p2*p2;

        e(ibd_lr, ibs0) = 0;
        e(ibd_lr, ibs1) = 0;
        e(ibd_lr, ibs2) = 1;        
    else
        e(ibd0, ibs0) = p1*p2;
        e(ibd0, ibs1) = p1*p1 + p2*p2;
        e(ibd0, ibs2) = p1*p2;

        e(ibd_l, ibs0) = 0;
        e(ibd_l, ibs1) = 1;
        e(ibd_l, ibs2) = 0;

        e(ibd_r, ibs0) = 0;
        e(ibd_r, ibs1) = 2*p1*p2;
        e(ibd_r, ibs2) = p1*p1 + p2*p2;

        e(ibd_lr, ibs0) = 0;
        e(ibd_lr, ibs1) = 0;
        e(ibd_lr, ibs2) = 1;        
    end
end

% mother missing
if( f_geno(1) ~= 0 && f_geno(2) ~= 0 && ( m_geno(1) == 0 || m_geno(2) == 0 ) )
    if( f_geno(1) == f_geno(2) )
        e(ibd0, ibs0) = 0;
        e(ibd0, ibs1) = 2*p1*p2;
        e(ibd0, ibs2) = p1*p1 + p2*p2;

        e(ibd_l, ibs0) = 0;
        e(ibd_l, ibs1) = 2*p1*p2;
        e(ibd_l, ibs2) = p1*p1 + p2*p2;

        e(ibd_r, ibs0) = 0;
        e(ibd_r, ibs1) = 0;
        e(ibd_r, ibs2) = 1;

        e(ibd_lr, ibs0) = 0;
        e(ibd_lr, ibs1) = 0;
        e(ibd_lr, ibs2) = 1;          
    else
        e(ibd0, ibs0) = p1*p2;
        e(ibd0, ibs1) = p1*p1 + p2*p2;
        e(ibd0, ibs2) = p1*p2;

        e(ibd_l, ibs0) = 0;
        e(ibd_l, ibs1) = 2*p1*p2;
        e(ibd_l, ibs2) = p1*p1 + p2*p2;

        e(ibd_r, ibs0) = 0;
        e(ibd_r, ibs1) = 1;
        e(ibd_r, ibs2) = 0;

        e(ibd_lr, ibs0) = 0;
        e(ibd_lr, ibs1) = 0;
        e(ibd_lr, ibs2) = 1;          
    end
end



% both missing %checked with summation to 1
if( (f_geno(1) == 0 || f_geno(2) == 0) && (m_geno(1) == 0 || m_geno(2) == 0) ) 
    e(ibd0, ibs0) = (p1*p1*p2*p2) + (p2*p2*p1*p1);
    e(ibd0, ibs1) = (2*p1*p2)*(p1*p1+p2*p2) + (p1*p1+p2*p2)*(2*p1*p2);
    e(ibd0, ibs2) = (2*p1*p2)^2 + (p1*p1)^2 + (p2*p2)^2;

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