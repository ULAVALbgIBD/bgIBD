function e = sibling_core_genotype(f_geno, m_geno, freq, emission_option)




a2g = emission_option.allele2genotype;
g2p = emission_option.genotype2pair;
a2i = emission_option.sibling_ibd_states;
nf = emission_option.nuclear_emission;

n_col = max(max(g2p));
n_row = max(max(max(max(a2i))));



p = freq(1:2);

e(1:n_row,1:n_col) = 0;

if( f_geno(1) ~= 0 && f_geno(2) ~= 0 && m_geno(1) ~= 0 && m_geno(2) ~= 0 )
    pair = code_pair([f_geno,m_geno], emission_option);
    e = nf{pair};
end
if( f_geno(1) == 0 || f_geno(2) == 0 )
    if( m_geno(1) ~= 0 && m_geno(2) ~= 0 )
        for i1 = 1:2
            for i2 = 1:2
                pair = code_pair([i1,i2,m_geno], emission_option);
                e = e + (p(i1)*p(i2)).*nf{pair};
            end
        end               
    end
end
if( f_geno(1) ~= 0 && f_geno(2) ~= 0 )
    if( m_geno(1) == 0 || m_geno(2) == 0 )
        for j1 = 1:2
            for j2 = 1:2
                pair = code_pair([f_geno,j1,j2], emission_option);
                e = e + (p(j1)*p(j2)).*nf{pair};
            end
        end
    end
end
if( f_geno(1) == 0 || f_geno(2) == 0 )
    if( m_geno(1) == 0 || m_geno(2) == 0 )
        for i1 = 1:2
            for i2 = 1:2
                for j1 = 1:2
                    for j2 = 1:2
                        pair = code_pair([i1,i2,j1,j2], emission_option);
                        e = e + (p(i1)*p(i2)*p(j1)*p(j2)).*nf{pair};
                    end
                end
            end
        end
    end
end

end


