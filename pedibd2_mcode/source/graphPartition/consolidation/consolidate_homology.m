function [chr_dist error] = consolidate_homology(spectrum, range, kinship2, debug, printid)

% let it consolidate if there is only one missing
chr_dist = [];
error = 0;

genotyped = range.family_range;
ngeno = length(genotyped);
if( ngeno <= 0 )
    error = 1;
    disp('no genotyped individuals');
    return;
end
nind = length(range.pedigree_range_full);
if( nind <= 0 || nind < ngeno || any(genotyped > nind) || any(genotyped < 1) )
    error = 1;
    disp('error in family structures');
    return;
end
family = range.structure;
[r, c] = size(family);
if( r <= 0 || r ~= nind || c < 12 )
    error = 1;
    disp('error in family structures');
    return;
end

[rows, cols] = size(spectrum);
if( rows <= 0 || cols <= 0 )
    error = 1;
    disp('error in chromosome spectrum');
    return;
end
if( cols ~= ngeno || rows ~= 2 * cols )
    error = 1;
    disp('error in chromosome spectrum');
    return;
end

[d1, d2, d3, d4] = size(kinship2);
if( d1 ~= ngeno || d2 ~= ngeno || d3 ~= 2 || d4 ~= 2 )
    error = 1;
    disp('error in kinship');
    return;
end

chr_affinity(1:2*ngeno,1:2*ngeno) = 0;
for i = 1:ngeno
    p_bit = (spectrum(2*i-1, 1:ngeno) ~= 0);
    m_bit = (spectrum(2*i, 1:ngeno) ~= 0);
    for j = 1:ngeno
        if( j  == i )
            chr_affinity(2*i-1,2*i-1) = ngeno;
            chr_affinity(2*i,2*i) = ngeno;
            continue;
        end        
        r_p_bit = (spectrum(j*2-1,1:ngeno) ~= 0);
        r_m_bit = (spectrum(j*2,1:ngeno) ~= 0);
        if( spectrum(i*2-1,j) ~= 0 )
            if( spectrum(j*2-1,i) ~= 0 )
                chr_affinity(i*2-1,j*2-1) = nnz(p_bit == r_p_bit);
            end
            if( spectrum(j*2,i) ~= 0 )
                chr_affinity(i*2-1,j*2) = nnz(p_bit == r_m_bit);
            end
        end
        if( spectrum(i*2,j) ~= 0 )
            if( spectrum(j*2-1,i) ~= 0 )
                chr_affinity(i*2,j*2-1) = nnz(m_bit == r_p_bit);
            end
            if( spectrum(j*2,i) ~= 0 )
                chr_affinity(i*2,j*2) = nnz(m_bit == r_m_bit);
            end         
        end
    end
end

homologous_chromosomes(1:2*ngeno,1:2*ngeno) = 0;
for i = 1:ngeno
    p_bit = (spectrum(2*i-1, 1:ngeno) ~= 0);
    m_bit = (spectrum(2*i, 1:ngeno) ~= 0);
    for j = 1:ngeno
        if( j  == i )
            homologous_chromosomes(2*i-1,2*i-1) = ngeno;
            homologous_chromosomes(2*i,2*i) = ngeno;
            continue;
        end
        r_p_bit = (spectrum(j*2-1, 1:ngeno) ~= 0);
        r_m_bit = (spectrum(j*2, 1:ngeno) ~= 0);
        if( spectrum(j*2-1,i) ~= 0 )
            if( spectrum(j*2,i) ~= 0 )
                if( spectrum(i*2-1,j) ~= 0 )
                    if( spectrum(i*2,j) ~= 0 )
                        if( nnz(p_bit == r_p_bit) + nnz(m_bit == r_m_bit) >= nnz(p_bit == r_m_bit) + nnz(m_bit == r_p_bit) )
                            homologous_chromosomes(2*i-1,j*2-1) = nnz(p_bit == r_p_bit);
                            homologous_chromosomes(2*i,j*2) = nnz(m_bit == r_m_bit);
                        else
                            homologous_chromosomes(2*i-1,j*2) = nnz(p_bit == r_m_bit);
                            homologous_chromosomes(2*i,j*2-1) = nnz(m_bit == r_p_bit);                    
                        end
                    else
                        if( nnz(p_bit == r_p_bit) >= nnz(p_bit == r_m_bit) )
                            homologous_chromosomes(2*i-1,j*2-1) = nnz(p_bit == r_p_bit);
                        else
                            homologous_chromosomes(2*i-1,j*2) = nnz(p_bit == r_m_bit);
                        end
                    end
                else
                    if( spectrum(i*2,j) ~= 0 )
                        if( nnz(m_bit == r_p_bit) >= nnz(m_bit == r_m_bit) )
                            homologous_chromosomes(2*i,j*2-1) = nnz(m_bit == r_p_bit);
                        else
                            homologous_chromosomes(2*i,j*2) = nnz(m_bit == r_m_bit);
                        end
                    else
                        % no homologous binding
                        if( debug == 1 )
                            disp('non-symmetric binding');
                        end
                    end
                end
            else
                if( spectrum(i*2-1,j) ~= 0 )
                    if( spectrum(i*2,j) ~= 0 )
                        if( nnz(p_bit == r_p_bit) >= nnz(m_bit == r_p_bit) )
                            homologous_chromosomes(2*i-1,j*2-1) = nnz(p_bit == r_p_bit);
                        else
                            homologous_chromosomes(2*i,j*2-1) = nnz(m_bit == r_p_bit);
                        end
                    else
                        homologous_chromosomes(2*i-1,j*2-1) = nnz(p_bit == r_p_bit);
                    end
                else
                    if( spectrum(i*2,j) ~= 0 )
                        homologous_chromosomes(2*i,j*2-1) = nnz(m_bit == r_p_bit);
                    else
                        % no homologous binding
                        if( debug == 1 )
                            disp('non-symmetric binding');
                        end
                    end
                end
            end
        else
            if( spectrum(j*2,i) ~= 0 )
                if( spectrum(i*2-1,j) ~= 0 )
                    if( spectrum(i*2,j) ~= 0 )
                        if( nnz(p_bit == r_m_bit) >= nnz(m_bit == r_m_bit) )
                            homologous_chromosomes(2*i-1,j*2) = nnz(p_bit == r_m_bit);
                        else
                            homologous_chromosomes(2*i,j*2) = nnz(m_bit == r_m_bit);
                        end
                    else
                        homologous_chromosomes(2*i-1,j*2) = nnz(p_bit == r_m_bit);
                    end
                else
                    if( spectrum(i*2,j) ~= 0 )
                        homologous_chromosomes(2*i,j*2) = nnz(m_bit == r_m_bit);
                    else
                        % no homologous binding
                        if( debug == 1 )
                            disp('non-symmetric binding');
                        end
                    end
                end
            else
                % no homologous binding
            end
        end
    end
end

% modify this to be used after pedigree merging
% heuristic merging identical homologous chromosomes
symHomo = (homologous_chromosomes + homologous_chromosomes')/2;
distHomo = (ngeno - symHomo);

T = 1:2*ngeno;
[T, error] = consolidate_threshold(T, ngeno - symHomo, 0.1);
if( error ~= 0 )
    disp('error in consolidating homologous chromosomes');
    return;
end


if( debug == 1 )
    figure;
    subplot(1,2,1);
    [X, perm] = sort(T);    
    spy(homologous_chromosomes(perm,perm), 1);
    for i = 1:2*ngeno
        for j = 1:2*ngeno
            if( homologous_chromosomes(perm(i),perm(j)) > 0 )
                text(i, j, num2str(ngeno - homologous_chromosomes(perm(i),perm(j))), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized', 'FontSize', 0.01);
            end
            if( chr_affinity(perm(i),perm(j)) > 0 )
%                 text(i, j, num2str(ngeno - chr_affinity(perm(i),perm(j))), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized', 'FontSize', 0.01);
            end            
        end
    end
    set(gca,'XTick',1:ngeno*2);
    set(gca,'XTickLabel',(num2cell(X)));
    set(gca,'YTick',1:ngeno*2);
    set(gca,'YTickLabel',(num2cell(perm)));    
end

if( debug == 1 )
    show_spec(1:ngeno*2+1,1:ngeno+1) = 0;
    show_spec(1:ngeno*2,1:ngeno) = spectrum;
    for i = 1:ngeno
        show_spec(i*2-1:i*2,ngeno+1) = printid(i);
    end
    show_spec(ngeno*2+1,1:ngeno) = printid(1:ngeno);
    show_spec(1:ngeno*2, ngeno+2) = T;
    show_spec
end



[T error] = consolidate_parentoffspring(T, family, kinship2, homologous_chromosomes, chr_affinity, debug);
if( error ~= 0 )
    disp('error in consilidating parent child');
    return;
end

[T error] = consolidate_siblings(T, family, kinship2, homologous_chromosomes, chr_affinity, debug);
if( error ~= 0 )
    disp('error in consolidating sibling chromosomes');
    return;
end

[T, error] = consolidate_threshold(T, ngeno - symHomo, 1);
if( error ~= 0 )
    disp('error in consolidating homologous chromosomes');
    return;
end

if( debug == 1 )
    subplot(1,2,2);
    [X, perm] = sort(T);
    spy(homologous_chromosomes(perm,perm), 1);
    for i = 1:2*ngeno
        for j = 1:2*ngeno
            if( homologous_chromosomes(perm(i),perm(j)) > 0 )
                text(i, j, num2str(ngeno - homologous_chromosomes(perm(i),perm(j))), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'FontUnits', 'normalized', 'FontSize', 0.01);
            end
        end
    end
    set(gca,'XTick',1:ngeno*2);
    set(gca,'XTickLabel',(num2cell(X)));
    set(gca,'YTick',1:ngeno*2);
    set(gca,'YTickLabel',(num2cell(perm)));    
end

if( debug == 1 )
    show_spec(1:ngeno*2+1,1:ngeno+1) = 0;
    show_spec(1:ngeno*2,1:ngeno) = spectrum;
    for i = 1:ngeno
        show_spec(i*2-1:i*2,ngeno+1) = printid(i);
    end
    show_spec(ngeno*2+1,1:ngeno) = printid(1:ngeno);
    show_spec(1:ngeno*2, ngeno+2) = T;
    show_spec
end

chr_dist = T;

end






























