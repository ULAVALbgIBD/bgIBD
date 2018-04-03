function [epl epl_log error] = generate_emission_all_loci(all_af, ms, te, markerlist, father_geno, mother_geno, sibpair, emission_option)


error = 0;
epl = [];

[nloc, fields] = size(markerlist);
if( nloc <= 0 || fields ~= 2 )
    disp('error in marker list');
    return;
end

[r, c] = size(all_af);
if( r <= 0 || r ~= nloc || c ~= 3 )
    disp('error in allele frequency');
    return;
end

[r, c] = size(mother_geno);
if( r <= 0 || r ~= nloc || c ~= 2 )
    disp('error in genotype data');
    return;
end

[r, c] = size(father_geno);
if( r <= 0 || r ~= nloc || c ~= 2 )
    disp('error in genotype data');
    return;
end

if( sibpair == 1 )
    fg = father_geno(1,1:2);
    mg = mother_geno(1,1:2);    
    [epl_loc error] = sibling_generate_emission(fg, mg, all_af(1,1:3), ms, te, emission_option);
else
    [epl_loc error] = generate_emission(all_af(1,1:3), ms, te, emission_option);
end
[n_ibd n_ibs] = size(epl_loc);
if( n_ibd <= 0 || n_ibs <= 0 )
    error = 1;
    disp('error in generating emission states');
    return;
end

epl = zeros(nloc, n_ibd, n_ibs);
epl_log = zeros(nloc, n_ibd, n_ibs);

res = 100;
flag(1:res+1) = false;
sib_flag(1:3*3*3*3*(res+1)) = false;
emission_static(1:res+1, 1:n_ibd, 1:n_ibs) = 0;
emission_sib_static(1:3*3*3*3*(res+1),1:n_ibd,1:n_ibs) = 0;
minor_allele(1:nloc,1) = round(res * min(all_af(1:nloc,1), all_af(1:nloc,2))) + 1;
if( sibpair == 1 )
    map(1:nloc) = (res+1)*3*3*3*father_geno(1:nloc,1) + (res+1)*3*3*father_geno(1:nloc,2) + (res+1)*3*mother_geno(1:nloc,1) + (res+1)*mother_geno(1:nloc,2) + minor_allele(1:nloc);
end

% each allele has 0,1,2 three states, plus res+1 state for general al-feq
% 3*3*3*3*(res+1)

for i = 1:nloc
    if( sibpair == 1 )
        fg = father_geno(i,1:2);
        mg = mother_geno(i,1:2);        
        if( ~sib_flag(map(i)) )
            [epl_loc error] = sibling_generate_emission(fg, mg, all_af(i,1:3), ms, te, emission_option);
            if( error ~= 0 )
                disp('error in generating emission states');
                return;
            end
            [nrows ncols] = size(epl_loc);
            if( nrows ~= n_ibd || ncols ~= n_ibs )
                error = 1;
                disp('error in generating emission states');
                return;
            end
            emission_sib_static(map(i), 1:n_ibd, 1:n_ibs) = epl_loc(1:n_ibd, 1:n_ibs);
            sib_flag(map(i)) = true;
        end
    else        
        if( ~flag(minor_allele(i)) )
            [epl_loc error] = generate_emission(all_af(i,1:3), ms, te, emission_option);
            if( error ~= 0 )
                disp('error in generating emission states');
                return;
            end
            [nrows ncols] = size(epl_loc);
            if( nrows ~= n_ibd || ncols ~= n_ibs )
                error = 1;
                disp('error in generating emission states');
                return;
            end            
            emission_static(minor_allele(i), 1:n_ibd, 1:n_ibs) = epl_loc(1:n_ibd, 1:n_ibs);
            flag(minor_allele(i)) = true;
        end
    end
end

if( sibpair == 1 )
    temp = emission_sib_static;
    valid = emission_sib_static > 0;
    temp(~valid) = -30;
    temp(valid) = reallog(emission_sib_static(valid));
    epl(1:nloc, 1:n_ibd, 1:n_ibs) = emission_sib_static(map(1:nloc), 1:n_ibd, 1:n_ibs);
    epl_log(1:nloc, 1:n_ibd, 1:n_ibs) = temp(map(1:nloc), 1:n_ibd, 1:n_ibs);
else
    temp = emission_static;
    valid = emission_static > 0;
    temp(~valid) = -30;
    temp(valid) = reallog(emission_static(valid));
    epl(1:nloc, 1:n_ibd, 1:n_ibs) = emission_static(minor_allele(1:nloc), 1:n_ibd, 1:n_ibs);
    epl_log(1:nloc, 1:n_ibd, 1:n_ibs) = temp(minor_allele(1:nloc), 1:n_ibd, 1:n_ibs);
end


end


