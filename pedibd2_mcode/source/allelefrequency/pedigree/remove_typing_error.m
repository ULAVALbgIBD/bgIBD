function [family_genotype mendel_error error] = remove_typing_error(genotype, family)

error = 0;
family_genotype = [];
mendel_error = [];
global debug_mode;

[nind c] = size( family );
if( nind <= 0 || c < 12 )
    error = 1;
    disp('error in family structures');
    return;
end

[nmarkers, d2, d3] = size(genotype);
if( d2 <= 0 || d2 ~= nind || nmarkers <= 0 || d3 ~= 2 )
    error = 1;
    disp('error in genotype data')
    return;
end

family_genotype = genotype;


error_dist = zeros(nmarkers, 1);
count_gt = 0;
count_all = 0;
for i = 1:nind
    father = family(i,3);
    mother = family(i,4);
    if( father == 0 && mother == 0 )
        continue;
    end
    f1 = zeros(nmarkers, 1);
    f2 = zeros(nmarkers, 1);
    m1 = zeros(nmarkers, 1);
    m2 = zeros(nmarkers, 1);
    if( father ~= 0 )
        f1(1:nmarkers) = genotype(1:nmarkers, father, 1);
        f2(1:nmarkers) = genotype(1:nmarkers, father, 2);
    end
    if( mother ~= 0 )
        m1(1:nmarkers) = genotype(1:nmarkers, mother, 1);
        m2(1:nmarkers) = genotype(1:nmarkers, mother, 2);
    end
    g1 = genotype(1:nmarkers, i, 1);
    g2 = genotype(1:nmarkers, i, 2);
    if( all(g1 == 0) && all(g2 == 0) )
        continue;
    end
    
    f1good = (f1 == 0) | (f2 == 0) | (g1 == 0) | (f1 == g1) | (f2 == g1);
    f2good = (f1 == 0) | (f2 == 0) | (g2 == 0) | (f1 == g2) | (f2 == g2);
    m1good = (m1 == 0) | (m2 == 0) | (g2 == 0) | (m1 == g2) | (m2 == g2);
    m2good = (m1 == 0) | (m2 == 0) | (g1 == 0) | (m1 == g1) | (m2 == g1);
    
    phase1good = f1good & m1good;
    phase2good = f2good & m2good;
    good = phase1good | phase2good;
    
    valid = ( g1 ~= 0 & g2 ~= 0 ) & ( ( f1 ~= 0 & f2 ~= 0 ) | ( m1 ~= 0 & m2 ~= 0 ) );
    
    count_all = count_all + nnz(valid);
    count_gt = count_gt + nnz(~good);
    error_dist(~good) = error_dist(~good) + 1;
    
end

family_genotype(error_dist > 0, 1:nind, 1:2) = 0;

mendel_error = error_dist > 0;

if( debug_mode == 1 )
    if( count_all > 0 )
        disp(['family ', num2str(family(1,12)), ': Mendelian inconsistency ', num2str(count_gt*100/count_all, '%.2f'), '%']);
        disp(['family ', num2str(family(1,12)), ': Loci removed ', num2str(nnz(error_dist > 0)*100/nmarkers, '%.2f'), '%']);
    end
end

if( count_all > 0 )
    if( count_gt/count_all > 0.05 )
        disp(['family ', num2str(family(1,12)), ': Mendelian inconsistency ', num2str(count_gt*100/count_all, '%.2f'), '%']);
        disp(['family ', num2str(family(1,12)), ': possible relationship mis-specification']);
        error = 1;
        return;
    end
    if( count_gt/count_all > 0.02 )
        disp(['family ', num2str(family(1,12)), ': Mendelian inconsistency ', num2str(count_gt*100/count_all, '%.2f'), '%']);
        disp(['family ', num2str(family(1,12)), ': Data quality check fails']);
        error = 1;
        return;
    end    
end

end









