function [T error] = consolidate_siblings(T, family, kinship2, homologous_chromosomes, chr_affinity, debug)

error = 0;

[nIND, c] = size(family);
if( nIND <= 0 || c < 12 )
    error = 1;
    disp('error in family structures');
    return;
end

isgeno = (family(1:nIND,7) == 1);
genotyped = find(isgeno);
ngeno = nnz(isgeno);

if( ngeno < 2 )
    return;
end

[d1 d2 d3 d4] = size(kinship2);
if( d1 ~= ngeno || d2 ~= d1 || d3 ~= 2 || d4 ~= 2 )
    error = 1;
    disp('error in kinship');
    reutrn;
end

[d1 d2] = size(homologous_chromosomes);
if( d1 ~= 2 * ngeno || d2 ~= d1 )
    error = 1;
    disp('error in chromosome homology');
    return;
end

[d1 d2] = size(chr_affinity);
if( d1 ~= 2 * ngeno || d2 ~= d1 )
    error = 1;
    disp('error in chromosome affinity');
    return;
end

% alleles cannot be combined, no inheritance paths
fconf(1:2*ngeno, 1:2*ngeno) = 0;
for i = 1:ngeno
    for j = 1:ngeno
        if( i == j )
            if( kinship2(i,i,1,2) == 0 )
                fconf(i*2-1,i*2) = 1;
                fconf(i*2, i*2-1) = 1;
            end
        else
            if( all(all(kinship2(i,j,1:2,1:2) == 0)) )
                for p = 1:2
                    for q = 1:2
                        fconf(i*2-2+p,j*2-2+q) = 1;
                    end
                end
            end            
        end
    end
end


full_sib(1:nIND, 1:nIND, 1:nIND) = 0;

for i = 1:nIND
    father = family(i,3);
    mother = family(i,4);
    if( family(i,7) > 0 && father ~= 0 && mother ~= 0 )
        full_sib(father, mother, i) = 1;
    end
end

is_sib(1:nIND, 1:nIND) = false;
for i = 1:nIND
    for j = 1:nIND
        sibship = (full_sib(i,j,1:nIND) > 0);
        if( any(sibship) )
            is_sib(sibship, sibship) = true;
        end
    end
end

% consolidate siblings, if more than 4 alleles
isgeno(1:nIND) = ( family(1:nIND,7) == 1 );
for i = 1:nIND
    for j = 1:nIND
        if( ~any(full_sib(i,j,1:nIND) > 0) )
            continue;
        end
        index = (full_sib(i,j,genotyped) > 0);
        if( ~any(index) )
            error = 1;
            disp('error in family indexing');
            return;
        end
        while(true)
            sibs(1:ngeno*2) = 0;
            for k = 1:ngeno
                if( index(k) > 0 )
                    sibs(k*2-1) = T(k*2-1);
                    sibs(k*2) = T(k*2);
                end
            end
            if( ~any(sibs > 0) )
                error = 1;
                disp('error in family indexing');
                return;
            end
            if( length(unique(sibs(sibs>0))) <= 4 )
                break;
            end
            if( nnz(index) < 3 )
                error = 1;
                disp('error in family indexing');
                return;
            end
            pmax = 0;
            qmax = 0;
            max_coherence = 0;
            paff = 0;
            qaff = 0;
            max_aff = 0;
            Tparent(1:4) = 0;
            if( isgeno(i) )
                ii = find(genotyped == i);
                if( isempty(ii) || length(ii) > 1 )
                    error = 1;
                    disp('error in family indexing');
                    return;
                end
                Tparent(1:2) = T(2*ii-1:2*ii);
            end
            if( isgeno(j) )
                jj = find(genotyped == j);
                if( isempty(jj) || length(jj) > 1 )
                    error = 1;
                    disp('error in family indexing');
                    return;
                end
                Tparent(3:4) = T(2*jj-1:2*jj);
            end
            if( all(Tparent ~= 0) )
                break;
            end
            for p = 1:ngeno
                if( index(p) <= 0 )
                    continue;
                end
                for q = p+1:ngeno
                    if( index(q) <= 0 )
                        continue;
                    end
                    for r = 1:2
                        if( any(Tparent == T(p*2+1-r)) )
                            continue;
                        end
                        for s = 1:2
                            if( any(Tparent == T(q*2+1-s)) )
                                continue;
                            end
                            if( T(p*2+1-r) == T(q*2+1-s) )
                                continue;
                            end
                            temp1 = (T == T(p*2+1-r));
                            temp2 = (T == T(q*2+1-s));
                            if( any(any(fconf(temp1, temp2) == 1)) )
                                continue;
                            end                            
                            temp3 = mean(mean(homologous_chromosomes(temp1,temp2)));
                            temp4 = mean(mean(chr_affinity(temp1,temp2)));
                            if( temp3 > max_coherence )
                                max_coherence = temp3;
                                pmax = p*2+1-r;
                                qmax = q*2+1-s;
                            end
                            if( temp4 > max_aff )
                                max_aff = temp4;
                                paff = p*2+1-r;
                                qaff = q*2+1-s;
                            end
                        end
                    end                
                end
            end
            if( pmax == 0 || qmax == 0 )
                pmax = paff;
                qmax = qaff;
            end
            if( pmax ~= 0 && qmax ~= 0 )
                c1 = T(pmax);
                c2 = T(qmax);
                if( c1 < c2 )
                    T( T == c2 ) = c1;
                else
                    T( T == c1 ) = c2;              
                end
            else
                if( debug == 1 )
                    disp('cannot consolidate siblings');
                end
                break;
            end
        end
    end
end

% combine sibling sharing, high-confidence sharing
for i = 1:nIND
    for j = 1:nIND
        if( i == j )
            continue;
        end
        if( ~is_sib(i,j) )
            continue;
        end
        if( ~isgeno(i) || ~isgeno(j) )
            continue;
        end
        fa = family(i,3);
        mo = family(i,4);
        if( fa == 0 || fa ~= family(j,3) )
            error = 1;
            disp('error in family indexing');
            return;
        end
        if( mo == 0 || mo ~= family(j,4) )
            error = 1;
            disp('error in family indexing');
            return;
        end
        Tparent(1:4) = 0;
        ifa = find(genotyped == fa);
        if( length(ifa) == 1 )
            Tparent(1:2) = T(2*ifa-1:2*ifa);
        end
        ima = find(genotyped == mo);
        if( length(ima) == 1 )
            Tparent(3:4) = T(2*ima-1:2*ima);
        end
        if( all(Tparent ~= 0) )
            break;
        end
        p = find(genotyped == i);
        q = find(genotyped == j);
        if( length(p) ~= 1 || length(q) ~= 1 )
            error = 1;
            disp('error in family indexing');
            return;
        end
        pmax = 0;
        qmax = 0;
        max_coherence = 0;
        for r = 1:2
            if( any(Tparent == T(p*2+1-r)) )
                continue;
            end
            for s = 1:2
                if( any(Tparent == T(q*2+1-s)) )
                    continue;
                end
                if( T(p*2+1-r) == T(q*2+1-s) )
                    continue;
                end
                if( chr_affinity(p*2+1-s, q*2+1-s) <= 0 )
                    continue;
                end
                temp1 = (T == T(p*2+1-r));
                temp2 = (T == T(q*2+1-s));
                if( any(any(fconf(temp1, temp2) == 1)) )
                    continue;
                end                            
                temp3 = mean(mean(homologous_chromosomes(temp1,temp2)));
                if( temp3 > max_coherence )
                    max_coherence = temp3;
                    pmax = p*2+1-r;
                    qmax = q*2+1-s;
                end
            end
        end                
        if( pmax ~= 0 && qmax ~= 0 && max_coherence >= ngeno/2)
            c1 = T(pmax);
            c2 = T(qmax);
            if( c1 < c2 )
                T( T == c2 ) = c1;
            else
                T( T == c1 ) = c2;              
            end
        end
    end
end



end
