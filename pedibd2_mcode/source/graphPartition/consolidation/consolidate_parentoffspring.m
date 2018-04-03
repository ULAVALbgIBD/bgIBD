function [T error] = consolidate_parentoffspring(T, family, kinship2, homologous_chromosomes, chr_affinity, debug)

error = 0;

[nind, c] = size(family);
if( nind <= 0 || c < 12 )
    error = 1;
    disp('error in family structures');
    return;
end

isgeno = (family(1:nind,7) == 1);
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

% fast access parent-child relationship
parent(1:nind,1:nind) = 0;
for i = 1:nind
    father = family(i,3);
    mother = family(i,4);
    if( father > 0 )
        parent(father, i) = -1;
    end
    if( mother > 0 )
        parent(mother, i) = 1;
    end
end
genoparent(1:ngeno,1:ngeno) = (parent(genotyped, genotyped) ~= 0);

founder(1:nind) = 0;
for i = 1:nind
    father = family(i,3);
    mother = family(i,4);
    if( father == 0 && mother == 0 )
        founder(i) = 1;
    end
end

% alleles cannot be combined, no inheritance paths
fconf(1:2*ngeno, 1:2*ngeno) = false;
for i = 1:ngeno
    for j = 1:ngeno
        if( i == j )
            if( kinship2(i,i,1,2) == 0 )
                fconf(i*2-1, i*2) = true;
                fconf(i*2, i*2-1) = true;
            end
        else
            if( all(all(kinship2(i,j,1:2,1:2) == 0)) )
                for p = 1:2
                    for q = 1:2
                        fconf(i*2-2+p, j*2-2+q) = true;
                    end
                end
            end            
        end
    end
end

% consolidate parent-child relationships
for i = 1:ngeno
    for j = i+1:ngeno
        if( genoparent(i,j) )
            g1(1:2*ngeno) = false;
            g1(2*i-1:2*i) = true;
            g2(1:2*ngeno) = false;
            g2(2*j-1:2*j) = true;
            [T error] = closest_link(T, g1, g2, homologous_chromosomes, fconf);
            if( error ~= 0 )
                disp('cannot combine homologous chromosomes, parent-child');
                return;
            end               
        end
    end
end

[ngs error] = genotyped_envelope(family);
if( error ~= 0 )
    disp('error in generating genotyped envelope');
    return;
end
[d1 d2] = size(ngs);
if( d1 ~= nind || d1 ~= d2 )
    error = 1;
    disp('error in generating genotyped envelope');
    return;
end


% map all individuals to genotyped individuals, both alleles
contractMap(1:nind) = 0;
contractMap(genotyped) = 1:ngeno;

ReP(1:ngeno,1:2*ngeno) = false;
ReM(1:ngeno,1:2*ngeno) = false;

for i = 1:nind
    if( isgeno(i) )
        father = family(i,3);
        mother = family(i,4);
        if( father ~= 0 )
            pshield = ngs(father, 1:nind);
            if( any(pshield) && all(isgeno(pshield)) )
                ReP(contractMap(i), contractMap(pshield)*2-1) = true;
                ReP(contractMap(i), contractMap(pshield)*2) = true;
            end
        end
        if( mother ~= 0 )
            mshield = ngs(mother, 1:nind);
            if( any(mshield) && all(isgeno(mshield)) )
                ReM(contractMap(i), contractMap(mshield)*2-1) = true;
                ReM(contractMap(i), contractMap(mshield)*2) = true;
            end
        end        
    end
end

% consolidating with ancestral parental envelopes, closed shield
for i = 1:ngeno
    g1(1:2*ngeno) = false;
    g1(i*2-1:i*2) = true;
    g2 = ReP(i,1:2*ngeno);
    if( any(g2) )
        [T error] = closest_link(T, g1, g2, homologous_chromosomes, fconf);
        if( error ~= 0 )
            disp('cannot combine homologous chromosomes');
            return;
        end        
    end
    g2 = ReM(i,1:2*ngeno);
    if( any(g2) )
        [T error] = closest_link(T, g1, g2, homologous_chromosomes, fconf);
        if( error ~= 0 )
            disp('cannot combine homologous chromosomes');
            return;
        end        
    end    
end




end
















