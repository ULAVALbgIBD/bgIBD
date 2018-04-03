function [ assignment, consistency, typing_error, conflict, error ] = assign_locus_flat( constraint, genotype )

assignment = [];
consistency = 0;
typing_error = 0;
error = 0;

[r1 c1] = size(constraint);
if( r1 <= 0 || c1 ~= 2 )
    error = 1;
    disp('error in global IBD');
    return;
end
[r2 c2] = size(genotype);
if( r2 <= 0 || c2 ~= 2 || r1 ~= r2 )
    error = 1;
    disp('error in genotype data');
    return;
end
nind = r1;

rep(1:2*nind) = 0;
allele(1:2*nind) = 0;

for i = 1:nind
    rep(i*2-1:i*2) = constraint(i,1:2);
end


% build homozygous
consistency = nind;
processed(1:nind) = false;
for i = 1:nind
    if( ~all(genotype(i,1:2) ~= 0) )
        continue;
    end
    p = genotype(i,1);
    m = genotype(i,2);
    if( p ~= m )
        continue;
    end
    group1 = (rep == rep(i*2-1));
    group2 = (rep == rep(i*2));    
    if( all(allele(group1) == p) || all(allele(group1) == 0) )
        allele(group1) = p;
    else
        consistency = i-1;
        break;
    end
    if( all(allele(group2) == m) || all(allele(group2) == 0) )
        allele(group2) = m;
    else
        consistency = i-1;
        break;
    end
    processed(i) = true;
end


% build heterozygous, one side fixed
changed = true;
while( consistency == nind && changed && any(~processed) )
    changed = false;
    for i = 1:nind
        if( processed(i) )
            continue;
        end
        if( ~all(genotype(i,1:2) ~= 0) )
            continue;
        end
        p = genotype(i,1);
        m = genotype(i,2);
        if( p == m )
            continue;
        end
        group1 = (rep == rep(i*2-1));
        group2 = (rep == rep(i*2));            
        if( all(allele(group1) == p) && all(allele(group2) == 0) )
            allele(group1) = p;
            allele(group2) = m;
            changed = true;
            processed(i) = true;
            continue;
        end
        if( all(allele(group1) == m) && all(allele(group2) == 0) )
            allele(group1) = m;
            allele(group2) = p;
            changed = true;
            processed(i) = true;
            continue;
        end
        if( all(allele(group1) == 0) && all(allele(group2) == p) )
            allele(group1) = m;
            allele(group2) = p;
            changed = true;
            processed(i) = true;
            continue;
        end
        if( all(allele(group1) == 0) && all(allele(group2) == m) )
            allele(group1) = p;
            allele(group2) = m;
            changed = true;
            processed(i) = true;
            continue;
        end
        if( all(allele(group1) == 0) && all(allele(group2) == 0) )
            continue;
        end
        if( all(allele(group1) == p) && all(allele(group2) == m) )
            continue;
        end
        if( all(allele(group1) == m) && all(allele(group2) == p) )
            continue;
        end
        consistency = i-1;
        break;
    end
end

% build heterozygous both sides unfixed
if( consistency == nind && any(~processed) )
    repibs = rep;
    offset(1:2*nind) = false;
    for i = 1:nind
        if( processed(i) )
            continue;
        end
        if( ~all(genotype(i,1:2) ~= 0) )
            continue;
        end
        p = genotype(i,1);
        m = genotype(i,2);
        if( p == m )
            continue;
        end
        group1 = (repibs == repibs(i*2-1));
        group2 = (repibs == repibs(i*2));  
        if( repibs(i*2-1) ~= repibs(i*2) )
            repibs(group2) = repibs(i*2-1);
            if( offset(i*2-1) == offset(i*2) )
                offset(group2) = ~offset(group2);
            end
        else
            if( offset(i*2-1) == offset(i*2) )
                consistency = i-1;
                break;
            end
        end
    end
end

if( consistency < nind )
    allele_vote(1:2*nind,1:2) = 0;
    for i = 1:nind
        if( ~all(genotype(i,1:2) ~= 0) )
            continue;
        end
        p = genotype(i,1);
        m = genotype(i,2);
        if( p ~= m )
            continue;
        end
        group1 = (rep == rep(i*2-1));
        group2 = (rep == rep(i*2));  
        allele_vote(group1, p) = allele_vote(group1, p) + 1;
        allele_vote(group2, m) = allele_vote(group2, m) + 1;
    end
    max_vote(1:2*nind) = 0;
    full_confidence(1:2*nind) = false;
    for i = 1:2*nind
        if( allele_vote(i,1) > allele_vote(i,2) )
            max_vote(i) = 1;
            if( allele_vote(i,2) == 0 && allele_vote(i,1) > 1 )
                full_confidence(i) = true;
            end
        end
        if( allele_vote(i,2) > allele_vote(i,1) )
            max_vote(i) = 2;
            if( allele_vote(i,1) == 0 && allele_vote(i,2) > 1 )
                full_confidence(i) = true;
            end
        end
    end
end

conflict(1:nind) = false;
assignment(1:nind,1:2) = zeros(nind,2);

if( consistency == nind )
    for i = 1:nind
        assignment(i,1:2) = allele(i*2-1:i*2);
    end
else
    % non-missing individual get assigned
    % even if locus is pronounced inconsistent
    for i = 1:nind
        if( all(genotype(i,1:2) ~= 0) )
            p = genotype(i,1);
            m = genotype(i,2);
            if( p == m )
                assignment(i,1) = p;
                assignment(i,2) = m;
            else
                phase1 = allele_vote(i*2-1,p) + allele_vote(i*2,m);
                phase2 = allele_vote(i*2-1,m) + allele_vote(i*2,p);
                if( phase1 >= phase2 )
                    assignment(i,1) = p;
                    assignment(i,2) = m;
                end
                if( phase2 > phase1 )
                    assignment(i,1) = m;
                    assignment(i,2) = p;
                end
            end
            if( all( assignment(i,1:2) ~= 0 ) )
                if( all(max_vote(2*i-1:2*i) ~= 0) )
                    if( ~all(assignment(i,1:2) == max_vote(2*i-1:2*i)) && ~all(assignment(i,2:-1:1) == max_vote(2*i-1:2*i)) )
                        typing_error = typing_error + 1;
                        conflict(i) = true;
                    end
                end
            end
        else
            if( all(full_confidence(2*i-1:2*i)) )
                % do not assign for erroneous loci
%                 assignment(i,1:2) = max_vote(2*i-1:2*i);
            end
        end
    end
end



end



















