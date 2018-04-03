function [assignment phased error] = phase_alleles(alleles_all, family, kinship2ex)

% swap phase for optimal configuration
% rename alleles by founders
% only processes genotyped individuals

error = 0;
assignment = [];
global debug_mode;

%%

[nIND c] = size(family);
if( nIND <= 0 || c < 12 )
    error = 1;
    disp('error in family structures');
    return;
end

isGENO(1:nIND) = (family(1:nIND,7) == 1);
nGENO = nnz(isGENO);

[d1 d2 d3 d4] = size(kinship2ex);
if( d1 ~= nIND || d2 ~= d1 || d3 ~= 2 || d4 ~= 2 )
    error = 1;
    disp('error in kinship');
    return;
end


for i = 1:nIND
    if( isGENO(i) )
        if( any(alleles_all(i,1:2) == 0) )
            error = 1;
            disp('genotyped individual not assigned');
            return;
        end
    end
end

%%

% erase all founder tagging, minus labels
assignment = alleles_all;

phased(1:nIND) = false;

% must explicitly extend kinship to all individuals
% expand kinship to all individuals
exp_kin = kinship2ex;


% switch paternal and maternal alleles to be the most probable
% find all IBD alleles in the pedigree
% find the optimal phase, based on maximum kinship strength
for i = 1:nIND
    if( all(assignment(i,1:2) ~= 0, 2) )
        p = assignment(i,1);
        m = assignment(i,2);
        % in case p and m are singeltons, two summation will be equal
        reP = assignment(1:nIND,1) == p | assignment(1:nIND,2) == p;
        reM = assignment(1:nIND,1) == m | assignment(1:nIND,2) == m;        
        sum11 = sum(exp_kin(i,reP,1,1)) + sum(exp_kin(i,reP,1,2));
        sum12 = sum(exp_kin(i,reP,2,1)) + sum(exp_kin(i,reP,2,2));
        sum21 = sum(exp_kin(i,reM,1,1)) + sum(exp_kin(i,reM,1,2));
        sum22 = sum(exp_kin(i,reM,2,1)) + sum(exp_kin(i,reM,2,2));
        if( sum11 + sum22 < sum12 + sum21 )
            assignment(i,1) = m;
            assignment(i,2) = p;
        end
    end
end


% determine that one side of allele is of determined parental source
% this step can be applied iteratively
for i = 1:nIND
    if( all(assignment(i,1:2) ~= 0, 2) )
        p = assignment(i,1);
        m = assignment(i,2);
        % the status array includes i itself
        reP = assignment(1:nIND,1) == p | assignment(1:nIND,2) == p;
        reM = assignment(1:nIND,1) == m | assignment(1:nIND,2) == m;
        kinP = exp_kin(i, reP, 1, 1) + exp_kin(i, reP, 1, 2);
        kinM = exp_kin(i, reM, 2, 1) + exp_kin(i, reM, 2, 2);
        kinP_alt = exp_kin(i, reM, 1, 1) + exp_kin(i, reM, 1, 2);
        kinM_alt = exp_kin(i, reP, 2, 1) + exp_kin(i, reP, 2, 2);        
        if( (any(kinP <= 0) || any(kinM <= 0)) && (any(kinP_alt <= 0) || any(kinM_alt <= 0)) )
            phased(i) = false;
            if( debug_mode == 1 )
                disp('           non-pedigree sharing');
                disp(['           ', num2str(i), ': no paths via both parents']);
            end
            continue;
        end
        if( any(kinP <= 0) || any(kinM <= 0) )
            assignment(i,1) = m;
            assignment(i,2) = p;
            phased(i) = true;
        end
        if( any(kinP_alt <= 0) || any(kinM_alt <= 0) )
            assignment(i,1) = p;
            assignment(i,2) = m;
            phased(i) = true;
        end

    end
end


% ancestor or self closed genotyped shield
% include self
% if self not genotyped, trace to parents until founders
% assuming all relevant relatives are explicitly listed
[ngs error] = genotyped_envelope(family);
if( error ~= 0 )
    disp('error in generating ancestral shield');
    return;
end
        
% this check is necessary beyond kinship2 check
% still not use all infomation
% a shield is a closed ancestral cover, formed by closest ancestors
% all distant allele sharing must pass
% a detailed check involving testing each allele sharing
for i = 1:nIND
    if( all(assignment(i,1:2) ~= 0) )
        father = family(i,3);
        mother = family(i,4);
        p = assignment(i,1);
        m = assignment(i,2);
        
        if( father ~= 0 && mother ~= 0 )
            pshield = ngs(father, 1:nIND);
            mshield = ngs(mother, 1:nIND);
            if( any(pshield) && all(isGENO(pshield)) && any(mshield) && all(isGENO(mshield)) )
                condition1 = any(any(assignment(pshield,1:2) == p));
                condition2 = any(any(assignment(pshield,1:2) == m));
                condition3 = any(any(assignment(mshield,1:2) == p));
                condition4 = any(any(assignment(mshield,1:2) == m));
                if( (~condition1 || ~condition4) && (~condition2 || ~condition3) )
                    phased(i) = false;
                    if( debug_mode == 1 )
                        disp('           non-pedigree sharing');
                        disp(['           ', num2str(i), ': sharing not observed in ancestors']);
                    end
                    continue;
                end
            end
        end        
        
        % this is different from shield check in pairwise IBD correction
        % given each individual, the shield could be shrinked to only
        % relevant ancestors
        % such that more powerful
        
        % all sharing through father must come through pshield
        % check if m exist in pshield, if not, p must
        % otherwise, error
        if( father ~= 0 )
            pshield = ngs(father,1:nIND);
            if( any(pshield) && all(isGENO(pshield)) )
                if( ~any(any(assignment(pshield,1:2) == m)) )
                    phased(i) = true;
                    assignment(i,1) = p;
                    assignment(i,2) = m;
                end
                if( ~any(any(assignment(pshield,1:2) == p)) )
                    phased(i) = true;
                    assignment(i,1) = m;
                    assignment(i,2) = p;
                end
            end
        end
        
        % all sharing through mother must come through mshield
        % check if p exist in mshield, if not, m must
        % otherwise, error
        if( mother ~= 0 )
            mshield = ngs(mother,1:nIND);           
            if( any(mshield) && all(isGENO(mshield)) )
                if( ~any(any(assignment(mshield,1:2) == p)) )
                    phased(i) = true;
                    assignment(i,1) = p;
                    assignment(i,2) = m;
                end
                if( ~any(any(assignment(mshield,1:2) == m)) )
                    phased(i) = true;
                    assignment(i,1) = m;
                    assignment(i,2) = p;
                end
            end
        end

    end
end

assignment = remapFOUNDER(assignment);

end













