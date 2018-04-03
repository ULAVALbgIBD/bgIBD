function [ngs error] = genotyped_ancestor(family)

error = 0;
ngs = false;

% nearest genotyped shield
% trace to nearest genotyped ancestral source, if none, trace to founders
% the shield is flawless
% assuming ancestor always preceding offspring

% not ALL individuals in the envelop are genotyped

[nIND c] = size(family);
if( nIND <= 0 || c < 12 )
    error = 1;
    disp('error in family structures');
    return;
end


isGENO(1:nIND) = (family(1:nIND, 7) == 1);
ngs(1:nIND,1:nIND) = false; 
for i = 1:nIND
    father = family(i,3);
    mother = family(i,4);
    if( father ~= 0 )
        if( isGENO(father) || ~any(ngs(father, 1:nIND)))
            ngs(i, father) = true;
        else
            ngs(i, 1:nIND) = ngs(i, 1:nIND) | ngs(father, 1:nIND);
        end
    end
    if( mother ~= 0 )
        if( isGENO(mother) || ~any(ngs(mother, 1:nIND)) )
            ngs(i, mother) = true;
        else
            ngs(i, 1:nIND) = ngs(i, 1:nIND) | ngs(mother, 1:nIND);
        end
    end
end



end