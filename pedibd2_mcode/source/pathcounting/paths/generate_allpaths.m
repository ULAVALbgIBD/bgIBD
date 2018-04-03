function [singleLINEAL, error] = generate_allpaths(family)

% freeALLELE
global freeALLELE;

% all inheritance path from i to j

[rows, cols] = size( family );
if( rows <= 0 || cols < 12 )
    error = 1;
    disp('error in family structure');
    return;
end
nIND = rows;

P = false(nIND, nIND);
for i = 1:nIND
    if( family(i, 3) > 0 )
        P(i, family(i,3)) = true;
    end
    if( family(i, 4) > 0 )
        P(i, family(i,4)) = true;
    end
end

% immediate predecessor if diffused from i to j
freeALLELE = zeros(nIND, nIND, 2);
for i = 1:nIND  
    % first argument does not change
    % for marking purposes
    freeALLELE(i,i,1:2) = -1;    % reference individual is multiple-sourced
    % identify all free individuals
    traverse(i, i, P, true);
end

singleLINEAL = freeALLELE;


clear global freeALLELE;

error = 0;


end



