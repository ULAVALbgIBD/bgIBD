


function moveCLOSEtoI(i, j)

global mostCLOSEtoI;
% most CLOSE to I
% mostCLOSEtoI(i,j,k)
% j selects at most one allele of i from k, 

nIND = size(mostCLOSEtoI,1);

for k = 1:nIND
    if( mostCLOSEtoI(i,j,k) )
        % update each of the closer relative k between
        % i and j
        % if k could be moved closer to i
        % recode it
        if( k ~= i )
            moveCLOSEtoI(i, k);
            closerR = find(mostCLOSEtoI(i, k, :));
            if( length(closerR) == 1 && closerR ~= i )
                if( ~mostCLOSEtoI(i,j,closerR) )
                    % if yes, already multiple selection
                    mostCLOSEtoI(i,j,k) = false;
                    mostCLOSEtoI(i,j,closerR) = true;
                end
            end
        end
    end
end

end


