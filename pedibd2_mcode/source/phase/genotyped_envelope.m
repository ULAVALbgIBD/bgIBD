function [ngs error] = genotyped_envelope(family)

error = 0;
ngs = false;

% ancestor or self closed genotyped shield
% include self
% if self not genotyped, trace to parents until founders
% assuming all relevant relatives are explicitly listed

% ALL individuals in the envelop are genotyped

[nind c] = size(family);
if( nind <= 0 || c < 12 )
    error = 1;
    disp('error in family structures');
    return;
end

isgeno(1:nind) = (family(1:nind,7) == 1);

ngs(nind,nind) = false;
for i = 1:nind
    if( isgeno(i) )
        ngs(i,i) = true;
    else
        father = family(i,3);
        mother = family(i,4);
        if( father ~= 0 && mother ~= 0 )
            ngs(i,1:nind) = ngs(father,1:nind) | ngs(mother,1:nind);
        else
            ngs(i,i) = true;
        end
    end
end



end