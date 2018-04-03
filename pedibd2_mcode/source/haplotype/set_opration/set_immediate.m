function set_immediate( index, anc, off )

global ancestral;
global complexity;

complexity = complexity + 1;

if( index < 0 )
    ancestral(-index,1) = anc;
    ancestral(-index,3) = off;
else
    ancestral(index,2) = anc;
    ancestral(index,4) = off;
end

end

