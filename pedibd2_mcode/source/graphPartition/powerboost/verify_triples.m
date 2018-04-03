function [error] = verify_triples(triples_view, nseg, ngeno)

error = 0;
global debug_mode;

if( ngeno < 1 )
    return;
end

if( isempty(triples_view) )
    error = 1;
    return;
end

num_triples = triples_view.num_triples;
triple_viewpair = triples_view.triple_viewpair;
ambig_triples = triples_view.ambig_triples;
deter_triples = triples_view.deter_triples;

[d1, d2, d3] = size(deter_triples);
if( d1 ~= ngeno || d1 ~= d2 || d2 ~= d3 )
    error = 1;
    disp('error in triplets');
    return;
end
[d1, d2, d3] = size(ambig_triples);
if( d1 ~= ngeno || d1 ~= d2 || d2 ~= d3 )
    error = 1;
    disp('error in triplets');
    return;
end
[d0, d1, d2, d3] = size(triple_viewpair);
if( d0 ~= nseg || d1 ~= ngeno || d1 ~= d2 || d2 ~= d3 )
    error = 1;
    disp('error in triplets');
    return;
end
if( max(max(max(ambig_triples))) > num_triples )
    error = 1;
    disp('error in triplets');
    return;
end

if( ngeno < 3 )
    return;
end

for i = 1:ngeno
    for j = 1:ngeno
        for k = 1:ngeno
            if( ambig_triples(i,j,k) > 0 )
                if( ambig_triples(j,i,k) <= 0 )
                    error = 1;
                    return;
                end
                if( ambig_triples(j,k,i) <= 0 )
                    error = 1;
                    return;
                end
                if( ambig_triples(k,i,j) <= 0 )
                    error = 1;
                    return;
                end
                if( ambig_triples(k,j,i) <= 0 )
                    error = 1;
                    return;
                end
            end
            if( ambig_triples(i,j,k) > num_triples )
                error = 1;
                return;
            end
        end
    end
end


for i = 1:ngeno
    for j = 1:ngeno
        for k = 1:ngeno
            if( ambig_triples(i,j,k) > 0 )
                status1 = triple_viewpair(1:nseg,i,j,k);
                if( any(triple_viewpair(1:nseg,i,k,j) ~= status1) )
                    error = 1;
                    return;
                end
                status2 = triple_viewpair(1:nseg,j,i,k);
                if( any(triple_viewpair(1:nseg,j,k,i) ~= status2) )
                    error = 1;
                    return;
                end
                status3 = triple_viewpair(1:nseg,k,i,j);
                if( any(triple_viewpair(1:nseg,k,j,i) ~= status3) )
                    error = 1;
                    return;
                end
                if( any(status1 ~= status2) || any(status2 ~= status3) || any(status1 ~= status3) )
                    if( debug_mode == 1 )
                        error = 1;
                        disp('non-symmetric triplets');
                        return;
                    end
                end
            end
        end
    end
end


end













