
function t_ibs = random_ibs(input, scale)

    sz = size(input);
    if( isempty(scale) )
        scale = 0.2;
    end
    t_ibs = double(input) + 1 * scale * (rand(sz(1), sz(2))-0.5);
    
end