function [pr, tpl, tpl_log, error] = generate_transition_all_loci(all_paths, parameters, de_ibd, de_tran)

error = 0;
tpl = [];
pr = [];

markerlist = parameters.sampled_markerlist; 

[nLOC, fields] = size(markerlist);

if( nLOC <= 0 || fields ~= 2 )
    disp('error in marker list');
    return;
end

mi(1:nLOC) = 0;
r(1:nLOC) = 0;

mi(1) = markerlist(1,2);
mi(2:nLOC) = markerlist(2:nLOC,2) - markerlist(1:nLOC-1,2);

% if there is genetic distance use genetic distance
if( isfield(parameters, 'genetic_map') )
    r = 0.01 * parameters.genetic_map;
else
    r(1:nLOC) = mi(1:nLOC) * 0.01 / 1000000;
end


r(1) = mean(r(2:end));

[pr, tpl_loc] = generate_transition_plus(all_paths, r(1), de_ibd, de_tran);

[nrows ncols] = size(tpl_loc);
if( nrows ~= ncols || nrows <= 0 )
    error = 1;
    disp('error in generating transition probabilities');
    return;
end
n_states = nrows;


flag(1:200) = false;
tpl_static(1:200, 1:n_states, 1:n_states) = 0;


tick(1:nLOC) = 0;
tick(mi <= 0) = 1;
tick(mi > 0) = round(10 * log10(mi(mi>0))) + 1;


for i = 1:nLOC
    % probability of change from marker i-1 to i
    if( tick(i) > 200 )
        error = 1;
        disp('error: unexpected large marker interval');
        return;
    end
    if( ~flag(tick(i)) )
        frac = (10^((tick(i)-1)/10))/(10^8);
        [pr, tpl_loc] = generate_transition_plus(all_paths, frac, de_ibd, de_tran);
        [nrows ncols] = size(tpl_loc);
        if( nrows ~= n_states || ncols ~= n_states )
            error = 1;
            disp('error in generating transition probabilites');
            return;
        end
        flag(tick(i)) = true;
        tpl_static(tick(i), 1:n_states, 1:n_states) = tpl_loc(1:n_states, 1:n_states);
    end
    
end

temp = tpl_static;
valid = tpl_static>0;
temp(~valid) = -30;
temp(valid) = reallog(tpl_static(valid));

tpl = tpl_static(tick(1:nLOC), 1:n_states, 1:n_states);
tpl_log = temp(tick(1:nLOC), 1:n_states, 1:n_states);


end







