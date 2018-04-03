

%%%%%%%%transition probability%%%%%%%%%%%%%%%%%%%%%%%%%%
function [input_pr input_tpl tpl_log status status_bypass error] = transition(input_allpaths, parameters)


input_pr = [];
input_tpl = [];
tpl_log = [];
status = [];
status_bypass = [];
error = 0;

sampled_markerlist = parameters.sampled_markerlist;

[nloc, fields] = size(sampled_markerlist);
if( ndims(sampled_markerlist) ~= 2 || nloc <= 0 || fields ~= 2 )
    disp('error in marker list');
    error = 1;
    return;
end

[npairs, fields] = size( input_allpaths );
if( ndims(input_allpaths) ~= 2 || npairs ~= 6 || fields < 1 )
    disp('error in inheritance paths');
    error = 1;
    return;
end


pt11 = (input_allpaths(2,:)); %pp
pt12 = (input_allpaths(3,:)); %pm
pt21 = (input_allpaths(4,:)); %mp
pt22 = (input_allpaths(5,:)); %mm

status(1:2,1:2) = 0;

% prio probability

status(1,1) = path_tran(pt11, 0);   %pp
status(1,2) = path_tran(pt12, 0);   %pm
status(2,1) = path_tran(pt21, 0);   %mp
status(2,2) = path_tran(pt22, 0);   %mm

sum_all = sum(sum(status));

input_pr = 0;
input_tpl = 0;

if( sum_all == 0 )
    status_bypass = 1;
    return;
end

% non-inbreeding parent-child bypass
if( status(1,1) + status(1,2) == 1 )
    if( status(2,1) + status(2,2) == 0 )
        status_bypass = 2;
        return;
    end
end

if( status(2,1) + status(2,2) == 1 )
    if( status(1,1) + status(1,2) == 0 )
        status_bypass = 2;
        return;
    end
end

if( status(1,1) + status(2,1) == 1 )
    if( status(1,2) + status(2,2) == 0 )
        status_bypass = 2;
        return;
    end
end

if( status(1,2) + status(2,2) == 1 )
    if( status(1,1) + status(2,1) == 0 )
        status_bypass = 2;
        return;
    end
end
    
status_bypass = 0;

min_prio = 1;
for i = 1:2
    for j = 1:2
        if( status(i,j) > 0 )
            min_prio = min(min_prio, status(i,j));
        end
    end
end
% input_tpl{i}, is the transition from locus i-1 to locus i

%%%%%%%%%%%%background ibd%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% maybe above foreground sharing level for distant relatives
bg_ibd = 0.01;

% prio probability determines number of lines in the background ibd

bg_tran = 400;    %generations apart

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% transition probability is not correct for parent-child relationship
[input_pr, input_tpl, tpl_log, error] = generate_transition_all_loci(input_allpaths, parameters, bg_ibd, bg_tran);





end


