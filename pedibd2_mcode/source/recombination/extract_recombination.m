

function [recombination error] = extract_recombination(sibship, family, intervals, parameters)

error = 0;
recombination = [];
global debug_mode;


status = sibship.sibship_status;

[nIND, fields] = size(family);
if( nIND <= 0 || fields < 12 )
    error = 1;
    disp('error in family structures');
    return;
end

[nSEG, c] = size( intervals );
if( nSEG <= 0 || c < 2 )
    error = 1;
    disp('error in global IBD');
    return;
end

[d1, d2, d3] = size(status);
if( d1 <= 0 || d1 ~= nIND || d2 ~= nIND || d3 <= 0 || d3 ~= nSEG )
    error = 1;
    disp('error in extracting recombination positions');
    return;
end


if( isempty(parameters) )
    error = 1;
    return;
end
[nLOC, cols] = size( parameters.sampled_markerlist );
if( cols ~= 2 || nLOC < nSEG )
    error = 1;
    return;
end
basepair(1:nLOC) = parameters.sampled_markerlist(1:nLOC, 2);

genotyped = family(1:nIND,7);



% smooth recombination breakpoints, close points are marked as invalid

% average sharing length of 100 meiosis apart
cutoff = 10^6;
nearby = 10 * cutoff;
% status 0,1 two states
% status -1 is an undefined state
% smooth the line
for i = 1:nIND
    for j = 1:nIND
        status_line(1:nSEG) = status(i,j,1:nSEG);
        if( ~any(status_line >= 0) )
            continue;
        end
        pre_seg = 1;
        for k = 1:nSEG
            if( status_line(k) >= 0 )
                first_seg = k;
                break;
            end
        end
        if( first_seg == nSEG )
            continue;
        end
        pre_seg = first_seg;
        for k = first_seg+1:nSEG
            cur_seg = k;
            pre_status = status_line(pre_seg);
            cur_status = status_line(cur_seg);            
            if( cur_status >= 0 && cur_status ~= pre_status )
                marker1 = intervals(pre_seg,1);
                marker2 = intervals(cur_seg-1,2);
                if( basepair(marker2) - basepair(marker1) < cutoff )
                    status_line(pre_seg:cur_seg-1) = -1;
                else
                    % let stay the same
                end
                pre_seg = cur_seg;
            end
        end
        valid_line(1:nSEG) = status(i,j,1:nSEG);
        status(i,j,valid_line>=0) = status_line(valid_line>=0);
    end
end


sibship_count = 0;
for i = 1:nIND
    for j = 1:nIND
        status_line(1:nSEG) = status(i,j,1:nSEG);
        if( ~any(status_line >= 0) )
            continue;
        end
        sibship_count = sibship_count + 1;
    end
end


if( debug_mode == 1 && sibship_count <= 10 )
    c = 0;
    for i = 1:nIND
        for j = 1:nIND
            status_line(1:nSEG) = status(i,j,1:nSEG);
            if( ~any(status_line >= 0) )
                continue;
            end           
            c = c + 1;
            subplot(sibship_count,1,c);
            temp(1:nSEG) = 0.5;
            temp(status_line(1:nSEG) >= 0) = status_line(status_line(1:nSEG) >= 0);
            for k = 1:nSEG
                line(intervals(k,1:2), [temp(k), temp(k)]);
            end
            text(intervals(fix((1+nSEG)/2)), 0.5, num2str(family(i,8)));
            ylim([-0.1,1.1]);      
        end
    end
end

list = [];
nrec = 0;
for i = 1:nIND
    for j = 1:nIND
        status_line(1:nSEG) = status(i,j,1:nSEG);
        if( ~any(status_line >= 0) )
            continue;
        end
        if( genotyped(i) == 0 && genotyped(j) == 0 )
            if( i >= j )
                continue;
            end
        end
        pre_seg = 1;
        for k = 2:nSEG
            cur_seg = k;
            pre_status = status_line(pre_seg);
            cur_status = status_line(cur_seg);
            if( cur_status >= 0 )
                if( pre_status >= 0 )
                    if( cur_status ~= pre_status )
                        marker1 = intervals(pre_seg,2);
                        marker2 = intervals(cur_seg,1);
                        nrec = nrec + 1;
                        list(nrec, 1) = marker1;
                        list(nrec, 2) = marker2;
                        if( genotyped(i) == 1 || genotyped(j) == 1 )
                            list(nrec, 3) = i;
                            list(nrec, 4) = 0;
                        else
                            list(nrec, 3) = i;
                            list(nrec, 4) = j;
                        end
                        if( debug_mode == 1 )
                            disp(['          ', num2str([marker1,marker2], '(%6d,%6d)'), ' recombination: gamete of ', num2str(family(i,8)), ' or ', num2str(family(j,8))]);
                        end                        
                    end
                end
                pre_seg = cur_seg;
            end
        end        
    end
end


valid = true(nSEG, 1);
for i = 1:nIND
    for j = 1:nIND
        status_line(1:nSEG) = status(i,j,1:nSEG);
        if( ~any(status_line >= 0) )
            continue;
            % whether it is a valid parent pair
        end
        for k = 1:nSEG
            if( status_line(k) < 0 )
                valid(k) = 0;
            end
        end        
    end
end

recombination.list = list;
recombination.valid = valid;



end




































