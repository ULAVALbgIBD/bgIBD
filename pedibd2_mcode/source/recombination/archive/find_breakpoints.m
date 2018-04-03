

function [breakpoints] = find_breakpoints(output_dec_state, sampled_markerlist)

%detect breakpoints

centromere = [12941,12942]; % fengyu's data, chromosome 6
centromere = [0,0]; % bypassing centromere filter
count = 0;
breakpoints = [];

for i = 1:length(output_dec_state(:,1))-1
    if( output_dec_state(i,1) ~= output_dec_state(i+1,1) )
        count = count + 1;
        breakpoints(count) = i;
        if( count >= 2 )
            % centromere filter
            if( breakpoints(count-1) <= centromere(1) && breakpoints(count)+1 >= centromere(2) )
                if( breakpoints(count) - breakpoints(count-1) < 400 )
                    count = count - 2;
                end
            end
        end
    end
end






%%







