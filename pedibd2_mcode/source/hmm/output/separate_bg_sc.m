

function [output_dec_state] = separate_bg_sc(temp, sampled_markerlist)



output_dec_state(1:length(temp),1:2) = 1;
for i = 1:length(temp);
    if( temp(i) == 2 )
        output_dec_state(i,1) = output_dec_state(i,1) + 1;
    end
    if( temp(i) == 3 )
        output_dec_state(i,2) = output_dec_state(i,2) + 1;
    end
end





%%







