function [prob_f tran_matrix] = transition_foreground(inheritance, markerlist)


% change accumulated in one column or one row
% ip is inheritance path

prob_f = inheritance.prob_f;
tran_list = inheritance.tran_list;

r(1:length(markerlist(:,1))) = 0;

for i = 1:length(markerlist(:,1))    
    if( i == 1 )
        r(i) = 0.01*(markerlist(i,2))/1000000;
    else
        r(i) = 0.01*(markerlist(i,2)-markerlist(i-1,2))/1000000;
    end
end

r(1) = mean(r(2:end));

for i = 1:length(r)
    tran_matrix{i} = transition_foreground_aux(prob_f, tran_list, r(i));
end

end










