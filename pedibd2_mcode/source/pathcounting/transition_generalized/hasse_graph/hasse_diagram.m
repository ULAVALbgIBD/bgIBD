
function new_edges = hasse_diagram(num)

% all subset combinations of size of num
% each set can choose to join a previous set, or emerge as a new one
% bell number

% the change number indicates the set representative
% the whole set is separated from the other one

if( num > 2 )
    old_edges = hasse_diagram(num - 1);
else
    new_edges(1,1:2) = [1,1];
    new_edges(1,3) = -2; % change at 2
    new_edges(1,4:5) = [1,2];
    return;
end

count = 0;
pro = old_edges(1:end,1:num-1);
change = old_edges(1:end,num);
suc = old_edges(1:end,num+1:num*2-1);

for i = 1:length(old_edges(:,1))
    for j = 1:num-1
        count = count + 1;
        new_pro(count,1:num-1) = pro(i,1:num-1);
        new_suc(count,1:num-1) = suc(i,1:num-1);
        new_pro(count,num) = new_pro(count,j);
        new_suc(count,num) = new_suc(count,j);
        new_change(count,1) = change(i,1);
        count = count + 1;
        new_pro(count,1:num) = new_suc(count-1,1:num);
        new_suc(count,1:num-1) = suc(i,1:num-1);
        new_suc(count,num) = num;
        new_change(count,1) = - num;
    end
    count = count + 1;
    new_pro(count,1:num-1) = pro(i,1:num-1);
    new_suc(count,1:num-1) = suc(i,1:num-1);    
    new_pro(count,num) = num;
    new_suc(count,num) = num;
    new_change(count,1) = change(i,1);
end


count = count + 1;
new_pro(count,1:num) = 1;
new_suc(count,1:num-1) = 1;
new_suc(count,num) = num;
new_change(count,1) = -num;


new_edges = [new_pro,new_change,new_suc];

