

function [output output_ibd1 error]= posterior_path(paths, ob)

error = 0;
output = [];
output_ibd1 = [];

global debug_mode;

oblist = ob;    %no need to map, ready states
I_prob = paths.pr;
T_prob = paths.tpl;
O_prob = paths.epl;
[d1, d2, d3] = size(T_prob);
if( d1 <= 0 || d2 <= 0 || d2 ~= d3 )
    error = 1;
    disp('error in transition states');
    return;
end
nloc = d1;
nstates = d2;
[d1, d2, d3] = size(O_prob);
if( d1 ~= nloc || d2 ~= nstates || d3 <= 0 )
    error = 1;
    disp('error in emission states');
    return;
end
len = length(I_prob);
if( len <= 0 || len ~= nstates )
    error = 1;
    disp('error in prior probability');
    return;
end
len = length(oblist);
if( len <= 0 || len ~= nloc )
    error = 1;
    disp('error in observation list');
    return;
end

hidden_states = [1:nstates];

[alpha error] = get_forward(I_prob, T_prob, O_prob, oblist, hidden_states);
if( error ~= 0 )
    disp('error in generating forward probability');
    return;
end
[beta error] = get_backward(I_prob, T_prob, O_prob, oblist, hidden_states);
if( error ~= 0 )
    disp('error in generating backward probability');
    return;
end



% [1p2p, 1p2m, 1m2p, 1m2m] included states
states = [
    1, 1, 1, 1; %1  %ibd0
    2, 1, 1, 1; %2  %ibd1_ll
    3, 1, 1, 1; %3
    1, 2, 1, 1; %4  %ibd1_lr
    1, 3, 1, 1; %5
    1, 1, 2, 1; %6  %ibd1_rl
    1, 1, 3, 1; %7
    1, 1, 1, 2; %8  %ibd1_rr
    1, 1, 1, 3; %9
    2, 1, 1, 2; %10 %ibd2
    3, 1, 1, 2; %11 %ibd1_rr
    2, 1, 1, 3; %12 %ibd1_ll
    3, 1, 1, 3; %13
    1, 2, 2, 1; %14
    1, 3, 2, 1; %15 %ibd1_rl
    1, 2, 3, 1; %16 %ibd1_lr
    1, 3, 3, 1; %17
    ];

ibd0 = [1,3,5,7,9,13,17];
ibd1 = [2,4,6,8,11,12,15,16];
ibd2 = [10,14];

ibd1_ll = [2,12];
ibd1_lr = [4,16];
ibd1_rl = [6,15];
ibd1_rr = [8,11];


output(1:nloc,1:3) = 0;
output_ibd1(1:nloc,1:2,1:2) = 0;

pos(1:nstates,1:nloc) = times(alpha(1:nstates,1:nloc) , beta(1:nstates,1:nloc));

pos(1:nstates,1:nloc) = pos(1:nstates,1:nloc)./repmat(sum(pos(1:nstates,1:nloc),1),nstates,1);

output(1:nloc,1) = sum(pos(ibd0,1:nloc), 1);
output(1:nloc,2) = sum(pos(ibd1,1:nloc), 1);
output(1:nloc,3) = sum(pos(ibd2,1:nloc), 1);

output_ibd1(1:nloc,1,1) = sum(pos(ibd1_ll,1:nloc), 1);
output_ibd1(1:nloc,1,2) = sum(pos(ibd1_lr,1:nloc), 1);
output_ibd1(1:nloc,2,1) = sum(pos(ibd1_rl,1:nloc), 1);
output_ibd1(1:nloc,2,2) = sum(pos(ibd1_rr,1:nloc), 1);



end





















