%%



deno = 1;
count = 0;
step = 200;
ld_smooth = [];
for i = 1:step:length(ld)-step
    count = count + 1;
    ld_smooth(count, 1) = sum(ld(i:i+step-1))/(step*deno);
    ld_smooth(count, 2) = markerlist(i+step/2, 2);
end

%%

subplot(2, 1, 1);
plot(markerlist(:,2), dec_state, 'r');
subplot(2, 1, 2);
plot(ld_smooth(:,2), ld_smooth(:,1));

%%

plotyy(markerlist(:,2), dec_state, ld_smooth(:,2), ld_smooth(:,1));

%%

% linkage disequilibrium

lld = get_ld_locus(data, member_list, 1, 2);
