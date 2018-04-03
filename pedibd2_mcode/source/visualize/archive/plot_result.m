
function plot_result(oblist, dec, sampled_markerlist)

scrsz = get(0,'ScreenSize');
figure('OuterPosition',[1 scrsz(4)/3 scrsz(3)*3/4 scrsz(4)/3]);
clear scrsz;

number = 1;

sbps = 1;
ebps = length(sampled_markerlist);


plot_number = 1;

x = sampled_markerlist(sbps:ebps, 2);
% y1 = input_oblist_genotype{s_id}(sbps:ebps);



y1 = oblist(sbps:ebps);
y2_1 = dec(sbps:ebps,1) - 1 + 0.2;
y2_2 = dec(sbps:ebps,2)- 1 - 0.2;

% y1 = input_oblist_genotype_sc{s_id_sc}(sbps:ebps);
% y2_1 = output_dec_state_sc{s_id_sc}(sbps:ebps, 1) - 1 + 0.2;
% y2_2 = output_dec_state_sc{s_id_sc}(sbps:ebps, 2) - 1 + 0.2;

% y2_2 = random_ref{iteration, s_id}(sbps:ebps) - 0.2;


for i = 1:1:number
    subplot(plot_number*number, 1, plot_number*(i - 1) + 1);
    hold on;
    a = (i-1)*length(y1)/number + 1;
    b = i*length(y1)/number;
    a = round(a);
    b = round(b);
    t_ibs = random_ibs(y1(a:b));
    plot(x(a:b), t_ibs, 'g.', 'MarkerSize', 2);
    plot(x(a:b), y2_1(a:b), ':b+', 'MarkerSize', 2);
    plot(x(a:b), y2_2(a:b), 'rx', 'MarkerSize', 0.5);
%     plot(x(a:b), y2_1(a:b), 'b');
%     plot(x(a:b), y2_2(a:b), 'r');
    ylim([-0.5 max([max(y1), max(y2_1), max(y2_1)]) + 0.5]);            
end

% title('IBD sharing between a pair of siblings');
% title(['IBD sharing between UM-', num2str(pedigree_all_missing(input_pairlist(s_id,1),8)), ' and UM-', num2str(pedigree_all_missing(input_pairlist(s_id,2),8))]);
% title(['IBD sharing between UM-', num2str(input_pairlist(s_id,1)), ' and UM-', num2str((input_pairlist(s_id,2)))]);
% title(['IBD sharing between UM-', num2str(pedigree_all_missing(input_pairlist(s_id,4),2)), ' and UM-', num2str(pedigree_all_missing(input_pairlist(s_id,5),2))]);


ylabel('IBS, IBD, bg-IBD');
ylabel('IBS, inferred IBD, actual IBD');

xlabel('Chromosomal location (bps)');
set(gca, 'YTick', [0,1,2]);
xlim([0,max(sampled_markerlist(:,2))]);



end





