
scrsz = get(0,'ScreenSize');
figure('OuterPosition',[1 scrsz(4)/3 scrsz(3)*3/4 scrsz(4)/3]);
clear scrsz;

number = 1;

sbps = 1;
ebps = length(sampled_markerlist);

% sbps = 12900;
% ebps = 13050;
% 
% sbps = 28000;
% ebps = 30000;

% sbps = 8500;
% ebps = 8700;

% sbps = 13000;
% ebps = 13400;

% sbps = 6000;
% ebps = 8000;



plot_number = 1;

x = sampled_markerlist(sbps:ebps, 2);
% y1 = input_oblist_genotype{s_id}(sbps:ebps);

y1 = input_oblist_ibs{s_id}(sbps:ebps);
y2_1 = output_dec_state{s_id}(sbps:ebps,1) - 1 + 0.2;
y2_2 = output_dec_state{s_id}(sbps:ebps,2)- 1 - 0.2;

% y1 = input_oblist_genotype_sc{s_id_sc}(sbps:ebps);
% y2_1 = output_dec_state_sc{s_id_sc}(sbps:ebps, 1) - 1 + 0.2;
% y2_2 = output_dec_state_sc{s_id_sc}(sbps:ebps, 2) - 1 + 0.2;

% y2_2 = random_ref{iteration, s_id}(sbps:ebps) - 0.2;


y3 = (output_del(:,sbps:ebps)');
lines = [1,2,3];
for i = 1:length(y3)
    temp = mean(y3(i,lines));
    y3(i,lines) = y3(i,lines) - temp;
end

line1(1:ebps-sbps+1) = 1;

temp(1:length(markerlist),1:length(lines)) = 0;
for i = 1:length(markerlist);
    temp(i,1:length(lines)) = input_epl{i}(lines,input_oblist_ibs{s_id}(i)+1);
end

interval = 10;

% temp(:,1) = smooth(temp(:,1), interval);
% temp(:,2) = smooth(temp(:,2), interval);
% temp(:,2) = smooth(temp(:,3), interval);

y4 = temp(sbps:ebps, :);
for i = 1:length(lines)
    y4(:,i) = smooth(log(y4(:,i)), interval);
end

y5 = smooth(sampled_af(sbps:ebps,2),interval); %major allele frequency



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
    ylim([-0.5 3.5]);  
    
    if( plot_number >= 2 )
        subplot(plot_number*number, 1, plot_number*(i - 1) + 2);
        plot(x(a:b), y3(a:b, lines));
    end
    
    if( plot_number >= 3 )
        subplot(plot_number*number, 1, plot_number*(i - 1) + 3);
        plot(x(a:b), y4(a:b, 1:length(lines)));
        hold on;
        plot(x(a:b), line1(a:b));        
    end
    
    if( plot_number >= 4 )
        subplot(plot_number*number, 1, plot_number*(i - 1) + 4);
        plot(x(a:b), y5(a:b));

    end       
end

title('IBD sharing between a pair of siblings');
title(['IBD sharing between UM-', num2str(pedigree_all_missing(input_pairlist(s_id,1),8)), ' and UM-', num2str(pedigree_all_missing(input_pairlist(s_id,2),8))]);
title(['IBD sharing between UM-', num2str(input_pairlist(s_id,1)), ' and UM-', num2str((input_pairlist(s_id,2)))]);
% title(['IBD sharing between UM-', num2str(pedigree_all_missing(input_pairlist(s_id,4),2)), ' and UM-', num2str(pedigree_all_missing(input_pairlist(s_id,5),2))]);


ylabel('IBS, IBD, bg-IBD');
ylabel('IBS, inferred IBD, actual IBD');

xlabel('Chromosomal location (bps)');
set(gca, 'YTick', [0,1,2]);
xlim([0,max(sampled_markerlist(:,2))]);

hold off;
clear x y1 y2 y2_1 y2_2 y3 y4 y5 line1 a b i sbps ebps t_ibs temp lines number interval;

%%

% temp = dec_state(10980:11020);
% plot(temp);


%%





