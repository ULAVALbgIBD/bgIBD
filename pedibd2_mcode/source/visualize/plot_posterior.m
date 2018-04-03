
function plot_posterior(geno, oblist, posterior, posteriorIBD1, viterbi, markerlist, id1, id2, limit, option)


count = 0;
ks = limit(1);
js = limit(2);

count = count + 1.5;
hold on;

global disp_mode;
if( disp_mode == 0 )
    coor = markerlist(:,1);
else
    coor = markerlist(:,2);
end
map = markerlist(:,2);

[temp1, temp2] = max(posterior(:,:), [], 2);

if( option == 0 )
    lim_co = coor(ks:js);
    lim_ob = oblist(ks:js);
    lim_vi = viterbi(ks:js,1);
    temp3 = lim_ob <= 1;
    plot(lim_co(temp3), random_ibs(lim_ob(temp3))-1, '.', 'MarkerSize', 3);
    temp4 = lim_vi == 2;
    plot(lim_co(temp4), lim_vi(temp4) - 2, 'gx');
    temp5 = lim_vi == 3;
    plot(lim_co(temp5), lim_vi(temp5) - 3, 'r+');
%     ylim([-0.6,1.1]);
    text(-10^6,0,[num2str(id1, '%02d'),' vs ',num2str(id2, '%02d')],'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    
else
    text((coor(ks)+coor(js))/2,count*3+2.5,[num2str(id1),'-',num2str(id2)]);  
    
    % posterior probability is corrected on the fly
    [smoothed_viterbi error] = smooth_viterbi1pair(oblist, viterbi, posterior, posteriorIBD1, map);
    
    plot(coor(ks:js), count*3 + 7 + (smoothed_viterbi(ks:js,1:2)), 'linewidth', 3);
    plot(coor(ks:js), count*3 + 4 + (viterbi(ks:js,1:2)), 'linewidth', 3);
    plot(coor(ks:js), count*3 + 3 * (posterior(ks:js,:)));
    plot(coor(ks:js), count*3 + random_ibs(oblist(ks:js), []), '.');
    if( js - ks < 500 )
        for i = ks:js
            for j = 1:4
                text(coor(i), j-0.1, num2str(geno(j,i,1)));
                text(coor(i), j+0.1, num2str(geno(j,i,2)));
            end
        end
    end
    ylim([0,15]);
end


end


%     ibs0 = 1;
%     ibs1 = 2;
%     ibs2_homo = 3;
%     ibs2_hete = 4;

