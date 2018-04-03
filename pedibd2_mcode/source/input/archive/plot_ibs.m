
function plot_ibs(input)

    t_ibs = input + 1 * 0.1 * (rand(1, length(input))-0.5);
    t_ibs = abs(t_ibs);
    plot((t_ibs), '.', 'MarkerSize', 1);
    %ylim([-2.5 2.5]);
    ylim([-0.5 2.5]);
    
end