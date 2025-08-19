function fig = make_linegraph_multiscale (contrast_results, p_iden) 
    fig = figure;
    p1 = plot(contrast_results(1, :), linewidth=3, color="#cc0066"); hold on;
    p2 = plot(contrast_results(2, :), linewidth=3, color="#9933ff");
    p3 = plot(contrast_results(3, :), linewidth=3, color="#0066cc"); 
    p4 = plot(contrast_results(4, :), linewidth=3, color="#00ccff"); 
    
    %stim_i = find(p_iden.stim(2, :));
    %x1 = stim_i(1);
    %x2 = stim_i(end);
    %p5 = patch([x1 x1 x2 x2], [20 18 18 20], "magenta", "FaceColor", "magenta", "FaceAlpha", .2, "LineStyle", "none");
    
    %rwd_i = find(p_iden.rwd);
    %x1 = rwd_i(1);
    %x2 = rwd_i(end);
    %p6 = patch([x1 x1 x2 x2], [20 18 18 20], "blue", "FaceColor", "blue", "FaceAlpha", .2, "LineStyle", "none"); hold off
    
    %h = [p1; p2; p3; p4; p5; p6];
    
    %legend(h, "10%", "50%", "90%", "cue", "reward");
end