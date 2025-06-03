p1 = plot(p10, linewidth=3, color="#a9a9a9"); hold on;
p2 = plot(p50, linewidth=3, color="#71797e");
p3 = plot(p90, linewidth=3, color="#000"); 

stim_i = find(p_iden.stim(2, :));
x1 = stim_i(1);
x2 = stim_i(end);
p4 = patch([x1 x1 x2 x2], [20 18 18 20], "magenta", "FaceColor", "magenta", "FaceAlpha", .2, "LineStyle", "none");

rwd_i = find(p_iden.rwd);
x1 = rwd_i(1);
x2 = rwd_i(end);
p5 = patch([x1 x1 x2 x2], [20 18 18 20], "blue", "FaceColor", "blue", "FaceAlpha", .2, "LineStyle", "none"); hold off

h = [p1; p2; p3; p4; p5];

legend(h, "10%", "50%", "90%", "cue", "reward");