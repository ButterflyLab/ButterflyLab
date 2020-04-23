function LT_scaling(startT, st, md, cut)

% N (log N)^i
i = [3, 3, 3, 3];

D = load(['./results/' startT '.mat']);

X = log2(D.N); X = X(st:end-cut);
Yff = log2(BF_trimdata(D.FacT_f)); Yff = Yff(st:end-cut);
Yfi = log2(BF_trimdata(D.FacT_i)); Yfi = Yfi(st:end-cut);
Yaf = log2(BF_trimdata(D.AppT_f)); Yaf = Yaf(st:end-cut);
Yai = log2(BF_trimdata(D.AppT_i)); Yai = Yai(st:end-cut);

c1 = {'#EDB120', 'g', 'r', 'b'};
mk = {'^', 'v', '<', '>'};
c2 = {'', 'm', 'c', 'k'};

pic = figure('visible', 'off');
hold on;
h(1) = plot(X, Yff, '-', 'Color', c1{1}, 'Marker', mk{1}, 'LineWidth', 3);
h(2) = plot(X, Yfi, '-', 'Color', c1{2}, 'Marker', mk{2}, 'LineWidth', 3);
h(3) = plot(X, Yaf, '-', 'Color', c1{3}, 'Marker', mk{3}, 'LineWidth', 3);
h(4) = plot(X, Yai, '-', 'Color', c1{4}, 'Marker', mk{4}, 'LineWidth', 3);
vec = X+(i(1)-1)*log2(X);
h(5) = plot(X, Yff(md)+vec-vec(md), '--', 'Color', c2{i(1)}, 'LineWidth',2);
vec = X+i(1)*log2(X);
h(6) = plot(X, Yff(md)+vec-vec(md), '--', 'Color', c2{i(1)+1}, 'LineWidth',2);
vec = X+(i(2)-1)*log2(X);
plot(X, Yfi(md)+vec-vec(md), '--', 'Color', c2{i(2)}, 'LineWidth',2);
vec = X+i(2)*log2(X);
plot(X, Yfi(md)+vec-vec(md), '--', 'Color', c2{i(2)+1}, 'LineWidth',2);
vec = X+(i(3)-1)*log2(X);
plot(X, Yaf(md)+vec-vec(md), '--', 'Color', c2{i(3)}, 'LineWidth',2);
vec = X+i(3)*log2(X);
plot(X, Yaf(md)+vec-vec(md), '--', 'Color', c2{i(3)+1}, 'LineWidth',2);
vec = X+(i(4)-1)*log2(X);
plot(X, Yai(md)+vec-vec(md), '--', 'Color', c2{i(4)}, 'LineWidth',2);
vec = X+i(4)*log2(X);
plot(X, Yai(md)+vec-vec(md), '--', 'Color', c2{i(4)+1}, 'LineWidth',2);

lb = {'N', 'N log_2(N)', 'N log_2^2(N)', 'N log_2^3(N)'};
legend(h,{'Fwd Fac','Inv Fac','Fwd App','Inv App',...
    lb{i(1)},lb{i(1)+1}},'Location','eastoutside');

axis tight;
axis square;
xlabel('log_2(N)'); ylabel('log_2(Time)');
set(gca, 'FontSize', 16);
b = get(gca);
set(b.XLabel, 'FontSize', 16); set(b.YLabel, 'FontSize', 16);
set(b.ZLabel, 'FontSize', 16); set(b.Title, 'FontSize', 16);
hold off;

saveas(gcf,['./results/' startT '.pdf'])

close(pic);

end
