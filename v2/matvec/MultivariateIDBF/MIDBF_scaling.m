function MIDBF_scaling(F, st, md, cut)

% N (log N)^i
if startsWith(F,'2D')
    i = [1, 1, 1, 0];
elseif startsWith(F,'3D')
    i = [1, 1, 2, 1];
elseif startsWith(F,'Hel')
    i = [0, 1, 1, 0];
end

D = load(['./results/' F '.mat']);
X = log2(D.Nsq); X = X(st:end-cut);
Yo = log2(BF_trimdata(D.T_order)); Yo = Yo(st:end-cut);
Yr = log2(BF_trimdata(D.T_rec)); Yr = Yr(st:end-cut);
Yf = log2(BF_trimdata(D.T_fac)); Yf = Yf(st:end-cut);
Ya = log2(BF_trimdata(D.T_app)); Ya = Ya(st:end-cut);

c1 = {'#EDB120', 'k', 'b', 'r'};
c2 = {'g', 'm', 'c'};

pic = figure('visible', 'off');
hold on;
h(1) = plot(X, Yo, '-o', 'MarkerSize', 10, 'Color', c1{1}, 'LineWidth', 3);
h(2) = plot(X, Yr, '-s', 'MarkerSize', 10, 'Color', c1{2}, 'LineWidth', 3);
h(3) = plot(X, Yf, '-^', 'MarkerSize', 10, 'Color', c1{3}, 'LineWidth', 3);
h(4) = plot(X, Ya, '-d', 'MarkerSize', 10, 'Color', c1{4}, 'LineWidth', 3);
vec = X+i(1)*log2(X);
h(5) = plot(X, Yo(md)+vec-vec(md), '--', 'Color', c2{i(1)+1}, 'LineWidth', 2);
vec = X+i(2)*log2(X);
h(6) = plot(X, Yr(md)+vec-vec(md), '--', 'Color', c2{i(2)+1}, 'LineWidth', 2);
vec = X+i(3)*log2(X);
h(7) = plot(X, Yf(md)+vec-vec(md), '--', 'Color', c2{i(3)+1}, 'LineWidth', 2);
vec = X+i(4)*log2(X);
h(8) = plot(X, Ya(md)+vec-vec(md), '--', 'Color', c2{i(4)+1}, 'LineWidth', 2);

% lgd = legend(h,{'Path Time','Recovery Time','Factorization Time',...
%     'Application Time','N','N log_2(N)','N log_2^2(N)'},...
%     'Location','eastoutside');

axis tight;
axis square;
xlabel('log_2(N)'); ylabel('log_2(\cdot)');
set(gca, 'FontSize', 16);
b = get(gca);
set(b.XLabel, 'FontSize', 16); set(b.YLabel, 'FontSize', 16);
set(b.ZLabel, 'FontSize', 16); set(b.Title, 'FontSize', 16);
hold off;

saveas(gcf,['./results/' F '.pdf'])

close(pic);

end
