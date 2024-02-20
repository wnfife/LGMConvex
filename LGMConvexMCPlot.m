if ~exist("solstruct", "var")
    rundir  = 'LGMConvex/results/';
    runname = 'N20_gaussianRidder';   % run to load
    load([rundir runname '/solstruct.mat']);
    load([rundir runname '/MC/MC_stats.mat']);
end

mcdir = [rundir runname '/MC/'];

opts = solstruct.opts;

saveplt = false;

% histogram of final marginals
f1          = figure(1);
f1.Position = [600, 200, 900, 800];
tiledlayout(2,2);
Nx   = 800;
xlbs = {'$x$ [DU]', '$y$ [DU]', ...
        '$\dot{x}$ [VU]', '$\dot{y}$ [VU]'};
XB(:,1) = zeros(opts.n,1) - 3.2*sqrt(diag(opts.Pf));
XB(:,2) = zeros(opts.n,1) + 3.2*sqrt(diag(opts.Pf));

XLIMS = [-3.03e-3, 3.03e-3; ...
         -5.4e-3, 5.4e-3; ...
         -3.2e-3, 3.2e-3; ...
         -2.4e-3, 2.4e-3; ...
         -1.7e-3, 1.7e-3; ...
         -1.1e-3, 1.1e-3;];

for i = 1:4
    nexttile;
    ei = reshape(MC_stats.MC_x.all(end,:,i), [], 1);
    histogram(ei, 70, ...
              'FaceColor', 'black', ...
              'Normalization', 'pdf'); hold on;
    xlabel(xlbs{i});
    %xlim(XLIMS(i,:));

    %cntout = sum(ei < -sqrt(10.83)*sqrt(opts.Pf(i,i))) + sum(ei > sqrt(10.83)*sqrt(opts.Pf(i,i)));
    %disp([num2str(100 - 100*(cntout/50000)) '% in 99.9% (' xlbs{i} ')']);

    Si = sqrt(opts.Pf(i,i));

    xline(solstruct.opt_m(i,end) - 3*std(ei), 'k--'); xline(solstruct.opt_m(i,end) + 3*std(ei), 'k--');
    xline(opts.xf(i) - 3*Si, 'r--'); xline(opts.xf(i) + 3*Si, 'r--')
    
    xt = linspace(opts.xf(i) - 3.5*Si, opts.xf(i) + 3.5*Si, Nx);
    yt = normpdf(xt, opts.xf(i), Si);
    plot(xt, yt, 'r-', 'LineWidth', 1);
end

if saveplt
    cleanfigure;
    toTikz([mcdir '/MC_marginals.tex']);
end


% plot x-y CL trajectories
f2          = figure(2);
f2.Position = [700, 200, 900, 800];
plot(MC_stats.MC_x.all(:,1:100:end,1), MC_stats.MC_x.all(:,1:100:end,2), 'Color', [rgb('DarkSlateGray'), 0.2] ); hold on;
xlabel('x-position [DU]'); ylabel('y-position [DU]');

% plotting covariance
theta = linspace(0, 2*pi, 100);
xt    = cos(theta);
yt    = sin(theta);

% plot solved for covariance at various times
for k = 1:3:size(solstruct.opt_P,3)
    Srr  = sqrt(11.98)*chol(solstruct.opt_P(1:2,1:2,k));
    epts = [xt;yt]' * squeeze(Srr);
    epts = bsxfun(@plus, epts, solstruct.opt_m(1:2,k)');
    
    figure(2);
    plot(epts(:,1), epts(:,2), 'k');
end

% plot final prescribed covariance 
Srr  = sqrt(11.98)*chol(opts.Pf(1:2,1:2));
epts = [xt;yt]' * squeeze(Srr);
epts = bsxfun(@plus, epts, opts.xf(1:2)');
plot(epts(:,1), epts(:,2), 'g--', 'LineWidth', 2);

if saveplt
    toTikz([mcdir '/MC_xCL.tex']);
end

% plot x-y CL filtered trajectories
f201          = figure(201);
f012.Position = [700, 200, 900, 800];
plot(MC_m(:,1:100:end,1), MC_m(:,1:100:end,2), 'Color', [rgb('DarkSlateGray'), 0.2] ); hold on;
xlabel('x-position [DU]'); ylabel('y-position [DU]');
% plotting covariance
theta = linspace(0, 2*pi, 100);
xt    = cos(theta);
yt    = sin(theta);
% plot solved for covariance at various times
for k = 1:3:size(solstruct.opt_P,3)
Srr  = sqrt(11.98)*chol(solstruct.opt_P(1:2,1:2,k));
epts = [xt;yt]' * squeeze(Srr);
epts = bsxfun(@plus, epts, solstruct.opt_m(1:2,k)');
figure(201);
plot(epts(:,1), epts(:,2), 'k');
end
% plot final prescribed covariance
Srr  = sqrt(11.98)*chol(opts.Pf(1:2,1:2));
epts = [xt;yt]' * squeeze(Srr);
epts = bsxfun(@plus, epts, opts.xf(1:2)');
plot(epts(:,1), epts(:,2), 'g--', 'LineWidth', 2);


% closed loop control
uCL = reshape(vecnorm(MC_stats.MC_uCL.all, 2, 3), opts.N, []);
f4 = figure(4);
f4.Position = [800, 200, 900, 800];

plot(uCL(:,1:100:end), 'Color', [rgb('DarkSlateGray'), 0.2], 'LineWidth', 0.4);
hold on;
yline(opts.umax, 'r--', 'LineWidth', 1.5);
xlabel('Node'); ylabel('$\| u_k \|$ [AU]'); hold on;
plot(vecnorm( solstruct.u_OL, 2, 1 ), 'k-', 'LineWidth', 1.5);

if saveplt
    cleanfigure;
    toTikz([mcdir '/MC_uCL.tex']);
end
