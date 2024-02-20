if ~exist("solstruct", "var")
    rundir  = 'LGMConvex/results/';
    runname = 'N20_gaussianKLD';   % run to load
    load([rundir runname '/solstruct.mat']);
end
opts = solstruct.opts;

saveplt = true;

% plotting covariance
theta = linspace(0, 2*pi, 100);

% plot x-y trajectory
f = figure(1);
f.Position = [300, 200, 900, 800];
plot(solstruct.opt_m(1,:), solstruct.opt_m(2,:), ...
     '+--', 'Color', 'black', 'DisplayName', '$m_x$'); hold on;
xlabel('x-position [DU]');
ylabel('y-position [DU]');

% plot final solved for covariance
xt   = cos(theta);
yt   = sin(theta);
Srr  = sqrt(11.98)*chol(solstruct.opt_P(1:2,1:2,end));
epts = [xt;yt]' * squeeze(Srr);
epts = bsxfun(@plus, epts, solstruct.opt_m(1:2,end)');
plot(epts(:,1), epts(:,2), 'k-', 'LineWidth', 2, 'DisplayName', '$P_{xx,N}$');

% plot final prescribed covariance 
Srr  = sqrt(11.98)*chol(opts.Pf(1:2,1:2));
epts = [xt;yt]' * squeeze(Srr);
epts = bsxfun(@plus, epts, opts.xf(1:2)');
plot(epts(:,1), epts(:,2), 'g--', 'LineWidth', 2, 'DisplayName', '$P_f$');
legend('Location','best'); legend('AutoUpdate','off');

% plot solved for covariance at various times
for k = 1:3:size(solstruct.opt_P,3)
    Srr  = sqrt(11.98)*chol(solstruct.opt_P(1:2,1:2,k));
    epts = [xt;yt]' * squeeze(Srr);
    epts = bsxfun(@plus, epts, solstruct.opt_m(1:2,k)');
    
    figure(1);
    plot(epts(:,1), epts(:,2), 'k');
end

if saveplt
    toTikz([rundir runname '/pos_opt_traj.tex']);
end


% plot control
f = figure(2);
f.Position = [1200, 200, 900, 800];
plot(solstruct.u_OL(1,:), solstruct.u_OL(2,:), ...
     '+--', 'Color', 'black', 'DisplayName', '$u_{OL}$'); hold on;
xlabel('$u_x$'); ylabel('$u_y$');
xlim([-2, 2]); ylim([-2, 2])

% plot control covariance
for k = 1:size(solstruct.Puu,3)
    Srr  = sqrt(11.98)*chol(solstruct.Puu(1:2,1:2,k));
    epts = [xt;yt]' * squeeze(Srr);
    epts = bsxfun(@plus, epts, solstruct.u_OL(1:2,k)');
    
    figure(2);
    plot(epts(:,1), epts(:,2), 'k');
end

% plot control bound
plot(2*xt, 2*yt, 'r--'); axis equal;

if saveplt
    toTikz([rundir runname '/cntrl_opt_traj.tex']);
end