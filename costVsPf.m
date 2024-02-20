clearvars; close all; clc;

% set of linearly spaced scalings of initial to final cov
Pscales = linspace(0.05, 1.5, 500); % in percent/100

% set initial cov
opts.P0 = diag([0.1; 0.1; 0.02; 0.02]);
S0      = chol(opts.P0);


% loop through convex solving and get cost
cost = zeros(length(Pscales),1);
for nn = 1:length(Pscales)
    ps = Pscales(nn);
    Sf = ps.*S0;
    opts.Pf = Sf*Sf';

    % run optimization
    LGMConvexDriver;

    % save cost if problem was solved (not infeasible/reachable)
    if cvx_infs
        continue
    else
        cost(nn) = cvx_optval;
    end

    % clear all vars except opts
    clearvars -except opts cost Pscales S0
end

% plot cost vs percent scaling
f = figure(102);
f.Position = [300, 200, 900, 800];
plot(Pscales, cost, 'k-', 'LineWidth', 0.8);
xline(0.101, 'r--', 'LineWidth', 1.4);
xlabel('\% $S_f$ / $S_0$'); ylabel('Cost (OL \& CL)');

