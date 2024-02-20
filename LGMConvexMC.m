clearvars; close all; clc;

% load case solution struct
rundir  = 'LGMConvex/results/';
runname = 'N20_gaussianRidder';
load([rundir runname '/solstruct.mat']);
opts = solstruct.opts;

% make MC directory
mcdir = [rundir runname '/MC/'];
mkdir(mcdir);

% number of samples
Nmc = 5000;

% estimation parameters
Hx  = [eye(2), zeros(2,2)];  % measurement of position
Svv = 0.05.*eye(2);          % measurement noise SRF

% storage
MC_x.all   = zeros(opts.N, Nmc, opts.n);
MC_xOL.all = zeros(opts.N, Nmc, opts.n);
MC_uCL.all = zeros(opts.N-1, Nmc, opts.m);
MC_m       = zeros(opts.N, Nmc, opts.n);

% sample GM
Xsamples = sampleGM(opts.GM.w, opts.GM.m, opts.GM.P, Nmc);

Sww = chol(opts.Pww)';

% start MC loop
for i = 1:Nmc
    % sample from GM
    x0i = Xsamples(:,i);

    % start sim loop
    MC_x.all(1,i,:)   = x0i;
    MC_xOL.all(1,i,:) = x0i;
    xkm1   = x0i;
    xkm1OL = x0i;
    yk     = x0i - solstruct.opt_m(:,1);
    
    [mkm1, Pkm1] = GMmeancov(opts.GM.w, opts.GM.m, opts.GM.P);
    MC_m(1,i,:)  = mkm1;

    ykm = mkm1 - solstruct.opt_m(:,1);
    for k = 1:opts.N
        % update mean
        Pkm = opts.Fx*Pkm1*opts.Fx' + opts.Fw*opts.Pww*opts.Fw';
        mkm = opts.Fx*mkm1;
        zk  = Hx*mkm + Svv*randn(2,1);
        K   = Pkm*Hx'/(Hx*Pkm*Hx' + Svv*Svv');
        mkp = mkm + K*(zk - Hx*mkm);
        M   = eye(4) - K*Hx;
        Pkp = M*Pkm*M' + K*(Svv*Svv')*K';

        % extract control
        u_OL = solstruct.u_OL(:,k);
        u_CL = u_OL + solstruct.U_gain(:,:,k)*yk;
        u_CLm = u_OL + solstruct.U_gain(:,:,k)*ykm;

        % update control feedback term
        yk = opts.Fx*yk;
        ykm = opts.Fx*ykm;

        % save control
        MC_uCL.all(k,i,:) = u_CL;

        % propagate sample
        wkm1 = Sww*randn(opts.n,1);
        xkOL = opts.Fx*xkm1OL + opts.Fu*u_OL + opts.Fw*wkm1;
        xkCL = opts.Fx*xkm1 + opts.Fu*u_CL + opts.Fw*wkm1;
        
        % update mean w control
        mkp = mkp + opts.Fu*u_CLm;

        % save samples
        MC_x.all(k+1, i, :)   = xkCL;
        MC_xOL.all(k+1, i, :) = xkOL;
        MC_m(k+1, i, :)       = mkp;

        % recursion
        xkm1   = xkCL;
        xkm1OL = xkOL;
        mkm1   = mkp;
        Pkm1   = Pkp;
    end
end

% save
MC_stats.MC_x   = MC_x;
MC_stats.MC_xOL = MC_xOL;
MC_stats.MC_uCL = MC_uCL;
save([mcdir 'MC_stats'], 'MC_stats');

LGMConvexMCPlot;
