% ==============================================================
% Linear Gaussian mixture steering convex formulation
% using first order markov feedback
% ==============================================================

clearvars; close all; clc;

%% parameter setup
rundir  = 'LGMConvex/results/';
runname = 'N20_heavyTail';

mkdir([rundir runname '/']);

% dimensions and nodes
opts.n = 4;
opts.m = 2;
opts.N = 20;
nU     = opts.N*opts.m;
nX     = (opts.N+1)*opts.n;

% time step
opts.dt = 0.25;

% initial and final Gaussians
opts.x0 = [2; 4; 3; 2];
opts.xf = [8; 2; 0; 0];
opts.P0 = diag([0.1; 0.1; 0.02; 0.02]);
opts.Pf = diag([0.1; 0.1; 0.01; 0.01]);

% GM approximation of initial Gaussian
[GM.w, GM.m, GM.P] = split2gm(1, opts.x0, opts.P0, [1], 3);

% UNCOMMENT if heavy-tailed Gaussian desired
GM.w(1) = GM.w(1) - 0.1; GM.w(3) = GM.w(3) - 0.1; GM.w = GM.w./sum(GM.w);
GM.m(1,1) = GM.m(1,1) - 0.3;   GM.m(1,3) = GM.m(1,3) + 0.3;
GM.P(1,1,1) = 6.2*GM.P(1,1,1); GM.P(1,1,3) = 6.2*GM.P(1,1,3);
GM.P(1,1,2) = 2*GM.P(1,1,2);

% UNCOMMENT if bimodal distribution (in x) desired
% GM.w = [0.5; 0.5];
% GM.m(:,1) = [1; opts.x0(2:end)];
% GM.m(:,2) = [3; opts.x0(2:end)];
% GM.m(:,3) = []; GM.P(:,:,3) = [];

% plot initial Gaussian and GM (for debugging)
plt_init = true;
if plt_init
    figure(101);
    yt = 0;
    xt = linspace(1, 3, 800);
    for i = 1:length(GM.w)
    yt = yt + GM.w(i)*normpdf(xt, GM.m(1,i), sqrt(GM.P(1,1,i)));
    end
    plot(xt, yt, 'k-');
    yg = normpdf(xt, opts.x0(1), sqrt(opts.P0(1,1)));
    hold on
    plot(xt, yg, 'Color', rgb('DodgerBlue'), 'LineWidth', 0.8);
    legend('GM', 'Gaussian');
    xlabel('$x$ [DU]'); ylabel('$p(x)$');
end

opts.GM = GM;

% mean and cov of initial GM
[mx0, Pxx0] = GMmeancov(GM.w, GM.m, GM.P);

% inverse of weighted sum of GM covariances
Lambda0 = inv(sum(GM.P.*reshape(GM.w,1,1,[]),3));

% SRF of initial covariance
S0 = chol(Pxx0)';

% LTI dynamics matrix
Fx      = eye(opts.n);
Fx(1,3) = opts.dt;
Fx(2,4) = opts.dt;
opts.Fx = Fx;

% constant control mapping matrix 
Fu      = [0, 0; ...
           0, 0; ...
           opts.dt, 0; ...
           0, opts.dt];
opts.Fu = Fu;

% process noise mapping matrix
Fw      = 0.1*eye(opts.n);
opts.Fw = Fw;

% constant process noise covariance
Pww      = 1e-4*eye(opts.n);
opts.Pww = Pww;

% block matrix form of dynamics
Fxbar = eye(opts.n);
Fubar = zeros(opts.n, opts.m*opts.N);
Fwbar = zeros(opts.n, opts.n*opts.N);
for k = 1:opts.N
    Fxbar = [Fxbar; Fx*Fxbar( end-opts.n + 1:end, : )];
    Fubar = [Fubar; Fx*Fubar( end-opts.n+1:end,1:(k-1)*opts.m ) Fu zeros(opts.n,(opts.N-k)*opts.m)];
    Fwbar = [Fwbar; Fw*Fwbar( end-opts.n+1:end,1:(k-1)*opts.n ) Fw zeros(opts.n,(opts.N-k)*opts.n)];
end
Pwwbar = kron(Pww, eye(opts.N));
Swwbar = chol(Pwwbar)';

% control magnitude chance constraint params
opts.beta = (100 - 99.99)/100;        % allowable violation probability
gamma     = sqrt(2*log(1/opts.beta)); % scaling factor in chance constraint
opts.umax = 2;                        % max CL control magnitude

% OL control cost scaling matrix
opts.R = eye(opts.N*opts.m);

% final time mapping matrix for state
EN = zeros(opts.n, nX);
EN(:,end-opts.n+1:end) = eye(opts.n);

% identity matrix of size nX = (N+1)*n
InX = eye(nX);

%% convex optimization
cvx_begin sdp                    % sdp = semidefinite programming
    cvx_solver mosek;            % fast large-scale solver

    % declare OL (V) and closed loop (G) optimization variables
    variable V(nU);          
    variable Garr(opts.m, opts.n, opts.N);
    expression G;
    expression U(opts.N); % closed loop magnitude chance constraint

    % build G (for now, need loop)
    G = blkdiag(Garr(:,:,1), Garr(:,:,2));
    for i = 3:opts.N
        G = blkdiag(G, Garr(:,:,i));
    end
    G = [G, zeros(nU, opts.n)];

    % build CL magnitude chance constraint evaluation
    for i = 1:opts.N
        EU = zeros(opts.m, opts.N*opts.m);  
        EU(:, (i-1)*opts.m + 1: i*opts.m) = eye(opts.m);
        U(i) = norm(EU*V) + gamma*(norm(S0*Fxbar'*G'*EU'));
    end

    % build final state constraints
    PHI = (InX + Fubar*G)*Fxbar;
    mN  = EN*(Fxbar*mx0 + Fubar*V);
    DM  = zeros(opts.n, opts.n);
    for k = 1:length(GM.w)
        DM = GM.w(k)*(GM.m(:,k) - mx0)*(GM.m(:,k) - mx0)';
    end
    C   = opts.Pf - EN*(Fxbar*DM*Fxbar' + Fwbar*Pwwbar*Fwbar')*EN'; % for schur compliment
    X   = [Lambda0, PHI'*EN'; ...
           EN*PHI, C];

    % final covariance square root factor
    Ssum0 = chol(sum(GM.P.*reshape(GM.w,1,1,[]),3))';
    SDM   = zeros(size(DM));
    try
        SDM = chol(DM)';
    catch
        SDM = zeros(size(DM));
    end
    SN = [Ssum0'*PHI'*EN'; SDM'*Fxbar'*EN'; Swwbar'*Fwbar'*EN'];

    % declare objective
    MV  = square_pos(norm( V'*opts.R ));           % norm of mean control
    TG  = square_pos(norm(S0'*Fxbar'*G', 'fro'));  % norm of control covariance
    minimize( TG + MV );
    
    % declare constraints
    subject to
        mN == opts.xf;     % final GM mean == xf
        X >= 0;            % final GM cov <= Pf
        %Scon <= rho_x;        
        U <= opts.umax;    % CL control magnitude (if assumed Gaussian) <= umax at each node
cvx_end

cvx_infs = false;
if strcmp(cvx_status, 'Infeasible')
    disp('INFEASIBLE --- NOT STORING OR PLOTTING');
    cvx_infs = true;
end

%% unpack, save, and plot
Pxxbar = PHI*sum(GM.P.*reshape(GM.w,1,1,[]),3)*PHI' + Fxbar*DM*Fxbar' + Fwbar*Pwwbar*Fwbar';
mxbar  = Fxbar*mx0 + Fubar*V;
Puubar = G*Fxbar*Pxx0*Fxbar'*G';

% decompose block matrices
u_OL_arr = reshape(V, opts.m, []);
G_arr    = zeros(opts.m,opts.n,opts.N);
Puu      = zeros(opts.m,opts.m,opts.N);
for k = 1:opts.N
    G_arr(:,:,k) = G( (k-1)*opts.m+1:k*opts.m, (k-1)*opts.n+1:k*opts.n );
    EU = zeros(opts.m, opts.N*opts.m);
    EU(:, (k-1)*opts.m + 1: k*opts.m) = eye(opts.m);
    Puu(:,:,k)   = EU*Puubar*EU';
end

% solution structure for saving and plotting
solstruct.opt_m  = zeros(opts.n, opts.N+1);
solstruct.opt_P  = zeros(opts.n, opts.n, opts.N+1);
solstruct.u_OL   = u_OL_arr;
solstruct.U_gain = G_arr;
solstruct.Puu    = Puu;
solstruct.opts   = opts;

solstruct.opt_m(:,1)   = mx0;
solstruct.opt_P(:,:,1) = Pxx0;
for k = 1:opts.N+1
    EX = zeros(opts.n, nX);
    EX(:, (k-1)*opts.n + 1: k*opts.n) = eye(opts.n);
    solstruct.opt_m(:,k)   = EX*mxbar;
    solstruct.opt_P(:,:,k) = EX*Pxxbar*EX';
end


if ~cvx_infs
    save( [rundir runname '/solstruct.mat'], "solstruct" );
    LGMConvexPlot;
end



