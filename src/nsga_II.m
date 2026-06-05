% =========================================================================
% nsga_II.m
%
% NSGA-II MULTI-OBJECTIVE OPTIMISATION — True Pareto Front
%
% Finds the TRUE Pareto front by optimising directly on mission metrics:
%   Objective 1: Settling time      
%   Objective 2: RMS control effort 
%
% Unlike the LQR alpha-sweep, NSGA-II searches over all
% possible Q and R diagonal entries 
%
% Design variables (12, in log10 scale):
%   x(1:4)   — Q diagonal for plant states [x, y, xdot, ydot]
%   x(5:10)  — Q diagonal for servo states
%   x(11:12) — R diagonal for control inputs [u1, u2]
%
% =========================================================================

clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

%% =========================================================================
%  SYSTEM SETUP
%% =========================================================================
[A, B, C, D, E]      = get_system_matrices();
w1                    = 1;
[G1, G2, Aa, Ba, Ea] = get_servo_matrices(A, B, C, E, w1);

p_obs = [-20+2j, -20-2j, -25+2j, -25-2j];
L     = place(A', C', p_obs)';

T             = 30;
w0            = [1; 1; 0; 0];
w_hat0        = [0; 0; 0; 0];
settle_thresh = 0.05;

% Store in global for objective function access
global CR3BP_SYS
CR3BP_SYS.A  = A;  CR3BP_SYS.B  = B;  CR3BP_SYS.C  = C;
CR3BP_SYS.E  = E;  CR3BP_SYS.G1 = G1; CR3BP_SYS.G2 = G2;
CR3BP_SYS.Aa = Aa; CR3BP_SYS.Ba = Ba; CR3BP_SYS.L  = L;
CR3BP_SYS.T  = T;  CR3BP_SYS.w0 = w0; CR3BP_SYS.w_hat0 = w_hat0;
CR3BP_SYS.settle_thresh = settle_thresh;
CR3BP_SYS.w1 = w1;

%% =========================================================================
%  LOAD LQR AND POLE PLACEMENT REFERENCE DATA
%% =========================================================================
if exist('../results/pareto_data.mat', 'file')
    pareto = load('../results/pareto_data.mat');
    fprintf('Loaded LQR Pareto data (%d points).\n\n', ...
            length(pareto.rms_u_valid));
else
    pareto = [];
    fprintf('WARNING: pareto_data.mat not found. Run script 06 first.\n\n');
end

%% =========================================================================
%  DESIGN VARIABLE BOUNDS (log10 scale)
%% =========================================================================
n_vars = 12;

lb = [log10(0.1),  log10(0.1),  log10(0.1),  log10(0.1),  ...  % q_plant
      log10(0.01), log10(0.01), log10(0.01), ...                % q_servo
      log10(0.01), log10(0.01), log10(0.01), ...
      log10(1e-6), log10(1e-6)];                                 % r_ctrl

ub = [log10(1000), log10(1000), log10(1000), log10(1000), ...   % q_plant
      log10(100),  log10(100),  log10(100),  ...                 % q_servo
      log10(100),  log10(100),  log10(100),  ...
      log10(10),   log10(10)];                                   % r_ctrl

%% =========================================================================
%  NSGA-II OPTIMISATION
%% =========================================================================
fprintf('=== NSGA-II Multi-Objective Optimisation ===\n');
fprintf('Design variables : %d (Q and R diagonal entries, log10 scale)\n', n_vars);
fprintf('Objectives       : settling time + RMS control effort\n');
fprintf('Population       : 100,  Generations: 100\n');

options = optimoptions('gamultiobj', ...
    'PopulationSize',    100,   ...
    'MaxGenerations',    100,   ...
    'ParetoFraction',    0.5,   ...
    'CrossoverFraction', 0.8,   ...
    'FunctionTolerance', 1e-4,  ...
    'Display',           'iter',...
    'PlotFcn',           @gaplotpareto);

t_start = tic;
[x_opt, f_opt, exitflag, output] = gamultiobj( ...
    @cr3bp_objectives, n_vars, [], [], [], [], lb, ub, options);
t_elapsed = toc(t_start);

fprintf('\n=== Optimisation Complete ===\n');
fprintf('Time       : %.1f minutes\n', t_elapsed/60);
fprintf('Generations: %d\n', output.generations);
fprintf('Solutions  : %d on Pareto front\n', size(f_opt,1));

%% =========================================================================
%  PROCESS RESULTS
%% =========================================================================

% Filter penalty solutions
valid = f_opt(:,1) < 90  & f_opt(:,2) < 90;
nsga2_settle_raw = f_opt(valid, 1);
nsga2_rms_raw    = f_opt(valid, 2);
x_valid          = x_opt(valid, :);

fprintf('Valid solutions: %d / %d\n', sum(valid), size(f_opt,1));
fprintf('Settling range : %.2f to %.2f s\n', min(nsga2_settle_raw), max(nsga2_settle_raw));
fprintf('RMS range           : %.2f to %.2f\n', min(nsga2_rms_raw),    max(nsga2_rms_raw));

% Sort by settling time so front goes from slow/cheap to fast/expensive — matches LQR front orientation
[nsga2_settle, sort_idx] = sort(nsga2_settle_raw, 'descend');
nsga2_rms  = nsga2_rms_raw(sort_idx);
x_sorted   = x_valid(sort_idx, :);

%% =========================================================================
%  KNEE POINT DETECTION
%  Front is now ordered: settle decreasing, rms increasing
%% =========================================================================
n_pts      = length(nsga2_settle);
efficiency = zeros(1, n_pts-1);
for ii = 1:n_pts-1
    d_settle = nsga2_settle(ii) - nsga2_settle(ii+1);   % settling improvement
    d_rms    = nsga2_rms(ii+1)  - nsga2_rms(ii);        % extra RMS cost
    if d_rms > 1e-6
        efficiency(ii) = d_settle / d_rms;
    end
end

% Knee detection using cumulative RMS change rate
% Normalise axes to [0,1] then find point of maximum curvature
% using the perpendicular distance from the diagonal of the bounding box
rms_min = min(nsga2_rms);   rms_max = max(nsga2_rms);
st_min  = min(nsga2_settle); st_max  = max(nsga2_settle);

rms_n = (nsga2_rms    - rms_min) / (rms_max - rms_min + 1e-10);
st_n  = (nsga2_settle - st_min)  / (st_max  - st_min  + 1e-10);

% Only consider the first 60% of points by RMS (exclude high-effort flat tail)
rms_cutoff = rms_min + 0.60 * (rms_max - rms_min);
search_pts = find(nsga2_rms <= rms_cutoff);
if isempty(search_pts)
    search_pts = 1:n_pts;
end

% Within search region, find point furthest from diagonal
p1 = [rms_n(search_pts(1)),   st_n(search_pts(1))];
p2 = [rms_n(search_pts(end)), st_n(search_pts(end))];
dv = p2 - p1;

dist = zeros(1, length(search_pts));
for ii = 1:length(search_pts)
    idx_i = search_pts(ii);
    dist(ii) = abs(dv(1)*(p1(2)-st_n(idx_i)) - (p1(1)-rms_n(idx_i))*dv(2)) / ...
               (norm(dv) + 1e-10);
end
[~, local_knee] = max(dist);
knee_idx = search_pts(local_knee);

nsga2_knee_settle = nsga2_settle(knee_idx);
nsga2_knee_rms    = nsga2_rms(knee_idx);
nsga2_knee_x      = x_sorted(knee_idx, :);

fprintf('\nNSGA-II Knee point:\n');
fprintf('  Settling time = %.2f s\n', nsga2_knee_settle);
fprintf('  RMS effort    = %.2f\n',   nsga2_knee_rms);

%% =========================================================================
%  SAVE
%% =========================================================================
if ~exist('../results','dir'); mkdir('../results'); end
save('../results/nsga2_results.mat', ...
     'nsga2_settle', 'nsga2_rms', 'x_sorted', ...
     'nsga2_knee_settle', 'nsga2_knee_rms', 'nsga2_knee_x', ...
     'output', 'exitflag', 't_elapsed', 'valid');
fprintf('\nResults saved to results/nsga2_results.mat\n');

%% =========================================================================
%  FIGURES
%% =========================================================================
if ~exist('../docs/figures','dir'); mkdir('../docs/figures'); end

%% ---- Figure 1: NSGA-II Front alone ------------------------------------
figure('Name','NSGA-II Pareto Front','NumberTitle','off', ...
       'Position',[50 50 900 650]);

plot(nsga2_rms, nsga2_settle, 'r-o', ...
     'LineWidth', 2, 'MarkerSize', 5, 'MarkerFaceColor', 'r', ...
     'DisplayName', 'NSGA-II Pareto Front');
hold on;

scatter(nsga2_knee_rms, nsga2_knee_settle, 260, [0.7 0 0], 'h', ...
        'LineWidth', 2.5, 'MarkerFaceColor', [0.7 0 0], ...
        'DisplayName', sprintf('Knee Point (%.1fs, %.1f)', ...
        nsga2_knee_settle, nsga2_knee_rms));

text(nsga2_knee_rms + (max(nsga2_rms)-min(nsga2_rms))*0.02, ...
     nsga2_knee_settle, ...
     sprintf('  Knee\n  (%.1fs, %.1f)', nsga2_knee_settle, nsga2_knee_rms), ...
     'Color', [0.7 0 0], 'FontSize', 9);

xlabel('RMS Control Effort  \surd(mean(u_1^2 + u_2^2))', 'FontSize', 12);
ylabel('Settling Time (s)', 'FontSize', 12);
title({'NSGA-II True Pareto Front', ...
       'Optimised directly on settling time and RMS control effort'}, ...
      'FontSize', 12);
legend('Location', 'best', 'FontSize', 10);
grid on;

saveas(gcf, '../docs/figures/09a_nsga2_pareto_front.png');
fprintf('Saved: 09a_nsga2_pareto_front.png\n');

%% ---- Figure 2: NSGA-II vs LQR comparison ------------------------------
if ~isempty(pareto)
    % Reload fresh to avoid any workspace variable shadowing
    p = load('../results/pareto_data.mat');

    figure('Name','NSGA-II vs LQR','NumberTitle','off','Position',[50 50 1000 700]);

    h1 = plot(p.rms_u_valid, p.settle_valid, 'b-o', ...
              'LineWidth',2,'MarkerSize',5,'MarkerFaceColor','b', ...
              'DisplayName','LQR Pareto Front (\alpha sweep)');
    hold on;

    h2 = plot(nsga2_rms, nsga2_settle, 'r-o', ...
              'LineWidth',2,'MarkerSize',5,'MarkerFaceColor','r', ...
              'DisplayName','NSGA-II True Pareto Front');

    h3 = scatter(p.pp_results(1).rms_u, p.pp_results(1).settle, 200, ...
                 [0.4 0 0],'s','LineWidth',2.5, ...
                 'DisplayName',sprintf('PP Set 1 (%.1fs, %.1f)', ...
                 p.pp_results(1).settle, p.pp_results(1).rms_u));

    h4 = scatter(p.pp_results(2).rms_u, p.pp_results(2).settle, 200, ...
                 [0.7 0 0.7],'s','LineWidth',2.5, ...
                 'DisplayName',sprintf('PP Set 2 (%.1fs, %.1f)', ...
                 p.pp_results(2).settle, p.pp_results(2).rms_u));

    h5 = scatter(p.bryson_u, p.bryson_settle, 220, [0 0.6 0],'^', ...
                 'LineWidth',2.5,'MarkerFaceColor',[0 0.6 0], ...
                 'DisplayName',sprintf('LQR Bryson (%.1fs, %.1f)', ...
                 p.bryson_settle, p.bryson_u));

    h6 = scatter(p.knee_u, p.knee_settle, 220, [1 0.5 0],'p', ...
                 'LineWidth',2.5,'MarkerFaceColor',[1 0.5 0], ...
                 'DisplayName',sprintf('LQR Knee (%.1fs, %.1f)', ...
                 p.knee_settle, p.knee_u));

    h7 = scatter(nsga2_knee_rms, nsga2_knee_settle, 260, [0.7 0 0],'h', ...
                 'LineWidth',2.5,'MarkerFaceColor',[0.7 0 0], ...
                 'DisplayName',sprintf('NSGA-II Knee (%.1fs, %.1f)', ...
                 nsga2_knee_settle, nsga2_knee_rms));

    xlabel('RMS Control Effort  \surd(mean(u_1^2 + u_2^2))','FontSize',12);
    ylabel('Settling Time (s)','FontSize',12);
    title({'NSGA-II vs LQR Pareto Front Comparison', ...
           'NSGA-II optimises Q/R directly — not restricted to \alpha scaling'}, ...
          'FontSize',11);
    legend([h1 h2 h3 h4 h5 h6 h7],'Location','best','FontSize',9);
    grid on;

    saveas(gcf,'../docs/figures/09b_nsga2_vs_lqr_comparison.png');
    fprintf('Saved: 09b_nsga2_vs_lqr_comparison.png\n');
end

%% =========================================================================
%  SUMMARY TABLE
%% =========================================================================
fprintf('\n=== Performance Comparison Summary ===\n');
fprintf('%-32s %-14s %-12s\n', 'Controller', 'Settling (s)', 'RMS Effort');
fprintf('%s\n', repmat('-', 1, 60));
if ~isempty(pareto)
    fprintf('%-32s %-14.2f %-12.2f\n', 'PP Set 1 (Conservative)', ...
            pareto.pp_results(1).settle, pareto.pp_results(1).rms_u);
    fprintf('%-32s %-14.2f %-12.2f\n', 'PP Set 2 (Aggressive)', ...
            pareto.pp_results(2).settle, pareto.pp_results(2).rms_u);
    fprintf('%-32s %-14.2f %-12.2f\n', 'LQR Bryson (alpha=1)', ...
            pareto.bryson_settle, pareto.bryson_u);
    fprintf('%-32s %-14.2f %-12.2f\n', 'LQR Knee Point', ...
            pareto.knee_settle, pareto.knee_u);
end
fprintf('%-32s %-14.2f %-12.2f\n', 'NSGA-II Knee Point', ...
        nsga2_knee_settle, nsga2_knee_rms);
fprintf('%-32s %-14.2f %-12.2f\n', 'NSGA-II Best Settling', ...
        min(nsga2_settle), nsga2_rms(end));
fprintf('%-32s %-14.2f %-12.2f\n', 'NSGA-II Best RMS', ...
        nsga2_settle(1), min(nsga2_rms));

fprintf('\n=== Key Finding ===\n');
if ~isempty(pareto)
    if min(nsga2_settle) < min(pareto.settle_valid)
        fprintf('NSGA-II achieves %.2fs settling vs LQR minimum %.2fs\n', ...
                min(nsga2_settle), min(pareto.settle_valid));
        fprintf('-> Controllers outside LQR family outperform LQR on settling time.\n');
    end
end
fprintf('All figures saved to docs/figures/\n');

%% =========================================================================
%  OBJECTIVE FUNCTION
%% =========================================================================
function f = cr3bp_objectives(x)
    global CR3BP_SYS

    q_plant = 10.^x(1:4);
    q_servo = 10.^x(5:10);
    r_ctrl  = 10.^x(11:12);

    Q = blkdiag(diag(q_plant), diag(q_servo));
    R = diag(r_ctrl);

    PENALTY = [100, 100];

    try
        Ka = lqr(CR3BP_SYS.Aa, CR3BP_SYS.Ba, Q, R);
    catch
        f = PENALTY; return;
    end

    cl_eigs = eig(CR3BP_SYS.Aa - CR3BP_SYS.Ba * Ka);
    if any(real(cl_eigs) >= 0)
        f = PENALTY; return;
    end

    max_eig = max(abs(real(cl_eigs)));
    dt      = min(0.005, 0.5 / max_eig);
    dt      = max(dt, 1e-5);

    try
        [w, ~, ~, u, t] = run_simulation( ...
            CR3BP_SYS.A,  CR3BP_SYS.B,  CR3BP_SYS.C, ...
            CR3BP_SYS.E,  CR3BP_SYS.G1, CR3BP_SYS.G2, ...
            Ka(:,1:4), Ka(:,5:10), CR3BP_SYS.L, ...
            CR3BP_SYS.w0, CR3BP_SYS.w_hat0, ...
            dt, CR3BP_SYS.T, CR3BP_SYS.w1);
    catch
        f = PENALTY; return;
    end

    if any(~isfinite(w(:))) || any(~isfinite(u(:)))
        f = PENALTY; return;
    end

    settle = compute_settling_time(t, w, CR3BP_SYS.settle_thresh);
    rms_u  = sqrt(mean(u(1,:).^2 + u(2,:).^2));

    if settle >= CR3BP_SYS.T * 0.95
        f = PENALTY; return;
    end

    % Soft penalty for excessive control effort
    max_u = max(abs(u(:)));
    if max_u > 400
        rms_u = rms_u * (1 + (max_u - 400)/400);
    end

    f = [settle, rms_u];
end