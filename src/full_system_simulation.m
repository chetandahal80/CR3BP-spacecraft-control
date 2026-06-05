% =========================================================================
% 08_full_system_simulation.m
%
% FULL CLOSED-LOOP SYSTEM SIMULATION AND COMPARISON
%
% Controllers compared:
%   1. Pole Placement Set 1  (conservative)
%   2. Pole Placement Set 2  (aggressive)
%   3. LQR — Bryson's Method (alpha=1)
%   4. LQR — Knee Point      (alpha scaling - optimal tradeoff from Pareto front)
%   5. LQR — NSGA-II Knee    (optimal tradeoff from NSGA-II optimisation)
%
% =========================================================================

clear; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

% System matrices
[A, B, C, D, E]      = get_system_matrices();
w1                    = 1;
[G1, G2, Aa, Ba, Ea] = get_servo_matrices(A, B, C, E, w1);

% Observer gain
p_obs = [-20+2j, -20-2j, -25+2j, -25-2j];
L     = place(A', C', p_obs)';

% Bryson's Q0, R0 baseline
x_max    = 0.5;  y_max    = 0.5;
xdot_max = 0.5;  ydot_max = 0.5;
u1_max   = 200;  u2_max   = 200;

Q0 = blkdiag(diag([1/x_max^2, 1/y_max^2, 1/xdot_max^2, 1/ydot_max^2]), ...
             1.0*eye(6));
R0 = diag([1/u1_max^2, 1/u2_max^2]);

%% =========================================================================
%  LOAD PARETO DATA
%% =========================================================================
if exist('../results/pareto_data.mat', 'file')
    pareto       = load('../results/pareto_data.mat');
    knee_alpha   = pareto.knee_alpha;
    knee_settle  = pareto.knee_settle;
    knee_u       = pareto.knee_u;
    rms_u_valid  = pareto.rms_u_valid;
    settle_valid = pareto.settle_valid;
    pp_results   = pareto.pp_results;
    fprintf('=== Pareto Data Loaded ===\n');
    fprintf('  Knee alpha         = %.6f\n', knee_alpha);
    fprintf('  Knee settling time = %.4f s\n', knee_settle);
    fprintf('  Knee RMS effort    = %.4f\n\n', knee_u);
else
    error('pareto_data.mat not found. Run script 06 first.');
end

%% =========================================================================
%  LOAD NSGA-II DATA
%% =========================================================================
nsga2_available = false;
if exist('../results/nsga2_results.mat', 'file')
    nsga2 = load('../results/nsga2_results.mat');
    nsga2_available = true;
    fprintf('=== NSGA-II Data Loaded ===\n');
    fprintf('  NSGA-II Knee settling time = %.4f s\n', nsga2.nsga2_knee_settle);
    fprintf('  NSGA-II Knee RMS effort    = %.4f\n\n', nsga2.nsga2_knee_rms);
else
    fprintf('WARNING: nsga2_results.mat not found. Run script 09 first.\n');
    fprintf('Continuing without NSGA-II controller.\n\n');
end

%% =========================================================================
%  DEFINE CONTROLLERS
%% =========================================================================
poles_1 = [-2.1587, -0.1, -0.5+1.8626j, -0.5-1.8626j, ...
            -2+2j, -2-2j, -3+1j, -3-1j, -4, -5];
poles_2 = [-3+3j, -3-3j, -4+4j, -4-4j, ...
            -3+1j, -3-1j, -4+2j, -4-2j, -5, -6];

Ka_pp1 = place(Aa, Ba, poles_1);
Ka_pp2 = place(Aa, Ba, poles_2);

% Load pre-computed LQR gains and dt 
if isfield(pareto, 'Ka_knee_saved') && isfield(pareto, 'Ka_bryson_saved')
    Ka_bryson     = pareto.Ka_bryson_saved;
    Ka_knee       = pareto.Ka_knee_saved;
    dt_lqr_bryson = pareto.dt_bryson;
    dt_lqr_knee   = pareto.dt_knee;
    fprintf('  Loaded LQR gains: dt_bryson=%.2e, dt_knee=%.2e\n\n', ...
            dt_lqr_bryson, dt_lqr_knee);
else
    Ka_bryson = lqr(Aa, Ba, Q0, R0);
    Ka_knee   = lqr(Aa, Ba, knee_alpha*Q0, (1/knee_alpha)*R0);
    cl_b = eig(Aa - Ba*Ka_bryson);
    dt_lqr_bryson = max(1e-5, min(0.005, 0.5/max(abs(real(cl_b)))));
    cl_k = eig(Aa - Ba*Ka_knee);
    dt_lqr_knee   = max(1e-5, min(0.005, 0.5/max(abs(real(cl_k)))));
end

% Build controller list
controllers  = {Ka_pp1, Ka_pp2, Ka_bryson, Ka_knee};
ctrl_labels  = {'PP Set 1 (Conservative)', ...
                'PP Set 2 (Aggressive)', ...
                'LQR — Bryson''s Method', ...
                sprintf('LQR — Knee Point (\\alpha=%.3f)', knee_alpha)};
colors       = {'b', 'r', [0 0.6 0], [1 0.5 0]};
lstyles      = {'-', '-', '-', '--'};

cl_pp1 = eig(Aa - Ba*Ka_pp1);
dt_pp1 = max(1e-5, min(0.005, 0.5/max(abs(real(cl_pp1)))));
cl_pp2 = eig(Aa - Ba*Ka_pp2);
dt_pp2 = max(1e-5, min(0.005, 0.5/max(abs(real(cl_pp2)))));
dt_list = {dt_pp1, dt_pp2, dt_lqr_bryson, dt_lqr_knee};

% Add NSGA-II knee controller if available
if nsga2_available
    x_knee   = nsga2.nsga2_knee_x;
    q_plant  = 10.^x_knee(1:4);
    q_servo  = 10.^x_knee(5:10);
    r_ctrl   = 10.^x_knee(11:12);
    Q_nsga   = blkdiag(diag(q_plant), diag(q_servo));
    R_nsga   = diag(r_ctrl);
    Ka_nsga  = lqr(Aa, Ba, Q_nsga, R_nsga);
    cl_nsga  = eig(Aa - Ba*Ka_nsga);
    dt_nsga  = max(1e-5, min(0.005, 0.5/max(abs(real(cl_nsga)))));

    controllers{end+1} = Ka_nsga;
    ctrl_labels{end+1} = sprintf('NSGA-II Knee (%.1fs, %.1f)', ...
                                  nsga2.nsga2_knee_settle, nsga2.nsga2_knee_rms);
    colors{end+1}      = [0.6 0 0.6];
    lstyles{end+1}     = ':';
    dt_list{end+1}     = dt_nsga;
end

n_ctrl = length(controllers);

%% =========================================================================
%  SIMULATE ALL CONTROLLERS
%% =========================================================================
T      = 30;
w0     = [1; 1; 0; 0];
w_hat0 = [0; 0; 0; 0];

results = cell(n_ctrl, 1);
u_store = cell(n_ctrl, 1);
t_store = cell(n_ctrl, 1);

fprintf('Simulating %d controllers...\n', n_ctrl);
for m = 1:n_ctrl
    Ka_m = controllers{m};
    Kx_m = Ka_m(:, 1:4);
    Ki_m = Ka_m(:, 5:10);
    dt_m = dt_list{m};

    [w_m, ~, ~, u_m, t_m] = run_simulation(A, B, C, E, G1, G2, ...
                                             Kx_m, Ki_m, L, ...
                                             w0, w_hat0, dt_m, T, w1);
    results{m} = w_m;
    u_store{m} = u_m;
    t_store{m} = t_m;
    fprintf('  %s — done (dt=%.2e)\n', ctrl_labels{m}, dt_m);
end
fprintf('\n');

%% =========================================================================
%  PERFORMANCE METRICS
%% =========================================================================
settle_thresh = 0.05;

fprintf('=== Full Comparison — Performance Metrics ===\n');
fprintf('%-40s %-12s %-10s %-10s %-12s\n', ...
        'Controller', 'Settle(s)', 'Max|u1|', 'Max|u2|', 'RMS|u|');
fprintf('%s\n', repmat('-', 1, 86));

metrics = struct();
for m = 1:n_ctrl
    t_m = t_store{m};
    w_m = results{m};
    u_m = u_store{m};

    metrics(m).settle = compute_settling_time(t_m, w_m, settle_thresh);
    metrics(m).max_u1 = max(abs(u_m(1,:)));
    metrics(m).max_u2 = max(abs(u_m(2,:)));
    metrics(m).rms_u  = sqrt(mean(u_m(1,:).^2 + u_m(2,:).^2));

    % Overshoot: max deviation after first 5% of simulation
    % (skips initial condition x0=1 which would always register as overshoot)
    skip_idx     = max(2, round(length(t_m) * 0.05));
    metrics(m).os_x = max(abs(w_m(1, skip_idx:end)));
    metrics(m).os_y = max(abs(w_m(2, skip_idx:end)));

    fprintf('%-40s %-12.3f %-10.2f %-10.2f %-12.3f\n', ...
            ctrl_labels{m}, metrics(m).settle, ...
            metrics(m).max_u1, metrics(m).max_u2, metrics(m).rms_u);
end

% Save metrics CSV
if ~exist('../results','dir'); mkdir('../results'); end
fid = fopen('../results/metrics_comparison.csv', 'w');
fprintf(fid, 'Controller,Settling_time_s,Max_u1,Max_u2,RMS_u,Overshoot_x,Overshoot_y\n');
for m = 1:n_ctrl
    fprintf(fid, '%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n', ...
            ctrl_labels{m}, metrics(m).settle, metrics(m).max_u1, ...
            metrics(m).max_u2, metrics(m).rms_u, ...
            metrics(m).os_x, metrics(m).os_y);
end
fclose(fid);
fprintf('\nMetrics saved to results/metrics_comparison.csv\n');

%% =========================================================================
%  FIGURES
%% =========================================================================
if ~exist('../docs/figures','dir'); mkdir('../docs/figures'); end

state_labels = {'x (m)', 'y (m)', 'dx/dt (m/s)', 'dy/dt (m/s)'};
state_titles = {'Position x', 'Position y', 'Velocity dx/dt', 'Velocity dy/dt'};

%% Figure 1: State Comparison
figure('Name','Full Comparison — States','NumberTitle','off', ...
       'Position',[30 30 1100 750]);
for s = 1:4
    subplot(2,2,s)
    for m = 1:n_ctrl
        plot(t_store{m}, results{m}(s,:), lstyles{m}, ...
             'Color', colors{m}, 'LineWidth', 2, ...
             'DisplayName', ctrl_labels{m});
        hold on;
    end
    yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel('Time (s)'); ylabel(state_labels{s});
    title(state_titles{s});
    legend('Location', 'best', 'FontSize', 7);
    grid on;
end
sgtitle('Full System — State Response Comparison (All Controllers)', ...
        'FontSize', 13, 'FontWeight', 'bold');
saveas(gcf, '../docs/figures/08a_full_state_comparison.png');
fprintf('Saved: 08a_full_state_comparison.png\n');

%% Figure 2: Control Effort
figure('Name','Full Comparison — Control','NumberTitle','off', ...
       'Position',[30 30 1000 600]);
for c = 1:2
    subplot(2,1,c)
    for m = 1:n_ctrl
        plot(t_store{m}, u_store{m}(c,:), lstyles{m}, ...
             'Color', colors{m}, 'LineWidth', 2, ...
             'DisplayName', ctrl_labels{m});
        hold on;
    end
    yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel('Time (s)'); ylabel(sprintf('u_%d', c));
    title(sprintf('Control Input u_%d', c));
    legend('Location', 'best', 'FontSize', 7);
    grid on;
end
sgtitle('Full System — Control Effort Comparison (All Controllers)', ...
        'FontSize', 13, 'FontWeight', 'bold');
saveas(gcf, '../docs/figures/08b_full_control_comparison.png');
fprintf('Saved: 08b_full_control_comparison.png\n');

%% Figure 3: State Norm Convergence
figure('Name','Convergence','NumberTitle','off','Position',[30 30 900 500]);
for m = 1:n_ctrl
    sn = vecnorm(results{m}(1:4,:)) + 1e-12;
    semilogy(t_store{m}, sn, lstyles{m}, 'Color', colors{m}, ...
             'LineWidth', 2.5, 'DisplayName', ctrl_labels{m});
    hold on;
end
yline(settle_thresh, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
text(T*0.97, settle_thresh*1.4, '5% threshold', 'FontSize', 8, ...
     'HorizontalAlignment', 'right', 'Color', [0.3 0.3 0.3]);
xlabel('Time (s)', 'FontSize', 11);
ylabel('||w|| (log scale)', 'FontSize', 11);
title('State Norm Convergence — All Controllers', 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 8);
grid on;
saveas(gcf, '../docs/figures/08c_convergence_comparison.png');
fprintf('Saved: 08c_convergence_comparison.png\n');

%% Figure 4: Pareto Front with all controllers marked
figure('Name','Pareto with all methods','NumberTitle','off', ...
       'Position',[30 30 950 680]);

plot(rms_u_valid, settle_valid, 'b-o', 'LineWidth', 2, ...
     'MarkerSize', 4, 'MarkerFaceColor', 'b', ...
     'DisplayName', 'LQR Pareto Front');
hold on;

marker_shapes = {'s', 's', '^', 'p', 'h'};
for m = 1:n_ctrl
    scatter(metrics(m).rms_u, metrics(m).settle, 200, ...
            colors{m}, marker_shapes{min(m,5)}, 'LineWidth', 2.5, ...
            'DisplayName', ctrl_labels{m});
end

all_u = [rms_u_valid, arrayfun(@(m) metrics(m).rms_u, 1:n_ctrl)];
all_s = [settle_valid, arrayfun(@(m) metrics(m).settle, 1:n_ctrl)];
all_u = all_u(isfinite(all_u) & all_u < 500);
all_s = all_s(isfinite(all_s) & all_s < 90);
x_pad = (max(all_u)-min(all_u))*0.12;
y_pad = (max(all_s)-min(all_s))*0.12;
if x_pad > 0 && y_pad > 0
    xlim([max(0,min(all_u)-x_pad), max(all_u)+x_pad]);
    ylim([max(0,min(all_s)-y_pad), max(all_s)+y_pad]);
end

xlabel('RMS Control Effort  \surd(mean(u_1^2 + u_2^2))', 'FontSize', 12);
ylabel('Settling Time (s)', 'FontSize', 12);
title({'LQR Pareto Front with All Controller Results', ...
       'Points on curve = Pareto optimal;  above/right = suboptimal'}, ...
      'FontSize', 11);
legend('Location', 'best', 'FontSize', 8);
grid on;
saveas(gcf, '../docs/figures/08d_pareto_all_controllers.png');
fprintf('Saved: 08d_pareto_all_controllers.png\n');

fprintf('\nAll figures saved to docs/figures/\n');
fprintf('=== Simulation Complete ===\n');