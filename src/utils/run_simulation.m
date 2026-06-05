function [w, w_hat, a, u_log, t] = run_simulation(A, B, C, E, G1, G2, ...
                                                    Kx, Ki, L, ...
                                                    w0, w_hat0, ...
                                                    dt, T, w1)
% RUN_SIMULATION  Forward Euler simulation of the full closed-loop system:
%                 plant + servo-compensator + Luenberger observer.
%
% Control law:        u = -Kx*w_hat - Ki*a
% Observer:           w_hat_dot = A*w_hat + B*u + L*(z - C*w_hat)
% Servo-compensator:  a_dot = G1*a + G2*z
% Real plant:         w_dot = A*w + B*u + E*d
%
% Inputs:
%   A, B, C, E     - plant matrices
%   G1, G2         - servo-compensator matrices
%   Kx             - 2x4 plant state feedback gain
%   Ki             - 2x6 servo-compensator feedback gain
%   L              - 4x2 observer gain
%   w0             - 4x1 true initial state
%   w_hat0         - 4x1 observer initial estimate
%   dt             - time step (s)
%   T              - simulation duration (s)
%   w1             - disturbance frequency (rad/s)
%
% Outputs:
%   w      - 4xN true state trajectory
%   w_hat  - 4xN estimated state trajectory
%   a      - 6xN servo-compensator state trajectory
%   u_log  - 2xN control input history
%   t      - 1xN time vector

    t = 0:dt:T;
    N = length(t);

    w      = zeros(4, N);
    w_hat  = zeros(4, N);
    a      = zeros(6, N);
    u_log  = zeros(2, N);

    w(:,1)     = w0;
    w_hat(:,1) = w_hat0;
    a(:,1)     = zeros(6,1);

    for k = 1:N-1
        % Disturbances
        d1 = 0.5 * sin(w1 * t(k));
        d2 = cos(w1 * t(k) + pi/4);
        d  = [d1; d2];

        % Measured output (position only)
        z     = C * w(:,k);
        z_hat = C * w_hat(:,k);

        % Control using estimated states + servo state
        u          = -Kx * w_hat(:,k) - Ki * a(:,k);
        u_log(:,k) = u;

        % Real plant (with disturbance)
        w_dot    = A * w(:,k) + B * u + E * d;
        w(:,k+1) = w(:,k) + dt * w_dot;

        % Servo-compensator (driven by measured output)
        a_dot    = G1 * a(:,k) + G2 * z;
        a(:,k+1) = a(:,k) + dt * a_dot;

        % Observer (disturbance-free model + correction)
        w_hat_dot    = A * w_hat(:,k) + B * u + L * (z - z_hat);
        w_hat(:,k+1) = w_hat(:,k) + dt * w_hat_dot;
    end

    % Fill last control input
    u_log(:,N) = -Kx * w_hat(:,N) - Ki * a(:,N);
end
