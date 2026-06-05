function [G1, G2, Aa, Ba, Ea] = get_servo_matrices(A, B, C, E, w1)
% GET_SERVO_MATRICES  Constructs the servo-compensator and augmented system
%                     matrices based on the Internal Model Principle.
%
% The disturbances and reference signal have Laplace denominators:
%   d1 ~ 1/(s^2 + w1^2)
%   d2 ~ 1/(s^2 + w1^2)
%   r  ~ 1/s  (constant reference = 0)
%
% Least common denominator: alpha(s) = s(s^2 + w1^2)
% Degree 3, so companion matrix beta is 3x3.
%
% Since system has 2 outputs (x, y), we need 2 channels:
%   G1 = blkdiag(beta, beta)  [6x6]
%   G2 = blkdiag(sigma, sigma) [6x2]
%
% Augmented state: wa = [w; a]  (10x1)
% Augmented system: wa_dot = Aa*wa + Ba*u + Ea*d
%
% Inputs:
%   A, B, C, E - system matrices from get_system_matrices()
%   w1         - disturbance frequency (rad/s)
%
% Outputs:
%   G1  - 6x6 servo-compensator dynamics matrix
%   G2  - 6x2 servo-compensator input matrix
%   Aa  - 10x10 augmented system matrix
%   Ba  - 10x2 augmented control input matrix
%   Ea  - 10x2 augmented disturbance matrix

    % Companion matrix of alpha(s) = s(s^2 + w1^2) = s^3 + w1^2*s
    % Coefficients: alpha_1=0, alpha_2=w1^2, alpha_3=0
    beta  = [0,  1,      0;
             0,  0,      1;
             0, -w1^2,   0];

    sigma = [0; 0; 1];

    % Two-channel servo-compensator (one per output)
    G1 = blkdiag(beta,  beta);
    G2 = blkdiag(sigma, sigma);

    % Augmented system matrices
    Aa = [A,       zeros(4,6);
          G2*C,    G1        ];

    Ba = [B;
          zeros(6,2)];

    Ea = [E;
          zeros(6,2)];
end
