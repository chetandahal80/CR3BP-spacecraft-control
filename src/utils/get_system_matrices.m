function [A, B, C, D, E] = get_system_matrices()
% GET_SYSTEM_MATRICES  Returns the CR3BP state-space matrices for the
%                      linearised spacecraft dynamics around the L2 point.
%
% Outputs:
%   A  - 4x4 system matrix
%   B  - 4x2 control input matrix
%   C  - 2x4 output (measurement) matrix
%   D  - 2x2 direct feedthrough matrix (zero)
%   E  - 4x2 disturbance input matrix
%
% System:  w_dot = A*w + B*u + E*d
%          z     = C*w
%
% State vector:  w = [x, y, x_dot, y_dot]^T
% Control input: u = [u1, u2]^T
% Disturbance:   d = [d1, d2]^T
% Output:        z = [x, y]^T
%
% Physical parameters (non-dimensionalised):
%   U_XX =  7.3809  
%   U_YY = -2.1904 

    Uxx =  7.3809;
    Uyy = -2.1904;

    A = [0,   0,   1,  0;
         0,   0,   0,  1;
         Uxx, 0,   0,  2;
         0,   Uyy, -2, 0];

    B = [0, 0;
         0, 0;
         1, 0;
         0, 1];

    C = [1, 0, 0, 0;
         0, 1, 0, 0];

    D = zeros(2, 2);

    E = [0, 0;
         0, 0;
         1, 0;
         0, 2];
end
