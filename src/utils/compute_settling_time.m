function st = compute_settling_time(t, w, threshold)
% COMPUTE_SETTLING_TIME  Computes settling time as the first time the state
%                        norm ||w|| drops below a threshold and stays there.
%
% Inputs:
%   t         - time vector [1 x N]
%   w         - state matrix [n x N]  
%   threshold - settling threshold (default: 0.01)
%
% Output:
%   st - settling time in seconds (returns t(end) if never settled)
%


    if nargin < 3
        threshold = 0.01;
    end

    state_norm = vecnorm(w);   
    idx        = find(state_norm < threshold, 1, 'first');

    if isempty(idx)
        st = t(end);
        warning('System did not settle within simulation time T=%.1f s', t(end));
    else
        st = t(idx);
    end
end