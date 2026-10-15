% COMPARE_N500_VS_N5000  Convergence check of Ehrenfest runs with N = 500 and N = 5000.
%
%   Compares the mean-field (Ehrenfest) runs in ../n500 and ../n5000 of the
%   same test case. Over the whole time interval it computes
%
%     max_t |Tr(rho_5000^2) - Tr(rho_500^2)|
%     max_t |<q>_5000 - <q>_500|
%     max_t |<p>_5000 - <p>_500|
%     max_t |Var(q)_5000 - Var(q)_500|,   Var(q) = <q^2> - <q>^2
%
%   and the relative differences  max_t |A_5000 - A_500| / max_t |A_5000|
%   for A = <q>, <p>, Var(q).
%
%   rho(t) = sum_k w_k rho_k(t) is the reduced density matrix and
%   <f>(t) = sum_k w_k f(q_k(t), p_k(t)). q and p are used as stored in the
%   result files (q_val, p_val), i.e. in the rescaled units of
%   ehrenfest_simulation.m (q = rescale1*q, p = rescale2*p).
%   All computations use explicit for loops (no vectorisation).

clear; clc

%% Load both runs
file_5000 = '../n5000/ehrenfest_simulation.mat';
file_500  = '../n500/ehrenfest_simulation.mat';

A = observables(file_5000);
B = observables(file_500);

% both runs must be on the same time grid
if A.n ~= B.n
    error('The two runs do not have the same number of time steps.');
end
for m = 1:A.n
    if abs(A.t(m) - B.t(m)) > 1e-12*abs(A.t(A.n))
        error('The two runs do not have the same time grid.');
    end
end

% maximum differences and maximum scales of the N = 5000 curves

%% Maximum differences and maximum values of the N = 5000 curves
d_pur = 0;  d_q = 0;  d_p = 0;  d_vq = 0;
s_q = 0;  s_p = 0;  s_vq = 0;
for m = 1:A.n
    d_pur = max(d_pur, abs(A.pur(m) - B.pur(m)));
    d_q   = max(d_q,   abs(A.mq(m)  - B.mq(m)));
    d_p   = max(d_p,   abs(A.mp(m)  - B.mp(m)));
    d_vq  = max(d_vq,  abs(A.vq(m)  - B.vq(m)));

    s_q  = max(s_q,  abs(A.mq(m)));
    s_p  = max(s_p,  abs(A.mp(m)));
    s_vq = max(s_vq, abs(A.vq(m)));
end

%% Relative differences and report
r_q  = d_q  / s_q;
r_p  = d_p  / s_p;
r_vq = d_vq / s_vq;

fprintf('N = %d vs N = %d,  T = %.4e,  %d time steps\n\n', A.n_b, B.n_b, A.t(A.n), A.n);
fprintf('                 absolute      relative\n');
fprintf('Tr(rho^2)      %.4e\n', d_pur);
fprintf('<q>            %.4e    %.4e\n', d_q,  r_q);
fprintf('<p>            %.4e    %.4e\n', d_p,  r_p);
fprintf('Var(q)         %.4e    %.4e\n', d_vq, r_vq);

%% ========================================================================
%  Local functions
%  ========================================================================

% -------------------------------------------------------------------------
function O = observables(file)
    % Time series of <q>, <p>, Var(q) and the purity Tr(rho^2) of one run.
    S = load(file, 'q_val', 'p_val', 'rho_val', 'w', 'l', 'dt', 'n_b');

    % ehrenfest_simulation.m saves before "l = l-1": stored steps are 1..l-1
    n   = S.l - 1;
    n_b = S.n_b;

    O.n   = n;
    O.n_b = n_b;
    O.t   = zeros(1,n);
    O.mq  = zeros(1,n);    % <q>(t)
    O.mp  = zeros(1,n);    % <p>(t)
    O.vq  = zeros(1,n);    % Var(q)(t)
    O.pur = zeros(1,n);    % Tr(rho^2)(t)

    for m = 1:n
        O.t(m) = (m-1)*S.dt;

        % weighted averages over the particles
        mq  = 0;
        mq2 = 0;
        mp  = 0;
        R   = zeros(2,2);  % rho(t) = sum_k w_k rho_k(t)
        for k = 1:n_b
            qk = S.q_val(k,m);
            pk = S.p_val(k,m);
            mq  = mq  + S.w(k)*qk;
            mq2 = mq2 + S.w(k)*qk^2;
            mp  = mp  + S.w(k)*pk;
            for i = 1:2
                for j = 1:2
                    R(i,j) = R(i,j) + S.w(k)*S.rho_val(i,j,k,m);
                end
            end
        end

        % purity Tr(rho^2) = sum_i sum_j rho_ij rho_ji
        pur = 0;
        for i = 1:2
            for j = 1:2
                pur = pur + R(i,j)*R(j,i);
            end
        end

        O.mq(m)  = mq;
        O.mp(m)  = mp;
        O.vq(m)  = mq2 - mq^2;
        O.pur(m) = real(pur);
    end
end
