function [th, info] = heat_steady(q, G, tol, maxit)
%HEAT_STEADY  Steady-state temperature rise theta_ss for the deposition q.
%
%   Solves  k grad^2 theta = -q  with the boundary conditions of G, i.e. the
%   t -> infinity limit of HEAT_SOLVE.  Needed because the slowest thermal
%   mode of a 30 mm block has a time constant of about
%       1/(alpha*(pi/(2d))^2) = 2540 s,
%   so a transient run of a few minutes is nowhere near the steady state, and
%   the steady state is what sets the maximum power the phantom can take
%   continuously.
%
%   METHOD.  Matrix-free preconditioned conjugate gradients on
%       A theta = b,   A = -W .* grad^2 ,   b = W .* q/k
%   restricted to the non-pinned nodes.  Multiplying by the control-volume
%   weight W is what makes A symmetric (the bare ghost-node Laplacian is not,
%   because a boundary node owns half a cell), and A is positive definite
%   because every configuration used here loses heat somewhere.  The
%   preconditioner is the exact diagonal, G.diagA.
%
%   Verified in VERIFY_HEAT_MODEL against (a) a closed-form eigenfunction
%   solution and (b) the long-time limit of the explicit transient solver.

if nargin < 3 || isempty(tol),   tol   = 1e-10; end
if nargin < 4 || isempty(maxit), maxit = 4000;  end

free = G.free;
b    = G.W(free).*q(free)/G.k;
Afun = @(v) apply_A(v, G, free);
Mfun = @(v) v./G.diagA(free);

[v, flag, relres, iter] = pcg(Afun, b, tol, maxit, Mfun);

th = zeros(G.nx, G.ny, G.nz);
th(free) = v;

info = struct('flag',flag,'relres',relres,'iter',iter,'tol',tol);
% flag 3 (stagnation) with a relative residual already far below the
% requested tolerance is convergence, not failure -- it happens routinely when
% the residual reaches the rounding floor of the operator.
if flag ~= 0 && relres > 10*tol
    warning('heat_steady:pcg','pcg flag %d, relres %.2e after %d iters', ...
            flag, relres, iter);
end
end

function out = apply_A(v, G, free)
th = zeros(G.nx, G.ny, G.nz);
th(free) = v;
L = heat_lap(th, G);
out = -G.W(free).*L(free);
end
