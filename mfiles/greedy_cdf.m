function [I_A_seq,X_n, radius]= Greedy_cdf(X_cand, X_Q, n, I_param)
% builds a design of size n (X_n) chosen amongst the Q_cand points in
% X_cand, greedily optimising the covering criterion w/ parameters I_param
% (q, b and B) and computed with parameters in the structure I_param
% Output:
% - I_A_seq is the sequence of values of the criterion
% - X_n are the indices of the points of X_cand that makeup the design
% - radius is the sequence of covering radius of the partial designs 
%   X_k = X_n(1:k)
% uses script ComputeI_A
%----

% Copyright 2021 CNRS [J. Rendas]
% See [A. Nogales Gómez, L. Pronzato and M.-J. Rendas. "Incremental space-filling design 
% based on coverings and spacings: improving upon low discrepancy sequences", 
% J. of Statistical Theory and Practice, 2021 (hal-02987983, arXiv:2106.05833)]
%
% Redistribution and use in source and binary forms, with or without modification, 
% are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, 
% this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice, 
% this list of conditions and the following disclaimer in the documentation and/or 
% other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
% EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
% SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
% TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR 
% BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN 
% ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH 
% DAMAGE.

% the code works directly with the powers of all distances
I_param.B      = I_param.B^(I_param.q+1); B_in = I_param.B;
I_param.b      = I_param.b^(I_param.q+1);
% extract dimensions
Q      = size(X_Q,1);    % size of the grid
C      = size(X_cand,1); % size of the set of possible design points
d      = size(X_Q,2);    % dimension of input space

if size(X_cand,2) ~= d
    display('incompatible dimensions of grid and design sets in F_greedy. Returning...')
    return
end

% compute distance matrix and its power just once
A_dist = pdist2(X_cand, X_Q).^(I_param.q+1);

% allocate space for speed
I_A_seq = NaN*ones(1,n);
radius  = NaN*ones(1,n);

% Initialisation : start with the best greedy
R            = inf(1,Q);
[C_F, X_n, ~]= Select_greedy(A_dist, R, I_param);   % add a point (forward)
R            = A_dist(X_n, : );                     % R is the line vector of distance of grid points to the design
I_A_seq  (1) = C_F;
radius(1)    = max(R);

% incrementally extend
i=1;
while (i < n) 
        I_param.B     = min(I_param.B, radius(i));             % clip for numerical stability
        i             = i+1;
        cand          = setdiff(1:C, X_n);                     % remaining candidate grid points
        [C_F, x_new]  = Select_greedy(A_dist(cand,:), R, I_param);  % select the next point
        X_n(end+1)    = cand(x_new);                                % extend the design
        R             = min(R, A_dist(X_n(end),:));                 % update distances from grid to the design
        radius(i)     = max(R)^(1/(I_param.q+1));                   % empirical coverign radius
        I_A_seq(i)    = C_F + (B_in - I_param.B)/(I_param.q+1);
end
%
end

function [I_A, x_new, I_all] = Select_greedy(A_dist, R, I_param)
% greedily selects a point (x_new) according to criterion I_A with
% parameters I_param. R is the current set of design-2-grid_points
% distances

R            = min(R, A_dist); 
I_all        = ComputeI_A(R,I_param);
[I_A, x_new] = max(I_all);
end




