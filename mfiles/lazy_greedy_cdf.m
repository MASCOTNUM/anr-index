function [i_a, X_n, radius]= Lazy_greedy_cdf(X_cand, X_Q, n, I_param)
% builds a design of size n (X_n) chosen amongst the Q points in
% X_cand, greedily optimisng the covering criterion w/ parameters I_param
% (q, b and B), which is computed in the set of points X_Q. I_param is a
% structure defining the parameters of I_A:
% I_param.q, I_param.b and I_param.B
%
% Output:
% - I_A_seq is the sequence of values of the criterion
% - X_n are the indices of the grid points that makeup the final design
% - radius is the sequence of covering radius of the partial designs
%       X_n(1:k), k \in [1,n]
% Uses script ComputeI_A.m
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
I_param.B      = I_param.B^(I_param.q+1);
I_param.b      = I_param.b^(I_param.q+1);

% extract dimensions
Q      = size(X_Q,1);    % size of the grid
Q_cand = size(X_cand,1); % size of the set of possible design points
d      = size(X_Q,2);    % dimension of input space

if size(X_cand,2) ~= d
    display('incompatible dimensions of grid and design sets in F_lazy_greedy. Returning...')
    return
else
    timebar= waitbar(0,'lazy greedy...'); 
end

% compute distance matrix and its power just once
A_dist = pdist2(X_cand, X_Q).^(I_param.q+1);

% allocate space for speed
i_a     = NaN*ones(1,n);
radius  = NaN*ones(1,n);

% Initialisation : start with the best greedy
[C_i, X_n, ~]  = Select_greedy( A_dist, inf(1,Q), I_param);  
R              = A_dist(X_n, : );  % line vector of distance of grid points to the design
i_a(1)         = C_i;
radius(1)      = max(R)^(1/(I_param.q+1));

delta        = NaN(Q_cand,1); % initialize the upper bound on the criterion increments (forces calculation for all)
% incrementally extend
i=1;   
cand          = setdiff(1:Q_cand, X_n);        
% work with normalized distances
while (i < n) 
    waitbar(i/n,timebar, ['design point ' num2str(i)]);
    i             = i+1;
    [C_i, x_new, delta_, I_param]  = Select_greedy_lazy(C_i, delta(cand), A_dist(cand,:), R, I_param); % select the next point
    delta(cand)   = delta_; % store the updated upper bounds
    
    X_n(i)        = cand(x_new); % augment the design
    cand(x_new)   = [ ];         % delete from the set of candidate points
    
    R             = min(R, A_dist(X_n(end),:));
    radius(i)     = max(R)^(1/(I_param.q+1));
    i_a(i)        = C_i;
end

delete(timebar)
%
end

%

function [I_A, x_new, delta, I_param] = Select_greedy_lazy(I_old, delta, A_dist, R, I_param)
% lazzy-greedily selects a point (x_new) in set X_cand according to criterion I_A with
% parameters I_param. 
% I_old is the current value of the criterion (in the current integration
% limits)
% delta is the current upper bound on the improvement 
% A_dist is the power of the distance matrix
% R is the current set of design-2-grid_points distances
% I_param gathers the criterion parameters

% uses UpdateR.m and ComputeI_A.m

% the integral is actually computed over a shrinking interval
B            = min(I_param.B, max(R)); OffSet = (I_param.B - B)/(I_param.q+1); I_param.B = B;
I_old        = I_old - OffSet;
% delta        = delta - OffSet;

if all(isnan(delta) )
    % first point
    R_           = min(R, A_dist); 
    I_all        = ComputeI_A(R_ ,I_param);
    [I_A, x_new] = max(I_all);
    delta        = I_all - I_old;
else
    % not the first design point, delta has been calculated
    LazyList    = 1:size(A_dist,1); % consider all candidates at first
    Updated     = [];
    
    while ~isempty(LazyList)
        to_try                = LazyList(1);
        R_                    = min(R, A_dist(to_try,:)); 
        I_new                 = ComputeI_A(R_,I_param);
        delta(to_try)         = I_new - I_old ;            % update upper bound for this element
        Updated(end+1)        = to_try;
        [delta_max, x_new]    = max(delta(Updated),[], 'omitnan'); % find new max of updated elements
        x_new                 = Updated(x_new);
        LazyList              = find(delta > delta_max); % gather remaining candidate for optimality
    end
    
    I_A   = delta_max + I_old ;
end

end

function [I_A, x_new, I_all] = Select_greedy(A_dist, R, I_param)
% greedily selects a point (x_new) according to criterion I_A with
% parameters I_param. R is the current set of design-2-grid_points
% distances

R_           = min(R, A_dist); 
I_all        = ComputeI_A(R_,I_param);
[I_A, x_new] = max(I_all);
end





