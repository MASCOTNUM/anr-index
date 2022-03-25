% This script calls the scripts F_greedy.m and F_lazy_greedy.m, which build
% designs that greedily optimize the criterion I_B,q in the paper
% [A. Nogales Gómez, L. Pronzato and M.-J. Rendas. "Incremental space-filling design 
% based on coverings and spacings: improving upon low discrepancy sequences", 
% J. of Statistical Theory and Practice, 2021 (hal-02987983, arXiv:2106.05833)]

% In this example the domain is the unit hypercube, of dimension d.
% By chosing sets X_Q and X_cand defined over other domains
% the scripts can be called in the same manner as here.
%-----
% uses scripts Covering_Radius.m and NumLabelPlot.m for graphic output
%-----

% Copyright 2021 CNRS [J. Rendas]
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

clear variables
close all

LAZY      =1;  % set LAZY = 0 to call standard greedy
%
d         = 2;                       % dimension of covered space
n         = 10*d;                     % maximum design size
V_d       = pi^(d/2)/(gamma(d/2+1));  % volume of the hypercube
% define the two sets as subsets of the Sobol sequence
Q         = 8*2048;       % size of X_Q
C         = floor(Q/2); % size of X_Cand, can be set to a different value
pS        = sobolset(d);
pSS       = scramble(pS,'MatousekAffineOwen');
k_frac    = 2; %
Complement=((fullfact(k_frac*ones(1,d))-1)/(k_frac-1)); % complement the sets with the corners of the hypercube
X_Q       = [pS(1:Q,:); Complement];  % include the vertices of the hypercube in X_Q
X_cand    = pS(1:C,:);
% define the parameters of the criterion
I_param.q = 10;  % should be in [5, 25]
if (1)  % to use the entire range of distance values use (0)  to set limits
    I_param.b = 0;
    I_param.B = sqrt(d)/2;
else % alternatively, set stricter limits
    nmin      = 5*d;  nmax = 20*d;
    I_param.b = (nmax*V_d)^(-1/d);
    I_param.B = sqrt(d)/2*(floor(nmin^(1/d)))^(-1);
end
% I_param.q = 5; Q= 2048; X_Q = pS(1:Q,:); X_cand = X_Q; % uncomment this line to replicate Figure 2 in the paper

%%

if LAZY
    % call lazy greedy implementation
    [i_a, design, cr]       = Lazy_greedy_cdf(X_cand, X_Q, n, I_param);
else
    % call standard greedy
    [i_a, design, cr]      = Greedy_cdf(X_cand, X_Q, n, I_param);
end

% on output
%   - i_a is the sequence of values of I_A(X_k), from k \in [1, n]
%   - design is the orddered sequence of indexes of the points of X_cand
%   retained as design points: X_k = design(1:k);
%   - cr is the sequence of values of the empirical values of the covering
%   radius: cr(k) = CR(X_k)

% compute the covering radius over a denser grid with 2^18 points
X_cr      = net(pSS,2^18)';
dense_cr  = Covering_Radius( X_cand(design,:)', 1:n, 0, X_cr, 2);  % lazy
% plot design only if d=2
if d==2
    h=figure(1);  % allows to choose the window where NumLabelPlot draws
    set(h, 'Position',[23 263 560 420])
    NumLabelPlot(X_Q, X_cand(design(1:n),1),X_cand(design(1:n),2), h, '*r')
    viscircles(X_cand(design(1:n),1:2), cr(end)*ones(n,1));
    plot([0 1 1 0 0],[0 0 1 1 0],'r', 'Linewidth',2)
    axis 'equal'
    axis tight; axis([0 1 0 1]);
    title(['final design $$(d=' num2str(d) ', q=' num2str(I_param.q) ', b= ' ...
        num2str(I_param.b) ', B= ' num2str(I_param.B)  '$$'],'Interpreter','latex')
    
    
end
% plot the trajectory of CR values

h=figure(2);set(h, 'Position',[592 264 560 420])
hold on
NORM = (V_d*(1:n)).^(-1/d);
plot(dense_cr./NORM,'b')
leg_str =['C= ' num2str(C) ', Q= ' num2str(Q) 'sc_sc'];
title(['cover radius $$(d=' num2str(d)  ', q= ' num2str(I_param.q) ', b= ' num2str(I_param.b) ', B= ' num2str(I_param.B) ')$$'],...
    'Interpreter','latex')
legend(leg_str, 'Location','NorthEast')
xlabel('iteration')
ylabel('normalized cover radius')


