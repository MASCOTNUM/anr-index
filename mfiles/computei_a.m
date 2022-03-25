function I_all = ComputeI_A(R_subset,I_param)
% computes criterion I_A (with parameters I_param) for distances R_subset
% all distances are already raised to q+1)
%-----

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

n  = size(R_subset, 1); % number of candidate design points
Q  = size(R_subset, 2); % grid size

B      = I_param.B;
b      = I_param.b;

% the lines below implement expression (14) in the paper
I_all  = (B- max(b,min(R_subset, B))* ones(Q,1)/Q)/(I_param.q+1);
end

