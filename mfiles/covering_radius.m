function [ CR ] = Covering_Radius( Xn, sequence, d_threshold, Xtest, k_fracfact )
% function [ CR ] = Covering_Radius( Xn, sequence, d_threshold, Xtest, k_fracfact )
% Computes the covering radius of the n point design Xn (a d*n matrix) in [0,1]^d
% sequence is a vector of m integers in {1,...,n}, m<=n
% CR is m row vector of covering radii for designs Xn(:,(1:i)), i in sequence
%   (if sequence = m, CR is just CR(Xn))
% for d<=d_threshold, the exact value is calculated (via Voronoi diagram)
% for d>d_threshold, CR is underestimated by evaluation on the Q grid points in Xtest 
% Xtest is a d*Q matrix, for instance of scrambled Sobol' points:
%   pSmm = sobolset(d); pSmm = scramble(pSmm,'MatousekAffineOwen');
%   Xtest=(net(pS,2^19))';
% To keep the computational cost reasonable, d_threshold should be <= 5 (this is enforced)
% When k_fracfact=k>0, Xtest is completed by a k^d fractional factorial
% design (==> k should be kept small form large d)
%-----

% Copyright 2021 CNRS [L. Pronzato]
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

d_threshold=min(d_threshold,5);
[d,n]=size(Xn);

if  k_fracfact>0
% if isempty(Xtest)==0 && k_fracfact>0
    % complete by a k^d fractional factorial design
    k=k_fracfact;
    Complement=((fullfact(k*ones(1,d))-1)/(k-1))';
    Xtest=[Xtest Complement];
    [~,N]=size(Xtest);
end
LCR=length(sequence);
CR=NaN(1,LCR);
timebar= waitbar(0,'covering radius...'); 
if d>d_threshold
    Dist_all=pdist2(Xn',Xtest');
end    
for i=1:LCR
    if mod(i,floor(LCR/10))<1e-2
        waitbar(i/LCR,timebar);
    end
    iCR=sequence(i);
    if d==1
        CR(i) = Minimax_d1( XX, 0,1 );
    elseif d<=d_threshold
        [ CR(i), ~, ~, ~ ] = voronoi_truncated01( Xn(:,1:iCR), [], [], 0 );
    else
        %Dist=pdist2(Xn(:,1:iCR)',Xtest');
        %CR(i)=max(min(Dist,[],1));
        CR(i)=max(min(Dist_all(1:iCR,:),[],1));
    end
end
close(timebar)
end
