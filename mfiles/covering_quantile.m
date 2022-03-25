function [ alphaCR ] = Covering_quantile( Xn, sequence, Xtest, k_fracfact, alpha )
% function [ alphaCR ] = Covering_quantile( Xn, sequence, Xtest, k_fracfact, alpha )
% Computes the alpha-quantiles of the distance of random point in in [0,1]^d
% to the designs Xm for m in sequence, a vector of m integers in {1,...,n}, m<=n
% alphaCR is m row vector of alpha-quantiles for designs Xn(:,(1:i)), i in sequence
%   (if sequence = n, CR is just CR(Xn))
% The alpha-quantiles are underestimated by evaluation on the Q grid points in Xtest 
% Xtest is a d*Q matrix, for instance of scrambled Sobol' points:
%   pSmm = sobolset(d); pSmm = scramble(pSmm,'MatousekAffineOwen');
%   Xtest=(net(pS,2^19))';
% When k_fracfact=k>0, Xtest is completed by a k^d fractional factorial design 
%      (==> k should be kept small form large d, k=2 is reasonable)
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

[d,n]=size(Xn);

if isempty(Xtest)==0 && k_fracfact>0
    % complete by a k^d fractional factorial design
    k=k_fracfact;
    Complement=((fullfact(k*ones(1,d))-1)/(k-1))';
    Xtest=[Xtest Complement];
end
[~,N]=size(Xtest);
index_left_alpha_quantile=max(ceil(alpha*N),1);
LalphaCR=length(sequence);
alphaCR=NaN(1,LalphaCR);
timebar= waitbar(0,'covering quantile...'); 
Dist_all=pdist2(Xn',Xtest');
for i=1:LalphaCR
    waitbar(i/LalphaCR,timebar);
    ialphaCR=sequence(i);
    if i==1
        D2Xn=min(Dist_all(1:ialphaCR,:),[],1);
    else
        D2Xn=min(D2Xn,Dist_all(ialphaCR,:));
    end    
    D2Xnsort=sort(D2Xn); 
    alphaCR(i)=D2Xnsort(index_left_alpha_quantile);
end
close(timebar)
end
