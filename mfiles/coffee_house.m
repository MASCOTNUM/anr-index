function X = coffee_house( Xcand, nmax, alpha, norepeat )
% function X = coffee_house( Xcand, nmax, alpha, norepeat   )
% = coffee-house design
% computes a design X (d*nmax) that should have 50% minimax and maximin efficiency 
% for all n <= nmax
% The design is searched in the "set" Xcand (d*Q) of candidate points
% No computation of Q*Q interdistances matrix is involved
% alpha >= 0 defines spacings
%          alpha=0 --> standard coffee house (maximin and minimax 50% efficient)
%          alpha=1 --> standard spacings
% (alpha corresponds to 1/beta in the paper "Incremental space-filling design based 
%   on coverings and spacings: improving upon low discrepancy sequences", by 
%   Nogales-Gómez, Pronzato and Rendas, submitted to the Journal of Statistical Theory 
%   and Practice (2021).           
% The first point X(:,1) is always the closest one in Xcand to ones(d,1)/2 (i.e., the 
%   closest to the center of the cube [0,1]^d 
% norepeat = 1 --> forces no repetitions of design points (which may happen when alpha>0 
%                  if all remaining candidates are on the boundary)
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

[d,Q]=size(Xcand);
IndX=zeros(1,nmax);

%------------------
% 1) Initialization
%------------------
Xall=NaN(d,nmax);
% initialize by "central point" (in case of a design [0,1]^d) 
Center=ones(d,1)/2; 
D_to_center=pdist2(Center',Xcand');
[~,i1]=min(D_to_center); X=Xcand(:,i1); IndX(1)=i1;

Xall(:,1)=X;
D_to_X=pdist2(Xcand',X');  % distances to design points

%--------------------
% 2) Start iterations
%--------------------
timebar= waitbar(0,'Coffee-house...'); 
for n=2:nmax
    waitbar(n/nmax,timebar);
    
    [~,ibest]=max(D_to_X);   
    if alpha>0
        if norepeat==0 
            [~,ibest]=max(min([D_to_X [Xcand' -Xcand'+1]/alpha],[],2)); 
            ibest=ibest(1);
        else
            Ind_free=setdiff(1:Q,IndX);
            Dum=[D_to_X [Xcand' -Xcand'+1]/alpha];
            Dum=Dum(Ind_free,:);
            [~,ibest]=max(min(Dum,[],2)); ibest=ibest(1);
            ibest=Ind_free(ibest);
        end
    end 
    D_to_X=min(D_to_X,pdist2(Xcand',Xcand(:,ibest)')); % update all D_to_X     
    x_last=Xcand(:,ibest); 
    Xall(:,n)=x_last;
    X=Xall(:,1:n);
    IndX(n)=ibest;
end
close(timebar) 
end
        
