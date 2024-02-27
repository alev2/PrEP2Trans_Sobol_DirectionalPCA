function dL = lege_eval_deriv(x,k,a,b)

% L = lege_eval_deriv(x,k,a,b) 
%
% returns the values of the derivative k-th Legendre polynomial orthoNORMAL in a,b w.r.t to rho=1/(b-a) in the points x
% ( x can be a matrix as well)
%
% N.B. the polynomials start from k=0: L_0(x) = 1, L_1(x) = (2*x-a-b)/(b-a)


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


% this function expresses L as a function of the standard Legendre polynomial (i.e. polynomials orthogonal in -1,1 w.r.t to rho=1 !! ),
% which are recursively calculated through the function standard_lege_eval, coded below in this .m file

% first compute the transformation of x in (a,b) to t in (-1,1)

t = ( 2*x - a - b ) / (b - a ) ;

% calculate the standard Legendre polynomials in t

dL = (2/(b-a))*standard_lege_eval_deriv(t,k,a,b);

% modify L to take into account normalizations

st_lege_norm = sqrt ( 2 / (2*k+1) );

% moreover, add an additional sqrt(2) to take into account a general interval (a,b) , not (-1,1)

dL = sqrt(2) * dL / st_lege_norm;









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function dL = standard_lege_eval_deriv(x,k,a,b)

% L = standard_lege_eval_deriv(x,k) 
%
% returns the values of the first derivative of the k-th standard Legendre polynomial (i.e. polynomials orthogonal and not orthonormal 
% in -1,1 w.r.t to rho=1 !! ) in the points x ( x can be a vector as well)
%
% N.B. the polynomials start from k=0: L_0(x) = 1, L_1(x) = x
%
% This is based on the recurrance relation: 
%
%   d/dx[L(N)] = N*L(N-1) + x*d/dx[L(N-1)]
%
% base steps

% read this as L(k-2)
L_2=ones(size(x));

% and this as L(k-1)
L_1=x;
dL_1=L_2;

if k==0
      L=L_2;
      dL=(b-a)/2;
      return
elseif k==1
      L=L_1;
      dL=dL_1;
      return
else
      % recursive step
      for ric=2:k
            %Compute L(N-1), which we will actually not use until the next k
            L = (2*ric - 1) / ric * x .* L_1 - (ric-1) / ric * L_2;
            %get d/dx[L(N)] = N*L(N-1)+d/dx(L(N-1))
            dL=ric*L_1+x.*dL_1;
            %Save these so we can calculate L(N-1) as needed
            L_2=L_1;
            L_1=L;            
            dL_1=dL;
      end

      return
end