function [Sob_i,Tot_Sob_i,Mean,Var,Vi,VTi] = compute_sobol_indices_from_sparse_grid_with_history(S,Sr,PrevVi,PrevVTi,PrevVar,nodal_values,domain,flags)


% COMPUTE_SOBOL_INDICES_FROM_SPARSE_GRID_WITH_HISTORY computes the modified Sobol indices (as introduced in Alexanderian et al., 2020)
%    of a scalar-valued function f:R^N -> R  in three steps: 
% 1) converts the sparse grid approximation of that function into its equivalent Polynomial Chaos Expansion (PCE);
%    this operation is performed by calling CONVERT_TO_MODAL
% 2) performs algebraic manipulations of the PCE coefficients.
% 3) the PCE coefficients are then integrated in time, and compared with
% the total variance over the time period
%
% [SOB_I,TOT_SOB_I,MEAN,VAR] = COMPUTE_SOBOL_INDICES_FROM_SPARSE_GRID_WITH_HISTORY(S,SR,PrevVi,PrevVTI, PREVVAR, NODAL_VALUES,DOMAIN,FLAGS) takes the
%       same inputs of CONVERT_TO_MODAL (see step 1 above) and returns the Sobol indices of the function F. 
%
%       In details:
%       1) SOB_I is the principal Sobol index of the i-th variable, y_i; i.e., the fraction of variability of F that
%          can be ascribed to y_i only
%       2) TOT_SOB_I is the total Sobol index of y_i; i.e., the fraction of variability of F that can be ascribed
%          to y_i either alone or combined with any other variable.
%       3) We need the HISTORIES of PrevVi, PrevVTi, and PrevVar so that we
%          can perform the necessary time integration in order to produce
%          the history-dependent sobol indices
%    
%       As by-products, the function returns also:
%       3) MEAN, the expected value of F, computed as the coefficient of the constant polynomial in the PCE
%       4) VAR, the variance of F, computed as sum of squares of PCE
%       coefficients minus MEAN^2 at each time step (a vector)
%       5) Vi, the first-order variation (a vector)
%       6) VTi, the total variation (a vector)
 
%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile, G. Porta
% See LICENSE.txt for license
%----------------------------------------------------


% first, call convert_to_modal
[gpce_coeffs,idx_set] = convert_to_modal(S,Sr,nodal_values,domain,flags);

N = size(Sr.knots,1);
Nb_coeff = length(gpce_coeffs);

% I create two matrices that I use to mark the coefficients that I need to sum to get the Sobol indices.
% Both have one column for each variable and one row for each gPCE coefficient.
%
% In the first one, I mark with 1 the entry (i,j) if the i-th coefficient must be used to compute
% the principal Sobol index of the j-th variable. 
%
% In the second one, I do the same but for the Total Sobol index
Fact_STi = zeros(Nb_coeff,N);
Fact_Si = zeros(Nb_coeff,N);

% loop on coefficients, for each one decide if I should mark it or not. Note that the first PCE coeff 
% gives the mean of the function and does not enter the Sobol index computation, so I can start from 
% the 2nd coefficient

for r=2:Nb_coeff

    % check the multi-idx of this coefficient and find which entries are non-zero
    Ind=find(idx_set(r,:)>0);     
    
    % the current coefficients is then to be used to compute the total sobol index
    % of all those variables
    Fact_STi(r,Ind) = 1;    
    
    % moreover, if there is only one non-zero entry in the multi-idx, than this coefficient
    % must be marked for use in the computation of the principal sobol index
    if length(Ind)==1
        Fact_Si(r,Ind) = 1;
    end
    
end

% compute mean and variance from the PCE coefficients
Mean = gpce_coeffs(1);
gpce_squared = (gpce_coeffs.^2);
TmpVar=sum(gpce_squared)-Mean.^2;

%here we append the variance to the previous variance
Var=[PrevVar TmpVar];

% now for each variable, sum all coefficients that have been marked as contributing to the Total Sobol index, 
% then normalize by variance
VTi=(Fact_STi')*gpce_squared;
%here we append the total variation to the previous total variation
VTi=[PrevVTi VTi];
%here we compute the generalized index
Tot_Sob_i=sum(VTi,2)/sum(Var);

% repeat for the principal Sobol index
Vi=(Fact_Si')*gpce_squared;
Vi=[PrevVi Vi];
Sob_i=sum(Vi,2)/sum(Var);