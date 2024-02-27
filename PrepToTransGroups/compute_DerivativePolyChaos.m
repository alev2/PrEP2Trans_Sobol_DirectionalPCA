function [deriv_modal_coeffs, deriv_K] = compute_DerivativePolyChaos(modal_coeffs,K,ind_dvar)

%This takes as an input a polynomial chaos expansion (in the form
%modal_coeffs) as well as its corresponding multi-index (K) and computes
%the derivative with respect to the variable indicated by the variable
%index ind_dvar. This variable should, in principal, be the i-th column of K.


%First, we look in the multi-index to find where this variable is non-zero,
%as well as storing the corresp. exponents
    [nzInds,~,nzPwrs] = find(K(:,ind_dvar)>0);

%get the derivative modal coefficients    
    deriv_modal_coeffs=zeros(size(modal_coeffs));
    deriv_modal_coeffs(nzInds)=modal_coeffs(nzInds);%.*K(nzInds,ind_dvar);

    deriv_K=zeros(size(K));
    deriv_K(nzInds,:)=K(nzInds,:);




end

