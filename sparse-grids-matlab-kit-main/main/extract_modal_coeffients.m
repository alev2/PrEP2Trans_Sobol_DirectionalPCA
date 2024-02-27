function [relevantCoeffs] = extract_modal_coeffients(modal_coeffs,K,var1_idx,var2_idx)

%This function takes a N-term polynomial chaos expansion in j-variables (modal_coeffs) and
%its corresponding multi-index (K, Nxj matrix), as well as
%at least one variable index (var1_idx), corresponding to a column of K. 
%One may optionally provide a second variable index (var2_idx).


%If given just one variable index, this function returns the PCE
%coefficients of all PCE terms containing this variable. These are
%organized in a five-column where the columns are given as:

%1. PCE_coeff: the coefficient of this term in the PCE expansion
%2. Var1_idx: the corresponding variable index (provided as an input)
%3. Var1_pwr: The exponent of var1 in this term in the PCE expansion
%4. Var2_idx: The index of the possible second variable in this term in the
%             PCE expansion. Note if this is the same as Var1, it's just
%             repeated a second time.
%5. Var2_pwr: The exponent of var2 in this PCE expansion term

%If given two variable indices, the output is basically the same, however,
%rather than finding and organizing all the var1-relevant terms, it only
%looks for those terms in which Var1 and Var2 interact directly.

%This function is quite useful when attempting to analyze the relative
%influences of particular terms.


%Note this only works for PCE expansions allowing for two-variable
%interactions.
if(max(sum(K~=0,2))>2)
    error('Incompatible PCE. At the moment this function only works with PCEs with up to two-variable interactions.');
end

%This overrides a possibly incorrect function call in which we don't have
%interaction terms (only principal ones), but the user specified
%intercation terms. In this case we override just provide variable 1.
if(max(sum(K~=0,2))==1)
    fprintf('Principal-only PCE. Only considering variable 1.');
    var2_idx=var1_idx;
end

%This is for the variable names in the output table
varNames={...
    'PCE Coeff',...
    'Var1_idx',...
    'Var1_Pwr',...
    'Var2_idx',...
    'Var2_Pwr'
    };


%Start with the 1 input variable case

if(nargin==3 || var1_idx==var2_idx)

    %First, we look in the multi-index to find where this variable is non-zero,
    %as well as storing the corresp. exponents.  
     
    [nzInds,~,nzPwrs] = find(K(:,var1_idx));

    %Here we get the interaction terms.
    Kpn=K;
    %We want to zero out the columns of var1 when there is an interaction
    %variable in addition to var1, but not those in which not those in which
    % var1 appears by itself.
    Kpn(:,var1_idx)=Kpn(:,var1_idx).*(Kpn(:,var1_idx)==sum(Kpn,2));
    Kpn=Kpn';
    [nzCompInds,~,nzCompPwrs]=find(Kpn(:,nzInds)); %find(Kpn(Kpn(:,nzInds)>0));

    %build the table and output
   relevantCoeffs=[modal_coeffs(nzInds) var1_idx.*(nzInds./nzInds) nzPwrs nzCompInds nzCompPwrs] ;   
   relevantCoeffs=array2table(relevantCoeffs,'VariableNames',varNames);

else

    %case 2: two variable inputs. This is much easier and self explanatory.
    nzInds = find(K(:,var1_idx).*K(:,var2_idx)~=0);
    relevantCoeffs=[modal_coeffs(nzInds) var1_idx*(nzInds./nzInds) K(nzInds,var1_idx) var2_idx*(nzInds./nzInds) K(nzInds,var2_idx)];
    relevantCoeffs=array2table(relevantCoeffs,'VariableNames',varNames);
end
