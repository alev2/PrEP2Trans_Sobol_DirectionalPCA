dominio=domain;
baseDomain=domain(1,:);
%baseDomain(end-1:end)=domain(2,end-1:end);
%baseDomain(end-5:end-4)=domain(2,end-5:end-4);
%dominio(:,end-1:end)=flipud(dominio(:,end-1:end));
%dominio(:,end-5:end-4)=flipud(dominio(:,end-5:end-4));

%evalDomain=interventionDirections(:,2)';
%evalDomain=directions(:,4)';
evalDomain=[0. 0. 0. 0. 1. 1. 1. 0. 0. 0.];

num_Yrs=length(yearInds);
consideredYears=yearRange(yearInds);

start_YrSum=2024;
end_YrSum=2031;

yr1=find(consideredYears==start_YrSum);
yr2=find(consideredYears==end_YrSum);


%indices: 
% 1. PrEP HETM
% 2. PrEP HETF
% 3. PrEP MSM
% 4. PrEP PWID
% 5. ART HET
% 6. ART MSM
% 7. ART PWID
% 8. Test HET
% 9. Test MSM
% 10. Test PWID
idx_dx=5;
%relevantInterval=linspace(dominio(1,idx_dx),dominio(2,idx_dx),101);
relevantInterval=evalDomain(idx_dx)*ones(2,1);% linspace(dominio(1,idx_dx),dominio(2,idx_dx),101);



polyChaosEval=0;
XX=Sr.knots;
output_Of_Interest=[];
output_Of_Interest_Deriv=[];
output_Of_Interest_DerivAnalytic=[];
output_Of_Interest_KnotPts=[];
       
[modal_coeffs1,K1] = convert_to_modal(S,Sr,sum(values_g1(yr1:yr2,:)),domain,'legendre');    
[modal_coeffs2,K2] = convert_to_modal(S,Sr,sum(values_g2(yr1:yr2,:)),domain,'legendre');    
[modal_coeffs3,K3] = convert_to_modal(S,Sr,sum(values_g3(yr1:yr2,:)),domain,'legendre');    
[modal_coeffs4,K4] = convert_to_modal(S,Sr,sum(values_g4(yr1:yr2,:)),domain,'legendre');    
[modal_coeffs5,K5] = convert_to_modal(S,Sr,sum(values_g5(yr1:yr2,:)),domain,'legendre');    
[modal_coeffs6,K6] = convert_to_modal(S,Sr,sum(values_g6(yr1:yr2,:)),domain,'legendre');    
[modal_coeffs7,K7] = convert_to_modal(S,Sr,sum(values_g7(yr1:yr2,:)),domain,'legendre');    
[modal_coeffs8,K8] = convert_to_modal(S,Sr,sum(values_g8(yr1:yr2,:)),domain,'legendre');    
[modal_coeffs9,K9] = convert_to_modal(S,Sr,sum(values_g9(yr1:yr2,:)),domain,'legendre');    
[modal_coeffs10,K10] = convert_to_modal(S,Sr,sum(values_g10(yr1:yr2,:)),domain,'legendre');    
[modal_coeffs11,K11] = convert_to_modal(S,Sr,sum(values_g11(yr1:yr2,:)),domain,'legendre');    
[modal_coeffs12,K12] = convert_to_modal(S,Sr,sum(values_g12(yr1:yr2,:)),domain,'legendre');    
[modal_coeffs13,K13] = convert_to_modal(S,Sr,sum(values_g13(yr1:yr2,:)),domain,'legendre');    
[modal_coeffs14,K14] = convert_to_modal(S,Sr,sum(values_g14(yr1:yr2,:)),domain,'legendre');    

[deriv_modal_coeffs1,deriv_K1]=compute_DerivativePolyChaos(modal_coeffs1,K1,idx_dx);
[deriv_modal_coeffs2,deriv_K2]=compute_DerivativePolyChaos(modal_coeffs2,K2,idx_dx);
[deriv_modal_coeffs3,deriv_K3]=compute_DerivativePolyChaos(modal_coeffs3,K3,idx_dx);
[deriv_modal_coeffs4,deriv_K4]=compute_DerivativePolyChaos(modal_coeffs4,K4,idx_dx);
[deriv_modal_coeffs5,deriv_K5]=compute_DerivativePolyChaos(modal_coeffs5,K5,idx_dx);
[deriv_modal_coeffs6,deriv_K6]=compute_DerivativePolyChaos(modal_coeffs6,K6,idx_dx);
[deriv_modal_coeffs7,deriv_K7]=compute_DerivativePolyChaos(modal_coeffs7,K7,idx_dx);
[deriv_modal_coeffs8,deriv_K8]=compute_DerivativePolyChaos(modal_coeffs8,K8,idx_dx);
[deriv_modal_coeffs9,deriv_K9]=compute_DerivativePolyChaos(modal_coeffs9,K9,idx_dx);
[deriv_modal_coeffs10,deriv_K10]=compute_DerivativePolyChaos(modal_coeffs10,K10,idx_dx);
[deriv_modal_coeffs11,deriv_K11]=compute_DerivativePolyChaos(modal_coeffs11,K11,idx_dx);
[deriv_modal_coeffs12,deriv_K12]=compute_DerivativePolyChaos(modal_coeffs12,K12,idx_dx);
[deriv_modal_coeffs13,deriv_K13]=compute_DerivativePolyChaos(modal_coeffs13,K13,idx_dx);
[deriv_modal_coeffs14,deriv_K14]=compute_DerivativePolyChaos(modal_coeffs14,K14,idx_dx);



%outputs (for mixing and prep, expanded output space)
%1. MSM incidence
%2. HETF incidence
%3. HETM incidence
%4. PWID incidence
%5. All incidence
%6. Spending
%7. Num PrEP MSM
%8. Num PrEP HETF
%9. Num PrEP HETM
%10. Num PrEP PWID
%11. Coverage MSM
%12. Coverage HETF
%13. Coverage HETM
%14. Coverage PWID

for j=1:length(relevantInterval)

     
    evalDomain(idx_dx)=relevantInterval(j);
    %evalDomain(idx_dx2)=relevantInterval(j);
    %evalDomain(idx_dx3)=relevantInterval(j);
     
    polyChaosEval=zeros(1,14);
    polyChaosEvalDeriv=zeros(1,14);
     
     for i=1:length(modal_coeffs1)          

       polyChaosEval(1)=polyChaosEval(1)+modal_coeffs1(i)*lege_eval_multidim(evalDomain',K1(i,:),domain(1,:),domain(2,:));
       polyChaosEval(2)=polyChaosEval(2)+modal_coeffs2(i)*lege_eval_multidim(evalDomain',K2(i,:),domain(1,:),domain(2,:));
       polyChaosEval(3)=polyChaosEval(3)+modal_coeffs3(i)*lege_eval_multidim(evalDomain',K3(i,:),domain(1,:),domain(2,:));
       polyChaosEval(4)=polyChaosEval(4)+modal_coeffs4(i)*lege_eval_multidim(evalDomain',K4(i,:),domain(1,:),domain(2,:));
       polyChaosEval(5)=polyChaosEval(5)+modal_coeffs5(i)*lege_eval_multidim(evalDomain',K5(i,:),domain(1,:),domain(2,:));
       polyChaosEval(6)=polyChaosEval(6)+modal_coeffs6(i)*lege_eval_multidim(evalDomain',K6(i,:),domain(1,:),domain(2,:));
       polyChaosEval(7)=polyChaosEval(7)+modal_coeffs7(i)*lege_eval_multidim(evalDomain',K7(i,:),domain(1,:),domain(2,:));
       polyChaosEval(8)=polyChaosEval(8)+modal_coeffs8(i)*lege_eval_multidim(evalDomain',K8(i,:),domain(1,:),domain(2,:));
       polyChaosEval(9)=polyChaosEval(9)+modal_coeffs9(i)*lege_eval_multidim(evalDomain',K9(i,:),domain(1,:),domain(2,:));
       polyChaosEval(10)=polyChaosEval(10)+modal_coeffs10(i)*lege_eval_multidim(evalDomain',K10(i,:),domain(1,:),domain(2,:));
       polyChaosEval(11)=polyChaosEval(11)+modal_coeffs11(i)*lege_eval_multidim(evalDomain',K11(i,:),domain(1,:),domain(2,:));
       polyChaosEval(12)=polyChaosEval(12)+modal_coeffs12(i)*lege_eval_multidim(evalDomain',K12(i,:),domain(1,:),domain(2,:));
       polyChaosEval(13)=polyChaosEval(13)+modal_coeffs13(i)*lege_eval_multidim(evalDomain',K13(i,:),domain(1,:),domain(2,:));
      polyChaosEval(14)=polyChaosEval(14)+modal_coeffs14(i)*lege_eval_multidim(evalDomain',K14(i,:),domain(1,:),domain(2,:));

       polyChaosEvalDeriv(1)=polyChaosEvalDeriv(1)+deriv_modal_coeffs1(i)*lege_eval_multidim_partialDer(evalDomain',deriv_K1(i,:),domain(1,:),domain(2,:),idx_dx);
       polyChaosEvalDeriv(2)=polyChaosEvalDeriv(2)+deriv_modal_coeffs2(i)*lege_eval_multidim_partialDer(evalDomain',deriv_K2(i,:),domain(1,:),domain(2,:),idx_dx);
       polyChaosEvalDeriv(3)=polyChaosEvalDeriv(3)+deriv_modal_coeffs3(i)*lege_eval_multidim_partialDer(evalDomain',deriv_K3(i,:),domain(1,:),domain(2,:),idx_dx);
       polyChaosEvalDeriv(4)=polyChaosEvalDeriv(4)+deriv_modal_coeffs4(i)*lege_eval_multidim_partialDer(evalDomain',deriv_K4(i,:),domain(1,:),domain(2,:),idx_dx);
       polyChaosEvalDeriv(5)=polyChaosEvalDeriv(5)+deriv_modal_coeffs5(i)*lege_eval_multidim_partialDer(evalDomain',deriv_K5(i,:),domain(1,:),domain(2,:),idx_dx);
       polyChaosEvalDeriv(6)=polyChaosEvalDeriv(6)+deriv_modal_coeffs6(i)*lege_eval_multidim_partialDer(evalDomain',deriv_K6(i,:),domain(1,:),domain(2,:),idx_dx);
       polyChaosEvalDeriv(7)=polyChaosEvalDeriv(7)+deriv_modal_coeffs7(i)*lege_eval_multidim_partialDer(evalDomain',deriv_K7(i,:),domain(1,:),domain(2,:),idx_dx);
       polyChaosEvalDeriv(8)=polyChaosEvalDeriv(8)+deriv_modal_coeffs8(i)*lege_eval_multidim_partialDer(evalDomain',deriv_K8(i,:),domain(1,:),domain(2,:),idx_dx);
       polyChaosEvalDeriv(9)=polyChaosEvalDeriv(9)+deriv_modal_coeffs9(i)*lege_eval_multidim_partialDer(evalDomain',deriv_K9(i,:),domain(1,:),domain(2,:),idx_dx);
       polyChaosEvalDeriv(10)=polyChaosEvalDeriv(10)+deriv_modal_coeffs10(i)*lege_eval_multidim_partialDer(evalDomain',deriv_K10(i,:),domain(1,:),domain(2,:),idx_dx);
       polyChaosEvalDeriv(11)=polyChaosEvalDeriv(11)+deriv_modal_coeffs11(i)*lege_eval_multidim_partialDer(evalDomain',deriv_K11(i,:),domain(1,:),domain(2,:),idx_dx);
       polyChaosEvalDeriv(12)=polyChaosEvalDeriv(12)+deriv_modal_coeffs12(i)*lege_eval_multidim_partialDer(evalDomain',deriv_K12(i,:),domain(1,:),domain(2,:),idx_dx);
       polyChaosEvalDeriv(13)=polyChaosEvalDeriv(13)+deriv_modal_coeffs13(i)*lege_eval_multidim_partialDer(evalDomain',deriv_K13(i,:),domain(1,:),domain(2,:),idx_dx);
       polyChaosEvalDeriv(14)=polyChaosEvalDeriv(14)+deriv_modal_coeffs14(i)*lege_eval_multidim_partialDer(evalDomain',deriv_K14(i,:),domain(1,:),domain(2,:),idx_dx);

     end

    output_Of_Interest=[output_Of_Interest;polyChaosEval];
    output_Of_Interest_DerivAnalytic=[output_Of_Interest_DerivAnalytic;polyChaosEvalDeriv];

end



% for j=1:length(Sr.knots)
% 
%      
%     evalDomain(idx_dx)=Sr.knots(j);
%     %evalDomain(idx_dx2)=relevantInterval2(j);
%     %evalDomain(idx_dx3)=relevantInterval3(j);
%      
%     polyChaosEval=zeros(1,12);
%      
%      for i=1:length(modal_coeffs1)          
%          
%        polyChaosEval(1)=polyChaosEval(1)+modal_coeffs1(i)*lege_eval_multidim(evalDomain',K1(i,:),domain(1,:),domain(2,:));
%        polyChaosEval(2)=polyChaosEval(2)+modal_coeffs2(i)*lege_eval_multidim(evalDomain',K2(i,:),domain(1,:),domain(2,:));
%        polyChaosEval(3)=polyChaosEval(3)+modal_coeffs3(i)*lege_eval_multidim(evalDomain',K3(i,:),domain(1,:),domain(2,:));
%        polyChaosEval(4)=polyChaosEval(4)+modal_coeffs4(i)*lege_eval_multidim(evalDomain',K4(i,:),domain(1,:),domain(2,:));
%        polyChaosEval(5)=polyChaosEval(5)+modal_coeffs5(i)*lege_eval_multidim(evalDomain',K5(i,:),domain(1,:),domain(2,:));
%        polyChaosEval(6)=polyChaosEval(6)+modal_coeffs6(i)*lege_eval_multidim(evalDomain',K6(i,:),domain(1,:),domain(2,:));
%        polyChaosEval(7)=polyChaosEval(7)+modal_coeffs7(i)*lege_eval_multidim(evalDomain',K7(i,:),domain(1,:),domain(2,:));
%        polyChaosEval(8)=polyChaosEval(8)+modal_coeffs8(i)*lege_eval_multidim(evalDomain',K8(i,:),domain(1,:),domain(2,:));
%        polyChaosEval(9)=polyChaosEval(9)+modal_coeffs9(i)*lege_eval_multidim(evalDomain',K9(i,:),domain(1,:),domain(2,:));
%        polyChaosEval(10)=polyChaosEval(10)+modal_coeffs10(i)*lege_eval_multidim(evalDomain',K10(i,:),domain(1,:),domain(2,:));
%        polyChaosEval(11)=polyChaosEval(11)+modal_coeffs11(i)*lege_eval_multidim(evalDomain',K11(i,:),domain(1,:),domain(2,:));
%        polyChaosEval(12)=polyChaosEval(12)+modal_coeffs12(i)*lege_eval_multidim(evalDomain',K12(i,:),domain(1,:),domain(2,:));
% 
% 
%      end
% 
%     output_Of_Interest_KnotPts=[output_Of_Interest_KnotPts;polyChaosEval];
% 
% end

dx=relevantInterval(2)-relevantInterval(1);

output_Of_Interest_Deriv=(output_Of_Interest(3:end,:)-output_Of_Interest(1:end-2,:))/(2*dx);
output_Of_Interest_Deriv2=(output_Of_Interest(3:end,:)-2*output_Of_Interest(2:end-1,:) + output_Of_Interest(1:end-2,:))/((relevantInterval(2)-relevantInterval(1))^2);

%plotPolyChaos