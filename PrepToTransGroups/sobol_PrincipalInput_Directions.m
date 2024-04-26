dominio=domain;
%baseDomain=domain(1,:);
%baseDomain(end-1:end)=domain(2,end-1:end);
%baseDomain(end-5:end-4)=domain(2,end-5:end-4);
%dominio(:,end-1:end)=flipud(dominio(:,end-1:end));
%dominio(:,end-5:end-4)=flipud(dominio(:,end-5:end-4));

%evalDomain= baseDomain;
%evalDomain=[.2 .2 1.5 1.25 .85 .85 .2 .2 1.25 1.25 .15 .15 .7 .7 .1 .1];
%evalDomain=[.1 .1 .1 .1 .1 .1 .1 .1 .1 .1 .9 .9 .1 .1 .9 .9];
evalDomain=[1 1 1 1 1]*.0;
%evalDomain=[0. 0. 0. 0. 0.];
%evalDomain=[.0          .01          1.1       0.9776        0.809        0.845        0.165        0.174          1.1          1.1       0.2234       0.2288       0.5679       0.5157      0.12718       0.1253];

% evalDomain([3 4 5 6 7 8 9 10 13 14] )=evalDomain([3 4 5 6 7 8 9 10 13 14])*1.1;
% evalDomain([11 12 15 16])=evalDomain([11 12 15 16])*.9;
% evalDomain([1 2])=.1;


polyChaosEval=0;
XX=Sr.knots;

numOutputs=4;
       

jacobianMatrix=[];

for i=1:numOutputs 
        
       [modal_coeffs,K]=convert_to_modal(S,Sr,reduced_Outputs(i,:),domain,'legendre'); 
       jacobian_CurRow=[];

       for idx_dx=1:length(evalDomain)
         
            jacobian_CurElem=0;
            [deriv_modal_coeffs,deriv_K]=compute_DerivativePolyChaos(modal_coeffs,K,idx_dx);

            for j=1:length(deriv_modal_coeffs)

                jacobian_CurElem=jacobian_CurElem+deriv_modal_coeffs(j)*lege_eval_multidim_partialDer(evalDomain',deriv_K(j,:),domain(1,:),domain(2,:),idx_dx);       

            end 
       
            jacobian_CurRow=[jacobian_CurRow jacobian_CurElem];
       
       end

       jacobianMatrix=[jacobianMatrix;jacobian_CurRow];

end

[Uu,Ss,Vv]=svd(jacobianMatrix);
Uk=Uu(:,1:numOutputs);
Sk=Ss(1:numOutputs,1:numOutputs);
%Sk(1,1)=1/Sk(1,1);
%Sk(2,2)=1/Sk(2,2);
%Sk(3,3)=1/Sk(3,3);

Vk=Vv(:,1:numOutputs);


principalDirections=Vk*(Sk\Uk')*(eye(numOutputs,numOutputs));

%principalDirections=jacobianMatrix\(eye(numOutputs,numOutputs));

