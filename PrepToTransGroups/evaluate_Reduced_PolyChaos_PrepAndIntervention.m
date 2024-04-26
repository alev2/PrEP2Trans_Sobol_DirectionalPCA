dominio=domain;
%baseDomain=domain(1,:);
%baseDomain(end-1:end)=domain(2,end-1:end);
%baseDomain(end-5:end-4)=domain(2,end-5:end-4);
%dominio(:,end-1:end)=flipud(dominio(:,end-1:end));
%dominio(:,end-5:end-4)=flipud(dominio(:,end-5:end-4));

rescaleIntervals=0;

% 
% mmm=[          
% 
% 
% 
% 
%          0.245007075227838
%           1.35360878781941
%         0.0377883736704299
%          -3.07375000984079
% 
% ];


evalDomain_Base=ones(1,10)*.5;
evalDomain_Base2=[0 0 0 0 1 1 1 0 0 0];
newDir=Vk*(Sk\Uk')*((UApx'*(stdMat\([150000; 24000; 14000; 12000;200000;425] - rowMeans(:,1)) )));
newDir2=Vk*(Sk\Uk')*((UApx'*(stdMat\([150000; 24000; 14000; 12000;200000;380] - rowMeans(:,1)) )));
newDir4=Vk*(Sk\Uk')*((UApx'*(stdMat\([130000; 22000; 12000; 11000;175000;475] - rowMeans(:,1)) )));
newDir5=Vk*(Sk\Uk')*((UApx'*(stdMat\([130000; 22000; 12000; 11000;175000;425] - rowMeans(:,1)) )));
newDir3=Vk*(Sk\Uk')*((UApx'*(stdMat\([120000; 18000; 10000; 10000;158000;490] - rowMeans(:,1)) )));
newDir6=Vk*(Sk\Uk')*((UApx'*(stdMat\([120000; 18000; 10000; 10000;158000;450] - rowMeans(:,1)) )));


%newDir=Vk*(Sk\Uk')*((UApx'*(stdMat\(AA(:,1) - rowMeans(:,1)) )));

%newDir=Vk*(Sk\Uk')*((UApx'*(stdMat\([400000; 67000; 36000; 31000; 534000] - rowMeans(:,1)) )));
%newDir=Vk*(Sk\Uk')*((UApx'*(stdMat\(mm - rowMeans(:,1)) )));

%newDir3=principalDirections(:,3)*5;  Vk*(Sk\Uk')*[1;0;0;0];

evalDomain=[evalDomain_Base' evalDomain_Base2' evalDomain_Base'+newDir evalDomain_Base'+newDir2 evalDomain_Base'+newDir4 evalDomain_Base'+newDir5 evalDomain_Base'+newDir3 evalDomain_Base'+newDir6];

evalDomain_Full=[];
if rescaleIntervals==1
    load('transformationIntervals.mat');
    for j=1:size(evalDomain,2)
        evalDomain_Temp=(valInNewInterval(evalDomain(:,j),zeros(16,1),ones(16,1),transformationIntervals(:,1),transformationIntervals(:,2)));
        evalDomain_Full=[evalDomain_Full evalDomain_Temp];
    end
else
    evalDomain_Full=evalDomain;
end



polyChaosEval=0;
XX=Sr.knots;

numOutputs=4;
       

jacobianMatrix=[];

[modal_coeffs1,K1]=convert_to_modal(S,Sr,reduced_Outputs(1,:),domain,'legendre'); 
[modal_coeffs2,K2]=convert_to_modal(S,Sr,reduced_Outputs(2,:),domain,'legendre'); 
[modal_coeffs3,K3]=convert_to_modal(S,Sr,reduced_Outputs(3,:),domain,'legendre'); 
[modal_coeffs4,K4]=convert_to_modal(S,Sr,reduced_Outputs(4,:),domain,'legendre'); 
%[modal_coeffs5,K5]=convert_to_modal(S,Sr,reduced_Outputs(5,:),domain,'legendre'); 

output_Of_Interest=[];

for j=1:size(evalDomain_Full,2)

        
        evalDomain_Cur=evalDomain_Full(:,j);
        polyChaosEval=zeros(4,1);
        
       for i=1:length(modal_coeffs1)
           
           polyChaosEval(1)=polyChaosEval(1)+modal_coeffs1(i)*lege_eval_multidim(evalDomain_Cur,K1(i,:),domain(1,:),domain(2,:));
           polyChaosEval(2)=polyChaosEval(2)+modal_coeffs2(i)*lege_eval_multidim(evalDomain_Cur,K2(i,:),domain(1,:),domain(2,:));
           polyChaosEval(3)=polyChaosEval(3)+modal_coeffs3(i)*lege_eval_multidim(evalDomain_Cur,K3(i,:),domain(1,:),domain(2,:));
           polyChaosEval(4)=polyChaosEval(4)+modal_coeffs4(i)*lege_eval_multidim(evalDomain_Cur,K4(i,:),domain(1,:),domain(2,:));
           %polyChaosEval(5)=polyChaosEval(5)+modal_coeffs5(i)*lege_eval_multidim(evalDomain_Cur,K5(i,:),domain(1,:),domain(2,:));

       end
       
       output_Of_Interest=[output_Of_Interest  polyChaosEval];

end

output_Of_Interest-output_Of_Interest(:,1)

checkSol= stdMat*UApx*output_Of_Interest+rowMeans(:,1:size(output_Of_Interest,2));
round(checkSol)
