dominio=domain;
%baseDomain=domain(1,:);
%baseDomain(end-1:end)=domain(2,end-1:end);
%baseDomain(end-5:end-4)=domain(2,end-5:end-4);
%dominio(:,end-1:end)=flipud(dominio(:,end-1:end));
%dominio(:,end-5:end-4)=flipud(dominio(:,end-5:end-4));

rescaleIntervals=0;


mmm=[          


         -7.59312365836971
         -7.31859778778378
          6.09893082455153
         -1.41493907814052

];


evalDomain_Base=[0.5 0.5 0.5 0.5 0.5];

%newDir=Vk*(Sk\Uk')*[0; 0; 0  ;0;];
newDir=Vk*(Sk\Uk')*mmm;

evalDomain=[evalDomain_Base' evalDomain_Base'+newDir];

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

output_Of_Interest=[];

for j=1:size(evalDomain_Full,2)

        
        evalDomain_Cur=evalDomain_Full(:,j);
        polyChaosEval=zeros(4,1);
        
       for i=1:length(modal_coeffs1)
            
           polyChaosEval(1)=polyChaosEval(1)+modal_coeffs1(i)*lege_eval_multidim(evalDomain_Cur,K1(i,:),domain(1,:),domain(2,:));
           polyChaosEval(2)=polyChaosEval(2)+modal_coeffs2(i)*lege_eval_multidim(evalDomain_Cur,K2(i,:),domain(1,:),domain(2,:));
           polyChaosEval(3)=polyChaosEval(3)+modal_coeffs3(i)*lege_eval_multidim(evalDomain_Cur,K3(i,:),domain(1,:),domain(2,:));
           polyChaosEval(4)=polyChaosEval(4)+modal_coeffs4(i)*lege_eval_multidim(evalDomain_Cur,K4(i,:),domain(1,:),domain(2,:));
               
       end
       
       output_Of_Interest=[output_Of_Interest  polyChaosEval];

end

output_Of_Interest-output_Of_Interest(:,1)


