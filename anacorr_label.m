function [ValRat, ValB, ValA1, ValA2, ValLik, ValRA, ValRAB, ValRABA, ValRLik, ValAB, ValABA, ValAA, SalRat, SalB, SalA1, SalA2, SalLik, SalRA, SalRAB, SalRABA, SalRLik, SalAB, SalABA, SalAA]=anacorr_label


[FileName,PathName] = uigetfile('.mat','Multiselect','on');
cd (PathName)
load (FileName)

% RespRat=pvalRat(:,1)<0.05;

[ValRat, ValB, ValA1, ValA2, ValLik, ValRA, ValRAB, ValRABA, ValRLik, ValAB, ValABA, ValAA]= anacorr_Label(pvalLabelVal);
[SalRat, SalB, SalA1, SalA2, SalLik, SalRA, SalRAB, SalRABA, SalRLik, SalAB, SalABA, SalAA]= anacorr_Label(pvalLabelSal);


nothing=0;



function [ValRat, ValB, ValA1, ValA2, ValLik, ValRA, ValRAB, ValRABA, ValRLik, ValAB, ValABA, ValAA]= anacorr_Label(pvalRat)

sigacc=0.06;
sigA2=0.1;

ValRat(:,1)= pvalRat(:,1)<sigacc; %& pvalRat(:,2)<sigacc;
ValB(:,1)=   pvalRat(:,3)<sigacc; %& pvalRat(:,4)<sigacc;
ValA1(:,1)=  pvalRat(:,5)<sigacc; %& pvalRat(:,6)<sigacc;
ValA2(:,1)=  pvalRat(:,7)<sigacc; %& pvalRat(:,8)<sigacc;
ValLik(:,1)= pvalRat(:,9)<sigacc; %& pvalRat(:,10)<sigacc;
ValRA(:,1)=  pvalRat(:,11)<sigacc; %& pvalRat(:,12)<sigacc;
ValRAB(:,1)= pvalRat(:,13)<sigacc; %& pvalRat(:,14)<sigacc;
ValRABA(:,1)= pvalRat(:,15)<sigacc; %& pvalRat(:,16)<sigacc;
ValRLik(:,1)= pvalRat(:,17)<sigacc; %& pvalRat(:,18)<sigacc;
ValAB(:,1)=  pvalRat(:,19)<sigacc; %& pvalRat(:,20)<sigacc;
ValABA(:,1)= pvalRat(:,21)<sigacc; %& pvalRat(:,22)<sigacc; 
ValAA(:,1)=  pvalRat(:,23)<sigacc; %& pvalRat(:,24)<sigacc;

ValRat(:,2)= pvalRat(:,2)<sigacc; %& pvalRat(:,2)<sigacc;
ValB(:,2)=   pvalRat(:,4)<sigacc; %& pvalRat(:,4)<sigacc;
ValA1(:,2)=  pvalRat(:,6)<sigacc; %& pvalRat(:,6)<sigacc;
ValA2(:,2)=  pvalRat(:,8)<sigacc; %& pvalRat(:,8)<sigacc;
ValLik(:,2)= pvalRat(:,10)<sigacc; %& pvalRat(:,10)<sigacc;
ValRA(:,2)=  pvalRat(:,12)<sigacc; %& pvalRat(:,12)<sigacc;
ValRAB(:,2)= pvalRat(:,14)<sigacc; %& pvalRat(:,14)<sigacc;
ValRABA(:,2)= pvalRat(:,16)<sigacc; %& pvalRat(:,16)<sigacc;
ValRLik(:,2)= pvalRat(:,18)<sigacc; %& pvalRat(:,18)<sigacc;
ValAB(:,2)=  pvalRat(:,20)<sigacc; %& pvalRat(:,20)<sigacc;
ValABA(:,2)= pvalRat(:,22)<sigacc; %& pvalRat(:,22)<sigacc; 
ValAA(:,2)=  pvalRat(:,24)<sigacc; %& pvalRat(:,24)<sigacc;

ValRat(:,3)= pvalRat(:,25)<sigacc; 
ValA1(:,3)=  pvalRat(:,27)<sigacc; 
ValB(:,3)=   pvalRat(:,26)<sigacc; 
ValA2(:,3)=  pvalRat(:,28)<sigacc; 
ValLik(:,3)= pvalRat(:,29)<sigacc; 
ValRA(:,3)=  pvalRat(:,30)<sigacc; 
ValRAB(:,3)= pvalRat(:,31)<sigacc; 
ValRABA(:,3)= pvalRat(:,32)<sigacc; 
ValRLik(:,3)= pvalRat(:,33)<sigacc; 
ValAB(:,3)=  pvalRat(:,34)<sigacc; 
ValABA(:,3)= pvalRat(:,35)<sigacc; 
ValAA(:,3)=  pvalRat(:,36)<sigacc; 
