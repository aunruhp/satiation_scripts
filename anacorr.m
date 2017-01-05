function anacorr

cd '/media/Projects/Alex/New Analysis/corrmatrices2'
load Tot_correlations_0_1000_0_2000_0_2000.mat

% RespRat=pvalRat(:,1)<0.05;

[ValOutRat, ValInRat, SalOutRat, SalInRat]= anacorr_Rat(pvalRat);
% [ValOutLikRat, ValInLikRat, SalOutLikRat, SalInLikRat]= anacorr_Rat(pvalLikRat);
% 
% [ValOutA1,ValA_A1,ValB_A1,ValAbsPrefA1,SalA_A1,SalB_A1,Diff_A1,AbsDiff_A1,Sum_A1,SalSum_A1,AbsSum_A1, ValC_A1,ValU_A1,SalC_A1,SalU_A1] = anacorr_AFC(pvalA1);
% [ValOutB,  ValA_B, ValB_B, ValAbsPrefB, SalA_B, SalB_B, Diff_B, AbsDiff_B, Sum_B, SalSum_B, AbsSum_B, ValC_B,ValU_B,SalC_B,SalU_B] = anacorr_AFC(pvalB);
% [ValOutA2,ValA_A2,ValB_A2,ValAbsPrefA2,SalA_A2,SalB_A2,Diff_A2,AbsDiff_A2,Sum_A2,SalSum_A2,AbsSum_A2, ValC_A2,ValU_A2,SalC_A2,SalU_A2] = anacorr_AFC(pvalA2);
% [ValOutLik, ValA_Lik, ValB_Lik, ValAbsPrefLik, SalA_Lik, SalB_Lik, Diff_Lik, AbsDiff_Lik, Sum_Lik, SalSum_Lik, AbsSum_Lik, ValC_lik,ValU_lik,SalC_lik,SalU_lik] = anacorr_AFC(pvalLikAFC);
% 
% [ValOutA1_B,ValA_A1_B,ValB_A1_B,ValAbsPrefA1_B,SalA_A1_B,SalB_A1_B,Diff_A1_B,AbsDiff_A1_B,Sum_A1_B,SalSum_A1_B,AbsSum_A1_B, ValC_A1_B,ValU_A1_B,SalC_A1_B,SalU_A1_B] = anacorr_AFC(pvalA1_B);
% [ValOutB_A2,ValA_B_A2,ValB_B_A2,ValAbsPrefB_A2,SalA_B_A2,SalB_B_A2,Diff_B_A2,AbsDiff_B_A2,Sum_B_A2,SalSum_B_A2,AbsSum_B_A2, ValC_B_A2,ValU_B_A2,SalC_B_A2,SalU_B_A2] = anacorr_AFC(pvalB_A2);
% [ValOutA1_Babs,ValA_A1_Babs,ValB_A1_Babs,ValAbsPrefA1_Babs,SalA_A1_Babs,SalB_A1_Babs,Diff_A1_Babs,AbsDiff_A1_Babs,Sum_A1_Babs,SalSum_A1_Babs,AbsSum_A1_Babs, ValC_A1_Babs,ValU_A1_Babs,SalC_A1_Babs,SalU_A1_Babs] = anacorr_AFC(pvalA1_Babs);
% [ValOutB_A2abs,ValA_B_A2abs,ValB_B_A2abs,ValAbsPrefB_A2abs,SalA_B_A2abs,SalB_B_A2abs,Diff_B_A2abs,AbsDiff_B_A2abs,Sum_B_A2abs,SalSum_B_A2abs,AbsSum_B_A2abs,ValC_B_A2abs,ValU_B_A2abs,SalC_B_A2abs,SalU_B_A2abs] = anacorr_AFC(pvalB_A2abs);
% 
% 
% 
% 
% 
% 
% 
% sum(ValOutRat), sum(ValInRat), sum(SalOutRat), sum(SalInRat)
% sum(ValOutLikRat), sum(ValInLikRat), sum(SalOutLikRat), sum(SalInLikRat)
% sum(ValOutA1),sum(ValA_A1),sum(ValB_A1),sum(ValAbsPrefA1),sum(SalA_A1),sum(SalB_A1),sum(Diff_A1),sum(AbsDiff_A1),sum(Sum_A1),sum(SalSum_A1),sum(AbsSum_A1), sum(ValC_A1), sum(ValU_A1), sum(SalC_A1), sum(SalU_A1)
% sum(ValOutA2),sum(ValA_A2),sum(ValB_A2),sum(ValAbsPrefA2),sum(SalA_A2),sum(SalB_A2),sum(Diff_A2),sum(AbsDiff_A2),sum(Sum_A2),sum(SalSum_A2),sum(AbsSum_A2), sum(ValC_B), sum(ValU_B), sum(SalC_B), sum(SalU_B)
% sum(ValOutB),sum(ValA_B),sum(ValB_B),sum(ValAbsPrefB),sum(SalA_B),sum(SalB_B),sum(Diff_B),sum(AbsDiff_B),sum(Sum_B),sum(SalSum_B),sum(AbsSum_B)           , sum(ValC_A2), sum(ValU_A2), sum(SalC_A2), sum(SalU_A2)
% sum(ValOutLik),sum(ValA_Lik),sum(ValB_Lik),sum(ValAbsPrefLik),sum(SalA_Lik),sum(SalB_Lik),sum(Diff_Lik),sum(AbsDiff_Lik),sum(Sum_Lik),sum(SalSum_Lik),sum(AbsSum_Lik), sum(ValC_lik), sum(ValU_lik), sum(SalC_lik), sum(SalU_lik)
% sum(ValOutA1_B),sum(ValA_A1_B),sum(ValB_A1_B),sum(ValAbsPrefA1_B),sum(SalA_A1_B),sum(SalB_A1_B),sum(Diff_A1_B),sum(AbsDiff_A1_B),sum(Sum_A1_B),sum(SalSum_A1_B),sum(AbsSum_A1_B), sum(ValC_A1_B), sum(ValU_A1_B), sum(SalC_A1_B), sum(SalU_A1_B)
% sum(ValOutB_A2),sum(ValA_B_A2),sum(ValB_B_A2),sum(ValAbsPrefB_A2),sum(SalA_B_A2),sum(SalB_B_A2),sum(Diff_B_A2),sum(AbsDiff_B_A2),sum(Sum_B_A2),sum(SalSum_B_A2),sum(AbsSum_B_A2), sum(ValC_B_A2), sum(ValU_B_A2), sum(SalC_B_A2), sum(SalU_B_A2)
% sum(ValOutA1_Babs),sum(ValA_A1_Babs),sum(ValB_A1_Babs),sum(ValAbsPrefA1_Babs),sum(SalA_A1_Babs),sum(SalB_A1_Babs),sum(Diff_A1_Babs),sum(AbsDiff_A1_Babs),sum(Sum_A1_Babs),sum(SalSum_A1_Babs),sum(AbsSum_A1_Babs), sum(ValC_A1_Babs), sum(ValU_A1_Babs), sum(SalC_A1_Babs), sum(SalU_A1_Babs)
% sum(ValOutB_A2abs),sum(ValA_B_A2abs),sum(ValB_B_A2abs),sum(ValAbsPrefB_A2abs),sum(SalA_B_A2abs),sum(SalB_B_A2abs),sum(Diff_B_A2abs),sum(AbsDiff_B_A2abs),sum(Sum_B_A2abs),sum(SalSum_B_A2abs),sum(AbsSum_B_A2abs), sum(ValC_B_A2abs), sum(ValU_B_A2abs), sum(SalC_B_A2abs), sum(SalU_B_A2abs)

nothing=0;


function [ValOutA1, ValAInA1,ValBinA1,ValAbsPrefA1,SalA_A1,SalB_A1,Diff_A1, AbsDiff_A1, Sum_A1, SalSum_A1, AbsSum_A1,ValC_A1,ValU_A1,SalC_A1,SalU_A1] = anacorr_AFC(pvalA1)

sigacc=0.05;
sigrej=0.05;

ValOutA1= pvalA1(:,4)<sigacc; % pvalA1(:,1)<sigacc) & pvalA1(:,2)>sigrej & pvalA1(:,3)>sigrej &
ValAInA1= pvalA1(:,8)<sigacc; % pvalA1(:,5)<sigacc) &  pvalA1(:,6)>sigrej & pvalA1(:,7)>sigrej &
ValBinA1=  pvalA1(:,12)<sigacc;% pvalA1(:,9)<sigacc)  & pvalA1(:,10)>sigrej & pvalA1(:,11)>sigrej & 
ValAbsPrefA1=  pvalA1(:,16)<sigacc; % pvalA1(:,13)<sigacc)  & pvalA1(:,14)>sigrej & pvalA1(:,15)>sigrej &
SalA_A1=  pvalA1(:,20)<sigacc; % pvalA1(:,17)<sigacc)  & pvalA1(:,18)>sigrej & pvalA1(:,19)>sigrej &
SalB_A1=  pvalA1(:,24)<sigacc; % pvalA1(:,21)<sigacc) & pvalA1(:,22)>sigrej & pvalA1(:,23)>sigrej &
Diff_A1= pvalA1(:,28)<sigacc; % pvalA1(:,25)<sigacc) & pvalA1(:,26)>sigrej & pvalA1(:,27)>sigrej &
AbsDiff_A1=  pvalA1(:,32)<sigacc; % pvalA1(:,29)<sigacc) & pvalA1(:,30)>sigrej & pvalA1(:,31)>sigrej & 
Sum_A1=  pvalA1(:,36)<sigacc; % pvalA1(:,33)<sigacc) & pvalA1(:,34)>sigrej & pvalA1(:,35)>sigrej &
SalSum_A1=  pvalA1(:,40)<sigacc;% pvalA1(:,37)<sigacc) & pvalA1(:,38)>sigrej & pvalA1(:,39)>sigrej &
AbsSum_A1=  pvalA1(:,44)<sigacc; % pvalA1(:,41)<sigacc) & pvalA1(:,42)>sigrej & pvalA1(:,43)>sigrej &
ValC_A1= (pvalA1(:,45)<0.05) & pvalA1(:,46)>0.05 & pvalA1(:,47)>0.05 & pvalA1(:,48)<0.05;
ValU_A1= (pvalA1(:,49)<0.05) & pvalA1(:,50)>0.05 & pvalA1(:,51)>0.05 & pvalA1(:,52)<0.05;
SalC_A1= (pvalA1(:,53)<0.05) & pvalA1(:,54)>0.05 & pvalA1(:,55)>0.05 & pvalA1(:,56)<0.05;
SalU_A1= (pvalA1(:,57)<0.05) & pvalA1(:,58)>0.05 & pvalA1(:,59)>0.05 & pvalA1(:,60)<0.05;



function [ValOutRat, ValInRat, SalOutRat, SalInRat]= anacorr_Rat(pvalRat)

sigacc=0.05;
sigrej=0.05;

ValOutRat= pvalRat(:,4)<sigacc; %(pvalRat(:,1)<sigacc) & pvalRat(:,2)>sigrej & pvalRat(:,3)>sigrej & 
ValInRat=  pvalRat(:,8)<sigacc; %(pvalRat(:,5)<sigacc) & pvalRat(:,6)>sigrej & pvalRat(:,7)>sigrej &
SalOutRat= pvalRat(:,12)<sigacc;%(pvalRat(:,9)<sigacc) & pvalRat(:,10)>sigrej & pvalRat(:,11)>sigrej & 
SalInRat= pvalRat(:,16)<sigacc; %(pvalRat(:,13)<sigacc) & pvalRat(:,14)>sigrej & pvalRat(:,15)>sigrej & 

