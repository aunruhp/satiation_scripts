function [ValRatOut, ValRatIn, SalRatOut, SalRatIn, ValAFCOut, ValA, ValB, SalAFCOut, SalA, SalB, ValDiff, SalDiff, ValSum, ValSumSal, ValSumAbs, ValChos, ValUnchos, SalChos, SalUnchos, ValAFCOutChan, ValAChan, ValBChan, SalAFCOutChan, SalAChan, SalBChan, ValDiffChan, SalDiffChan, ValSumChan, ValSumSalChan, ValSumAbsChan, ValChosChan, ValUnchosChan, SalChosChan, SalUnchosChan]=anacorr_pool


[FileName,PathName] = uigetfile('.mat','Multiselect','on');
cd (PathName)
load (FileName)

% RespRat=pvalRat(:,1)<0.05;

[ValRatOut(:,1), ValRatIn(:,1), SalRatOut(:,1), SalRatIn(:,1)]= corr_Rat(pvalRat);
[ValRatOut(:,2), ValRatIn(:,2), SalRatOut(:,2), SalRatIn(:,2)]= corr_Rat(pvalLikRat);

[ValAFCOut(:,1), ValA(:,1), ValB(:,1), SalAFCOut(:,1), SalA(:,1), SalB(:,1), ValDiff(:,1), SalDiff(:,1), ValSum(:,1), ValSumSal(:,1), ValSumAbs(:,1), ValChos(:,1), ValUnchos(:,1), SalChos(:,1), SalUnchos(:,1)]= corr_AFC(pvalA1);
[ValAFCOut(:,2), ValA(:,2), ValB(:,2), SalAFCOut(:,2), SalA(:,2), SalB(:,2), ValDiff(:,2), SalDiff(:,2), ValSum(:,2), ValSumSal(:,2), ValSumAbs(:,2), ValChos(:,2), ValUnchos(:,2), SalChos(:,2), SalUnchos(:,2)]= corr_AFC(pvalB);
[ValAFCOut(:,3), ValA(:,3), ValB(:,3), SalAFCOut(:,3), SalA(:,3), SalB(:,3), ValDiff(:,3), SalDiff(:,3), ValSum(:,3), ValSumSal(:,3), ValSumAbs(:,3), ValChos(:,3), ValUnchos(:,3), SalChos(:,3), SalUnchos(:,3)]= corr_AFC(pvalA2);
[ValAFCOut(:,4), ValA(:,4), ValB(:,4), SalAFCOut(:,4), SalA(:,4), SalB(:,4), ValDiff(:,4), SalDiff(:,4), ValSum(:,4), ValSumSal(:,4), ValSumAbs(:,4), ValChos(:,4), ValUnchos(:,4), SalChos(:,4), SalUnchos(:,4)]= corr_AFC(pvalLikAFC);

[ValAFCOutChan(:,1), ValAChan(:,1), ValBChan(:,1), SalAFCOutChan(:,1), SalAChan(:,1), SalBChan(:,1), ValDiffChan(:,1), SalDiffChan(:,1), ValSumChan(:,1), ValSumSalChan(:,1), ValSumAbsChan(:,1), ValChosChan(:,1), ValUnchosChan(:,1), SalChosChan(:,1), SalUnchosChan(:,1)]= corr_AFC(pvalA1_B);
[ValAFCOutChan(:,2), ValAChan(:,2), ValBChan(:,2), SalAFCOutChan(:,2), SalAChan(:,2), SalBChan(:,2), ValDiffChan(:,2), SalDiffChan(:,2), ValSumChan(:,2), ValSumSalChan(:,2), ValSumAbsChan(:,2), ValChosChan(:,2), ValUnchosChan(:,2), SalChosChan(:,2), SalUnchosChan(:,2)]= corr_AFC(pvalB_A2);
[ValAFCOutChan(:,3), ValAChan(:,3), ValBChan(:,3), SalAFCOutChan(:,3), SalAChan(:,3), SalBChan(:,3), ValDiffChan(:,3), SalDiffChan(:,3), ValSumChan(:,3), ValSumSalChan(:,3), ValSumAbsChan(:,3), ValChosChan(:,3), ValUnchosChan(:,3), SalChosChan(:,3), SalUnchosChan(:,3)]= corr_AFC(pvalA1_Babs);
[ValAFCOutChan(:,4), ValAChan(:,4), ValBChan(:,4), SalAFCOutChan(:,4), SalAChan(:,4), SalBChan(:,4), ValDiffChan(:,4), SalDiffChan(:,4), ValSumChan(:,4), ValSumSalChan(:,4), ValSumAbsChan(:,4), ValChosChan(:,4), ValUnchosChan(:,4), SalChosChan(:,4), SalUnchosChan(:,4)]= corr_AFC(pvalB_A2abs);


nothing=0;



function [ValRatOut, ValRatIn, SalRatOut, SalRatIn]= corr_Rat(pvalRat)

sigacc=0.05;


ValRatOut(:,1)= pvalRat(:,1)<sigacc & pvalRat(:,2)>sigacc;
ValRatIn(:,1)=  pvalRat(:,3)<sigacc & pvalRat(:,4)>sigacc;
SalRatOut(:,1)= pvalRat(:,5)<sigacc & pvalRat(:,6)>sigacc;
SalRatIn(:,1)=  pvalRat(:,7)<sigacc & pvalRat(:,8)>sigacc;
 



function [ValAFCOut, ValA, ValB, SalAFCOut, SalA, SalB, ValDiff, SalDiff, ValSum, ValSumSal, ValSumAbs, ValChos, ValUnchos, SalChos, SalUnchos]= corr_AFC(pvalAFC)

sigacc=0.05;
sigA2=0.05;

ValAFCOut(:,1)= pvalAFC(:,1)<sigacc &  pvalAFC(:,2) >sigacc;
ValA(:,1)=      pvalAFC(:,3)<sigacc &  pvalAFC(:,4) >sigacc;
ValB(:,1)=      pvalAFC(:,5)<sigacc &  pvalAFC(:,6) >sigacc;
SalAFCOut(:,1)= pvalAFC(:,7)<sigacc &  pvalAFC(:,8) >sigacc;
SalA(:,1)=      pvalAFC(:,9)<sigacc &  pvalAFC(:,10)>sigacc;
SalB(:,1)=      pvalAFC(:,11)<sigacc & pvalAFC(:,12)>sigacc;
ValDiff(:,1)=   pvalAFC(:,13)<sigacc & pvalAFC(:,14)>sigacc;
SalDiff(:,1)=   pvalAFC(:,15)<sigacc & pvalAFC(:,16)>sigacc;
ValSum(:,1)=    pvalAFC(:,17)<sigacc & pvalAFC(:,18)>sigacc;
ValSumSal(:,1)= pvalAFC(:,19)<sigacc & pvalAFC(:,20)>sigacc;
ValSumAbs(:,1)= pvalAFC(:,21)<sigacc & pvalAFC(:,22)>sigacc; 
ValChos(:,1)=   pvalAFC(:,23)<sigacc & pvalAFC(:,24)>sigacc;
ValUnchos(:,1)= pvalAFC(:,25)<sigacc & pvalAFC(:,26)>sigacc;
SalChos(:,1)=   pvalAFC(:,27)<sigacc & pvalAFC(:,28)>sigacc;
SalUnchos(:,1)= pvalAFC(:,29)<sigacc & pvalAFC(:,30)>sigacc;

