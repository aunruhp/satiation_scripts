function [ZPF, p_ZPF_2tail, p_ZPF_1tail] = ZPF_test(x,y, type) 

%% Corraltion-Matrix for modified Pearson-Filon statisti (ZPF)
    % include crossed correlations for a correlation-matrix
   
    
%% Random effect analysis
% in this analysis first the contrast is made and afterwards the
% second-level analysis. One common way to compare 2 R-pearson correlations
% is Fisher-Ztest, which can yield the statistic and CI as well, and the
% SD. With this I could d second level analysis. However it is for
% non-reoeated measures...
% solution Comparing Correlated but Nonoverlapping Correlation Coefficients
% allows comparing 2 vaariables in 2 different time points
%
% modified Pearson-Filon statistic, the ZPF statistic is superior to the PF statistic:

% Raghunathan, T. E., Rosenthal, R., and Rubin, D. (1996). Comparing correlated but nonoverlapping correlations. Psychol. Methods 1, 178–183
% which is a recommended test statistic for two, non-overlapping dependent correlations (Krishnamoorthy and Xia, 2007).
%code copied from http://core.ecu.edu/psyc/wuenschk/SAS/SAS-Programs.htm.
%some modifications
% this tests allow to measure the correlation of 2 variables in 2 time
% points to see if it changes, this is an example of non-overlapping
% correlation

% orginal paper of the method Dunn, O.J.,  &  Clark,  V. A.  (1969).  Correlation   coeffi-cients  measured  on  the  same  individuals. Journal  ofthe American  Statistical  Association, 64, 366-377
% Once the levels  were  adjusted, the  power  of ZPF   was  only slightly better  than that  of  PF.   Nev- ertheless,   the   main advantage  of   using  ZPF   is  that the levels   do  not   need  adjustment  
% Pair An identification number for the pair of correlations being compared
% N    The number of subjects.
% A RT time 1
% B value time 1
% X RT time 2
% Y value time 2
% AB   Correlation between variables A and B at time 1.
% AX   Correlation between variable A at time 1 and A at time 2.
% AY   Correlation between variable A at time 1 and B at time 2.
% BX   Correlation between variable B at time 1 and A at time 2.
% BY   Correlation between variable B at time 1 and B at time 2.
% XY   Correlation between variables A and B at time 2.

%  We want to compare AB with XY.     
    
% for independent correlations, between subjects for example, Fisher
% transform of r to Z, and then Z test. But here are dependent, within!
  
  

A=[x, y];

if type==0
    B=corrcoef(A);
elseif type==1
    B=spearmanMatrix(A);
end


AX= B(1,2); AB=B(1,3); AY=B(1,4);
BX= B(2,3); XY=B(2,4); BY=B(3,4);
N=size(A,1); n=N; % this test is used for persons, thus, each person is a datapoint in a vector of each variabel ,in my case each stimulus would be a person


k = (AX-BX*AB)*(BY-BX*XY)+(AY-AX*XY)*(BX-AX*AB)+(AX-AY*XY)*(BY-AY*AB)+(AY-AB*BY)*(BX-BY*XY);
% PF = (AB-XY)*SQRT(N)/SQRT((1-AB**2)**2+(1-XY**2)**2-k);
PF = (AB-XY)*sqrt(N)/sqrt((1-AB^2)^2+(1-XY^2)^2-k);
% p_PF = 2*probnorm(-1*ABS(PF));
% p_PF = 2*normcdf(-1*ABS(PF));
p_PF = 2*normcdf(-1*abs(PF));
ZAB = .5*log((1+AB)/(1-AB)); % fisher Z transformation
ZXY = .5*log((1+XY)/(1-XY)); % fisher Z transformation
% ZPF = sqrt((n-3)/2)*(ZAB-ZXY)/sqrt(1-(k/(2*(1-AB**2)*(1-XY**2))));
ZPF = sqrt((n-3)/2)*(ZAB-ZXY)/sqrt(1-(k/(2*(1-AB^2)*(1-XY^2))));
% p_ZPF = 2*probnorm(-1*ABS(ZPF));
p_ZPF_2tail = 2*normcdf(-1*abs(ZPF)); % 2-tailed pval

if ZPF>0
    p_ZPF_1tail = 1-normcdf(ZPF); % 2-tailed pval
elseif ZPF<0
    p_ZPF_1tail = normcdf(ZPF); % 2-tailed pval
end



