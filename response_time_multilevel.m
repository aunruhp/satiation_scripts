function [H,P,CI, stats]=response_time_multilevel

dbstop if error

type=0;  % type 0 means Pearsons matrix of correlation, if 1 means Spearmans

f={};
try
    cd('/media/Projects/Alex/Satiation Analysis');      %/Volumes/MTL/CS4 output')   %orginal form michael '/Volumes/Barnacle/Data/Satiation_data/Satiation_session_data'
catch
    cd('/media/Projects/Alex/');      %/Volumes/MTL/CS4 output')
end

DirName= dir;
for a = 1:length(DirName)
    if regexp(DirName(a).name, 'Sat_+[0-9][0-9][0-9]')
               f{end+1,1} = [DirName(a).name];
    end
end



for i=1:length(f)
   
        [rho2AFCabs{i}, rhoRat{i}, pval2AFCabs{i}, pvalRat{i}, ZPF{i}, pvalZPF{i}, ZPFRat{i}, pvalZPFrat{i},lm_2afc{i}, lm_rat{i}, lm_interact{i}, lm_MagValen{i}, lm_salience{i}, tbetas{i}, pbetas{i}, RT{i}, val{i},  pottlm{i}, compFval{i},interPval{i}, RT2AFC{i}, diff_vall2AFC{i}] = response_time_ana(f{i},type);

end


h=1; 
for i=[1:7, 9:11]     % skip patient 39 because he did not eat food, not used in any test of this type

    % 2AFC correlations
    
    for j=1:length(ZPF{i})
%        eval(['t1' num2str(j) '(h) = ZPF{i}(j);'])
       eval(['t2' num2str(j) '(h) = rho2AFCabs{i}(j,1);'])
       eval(['t3' num2str(j) '(h) = rho2AFCabs{i}(j,2);'])
    end
    
    
%     
%     for j=1:length(tbetas{i})
%         
%      % 2AFC Regression
% %        eval(['t4' num2str(j) '(h) = tbetas{i}(j) ;'])
%     end
    
    for j=1:2
     
     % 2AFC Regression

       eval(['t6' num2str(j) '(h) = lm_2afc{i}{j}.Coefficients.Estimate(2);'])
     
     % Rating Regression  

       eval(['t9' num2str(j) '(h) = lm_rat{i}{j}.Coefficients.Estimate(2);'])
       eval(['t10' num2str(j) '(h) = lm_rat{i}{j}.Coefficients.Estimate(3);'])
       eval(['t14' num2str(j) '(h) = lm_interact{i}{j}.Coefficients.Estimate(4);'])
 
     % Magnitude and Valence Regression  
       eval(['t16' num2str(j) '(h) = lm_MagValen{i}{j}.Coefficients.Estimate(2);'])

       eval(['t18' num2str(j) '(h) = lm_MagValen{i}{j}.Coefficients.Estimate(3);'])

       eval(['t20' num2str(j) '(h) = lm_MagValen{i}{j}.Coefficients.Estimate(4);'])
     
     %  only salience, also called magnitude
     
       eval(['t25' num2str(j) '(h) = lm_salience{i}{j}.Coefficients.Estimate(2);'])
     
     % Rating correlations, 1 difference then each correlation
      
%        eval(['t11' num2str(j) '(h) = ZPFRat{i}(j);'])
       eval(['t12' num2str(j) '(h) = rhoRat{i}(j,1);'])
       eval(['t13' num2str(j) '(h) = rhoRat{i}(j,2);'])
       

%      sum only in the model if difference not included, moreoer, quite
%      collinear in collintest, and correlted 0.6, way way more than I
%      expected

    end
    h=h+1;
end


% 2AFC Correlation tests

for i=1:length(ZPF{1}) 

       
        eval(['[P_2corrB{i},H_2corrB{i}] =ttest(t2' num2str(i) ');']);
        eval(['[P_2corrA{i},H_2corrA{i}] =ttest(t3' num2str(i) ');']);
        eval(['[P_2corrPairedT{i},H_2corrPairedT{i}] =ttest(t2' num2str(i) ', t3' num2str(i) ');']); % paired ttest of Spearmann correlation
end


% Rating correlation tests

for i=1:2
   
%     eval(['[P_RcorrAB{i}, H_RcorrAB{i}]=ttest(t11' num2str(i) ');']);
    eval(['[P_RcorrB{i}, H_RcorrB{i}]  =ttest(t12' num2str(i) ');']);  % group level correlation
    eval(['[P_RcorrA{i}, H_RcorrA{i}]  =ttest(t13' num2str(i) ');']);
    eval(['[P_RcorrPairedT{i}, H_RcorrPairedT{i}]  =ttest(t12' num2str(i) ', t13' num2str(i) ');']);  % paired ttest of Spearmann correlation
    
end



% Regression tests

% for i=1:length(tbetas{1})
%     
%        eval(['[H_betasAB{i},P_betasAB{i},CI_betasAB{i}, stats_betasAB{i}]=ttest(t4' num2str(i) ');']);
% 
% end


% Group lelvel Analysis out of subjectwise regression

for i=1:2    
    
    % Difference value in 2AFC

    eval(['[H_2beta{i}, P_2beta{i}]=ttest(t6' num2str(i) ');']);
    
    % Value and Salience in rating

    eval(['[H_RbetaVal{i}, HPRbetaVal{i}]=ttest(t9' num2str(i) ');']);
    eval(['[H_RbetaSal{i}, P_RbetaSal{i}]=ttest(t10' num2str(i) ');']);
    eval(['[H_RbetaInteraction{i}, P_RbetaInteraction{i}]=ttest(t14' num2str(i) ');']);
    eval(['[H_RbetaSalienceOnly{i}, P_RbetaSalienceOnly{i}]=ttest(t25' num2str(i) ');']);
    
    % Magnitude and valence in rating
    
    eval(['[H_betaMag{i}, P_betaMag{i}]=ttest(t16' num2str(i) ');']);
    eval(['[H_betaValence{i}, P_betaValence{i}]=ttest(t18' num2str(i) ');']);
    eval(['[H_betaMagValInter{i}, P_betaMagValInter{i}]=ttest(t20' num2str(i) ');']);

end

% Pair ttests of betas, group level difference before vs after
     
[P_2betaPairedT, H_2betaPairedT]=ttest(t61, t62);
[P_RbetaValPairedT, H_RbetaValPairedT]=ttest(t91, t92);
[P_RbetaSalPairedT, H_RbetaSalPairedT]=ttest(t101, t102);
[P_betaMagPairedT, H_betaMagPairedT]=ttest(t161,t162);
[P_betaValencePairedT, H_betaValencePairedT]=ttest(t181,t182);



%% Potthoff analysis of Before vs After

% this are repeated measures, so, using this to make test of Dependent
% variables. However, Potthoff is designed for 2 independent measures
% for the comparisson, this increasis artificially the DOF, but even 
% though negative, so more sensitive and nothing found

i=1;
for h=[1:7, 9:11]
    
Diff_Potthoffpvals{i}=pottlm{h}{1}{2}.Coefficients.pValue(3:4);
Value_Potthoffpvals{i}=pottlm{h}{2}{2}.Coefficients.pValue(3:4);
Sal_Potthoffpvals{i}=pottlm{h}{3}{2}.Coefficients.pValue(3:4);
Ftestpval_diff{i}=compFval{h}{1}(2);
Ftestpval_val{i}=compFval{h}{2}(2);
Ftestpval_sal{i}=compFval{h}{3}(2);

betaInter_diff(i)=pottlm{h}{1}{2}.Coefficients.Estimate(3);
betaSlope_diff(i)=pottlm{h}{1}{2}.Coefficients.Estimate(4);
betaInter_val(i)=pottlm{h}{2}{2}.Coefficients.Estimate(3);
betaSlope_val(i)=pottlm{h}{2}{2}.Coefficients.Estimate(4);
betaInter_sal(i)=pottlm{h}{3}{2}.Coefficients.Estimate(3);
betaSlope_sal(i)=pottlm{h}{3}{2}.Coefficients.Estimate(4);
i=i+1;
end

[H_ttestPottDiffInter, P_ttestPottDiffInter]=ttest(betaInter_diff);  % this one is positive, but only 2 patients had signifficant change in intercept though
[H_ttestPottDiffSlope, P_ttestPottDiffSlope]=ttest(betaSlope_diff);

[H_ttestPottValInter, P_ttestPottValInter]=ttest(betaInter_val);
[H_ttestPottValSlope, P_ttestPottValSlope]=ttest(betaSlope_val);

[H_ttestPottSalInter, P_ttestPottSalInter]=ttest(betaInter_sal);
[H_ttestPottSalSlope, P_ttestPottSalSlope]=ttest(betaSlope_sal); % this one is positive, but only 2 patients had signifficant change in intercept though

% if slope changes, then different effect on RT before an after, 
% if intercept changes not very meaningful.. Or maybe it could 
% make harder interpretation


%% Potthoff analsyis of Magnitude and Valence

% in this case, it is for the analysis of one dataset only, so Before and
% After are NOT mixed for regression. In first step magnitude is incldued, 
% and in a second step, 
i=1;
for h=[1:7, 9:11]
    
ValenceBefore_Potthoffpvals{i}=pottlm{h}{4}{2}.Coefficients.pValue(3:4);
ValeceAfter_Potthoffpvals{i}=pottlm{h}{4}{4}.Coefficients.pValue(3:4);

Ftestpval_valenceBefore{i}=compFval{h}{4}{1}(2);
Ftestpval_valenceAfter{i}=compFval{h}{4}{2}(2);

MagnitudeBetaBef(i)=pottlm{h}{4}{1}.Coefficients.Estimate(2);  % betas of magnitude when alone, so the same  as I wanted to do with salience alone
MagnitudeBetaAft(i)=pottlm{h}{4}{3}.Coefficients.Estimate(2);

betaInter_valenBef(i)=pottlm{h}{4}{2}.Coefficients.Estimate(3);
betaSlope_valenBef(i)=pottlm{h}{4}{2}.Coefficients.Estimate(4);
betaInter_valenAft(i)=pottlm{h}{4}{4}.Coefficients.Estimate(3);
betaSlope_valenAft(i)=pottlm{h}{4}{4}.Coefficients.Estimate(4);
i=i+1;
end


[H_ttestPott_ValenBefInter, P_ttestPott_ValenBefInter]=ttest(betaInter_valenBef);  
[H_ttestPott_ValenBefSlope, P_ttestPott_ValenBefSlope]=ttest(betaSlope_valenBef);

[H_ttestPott_ValenAftInter, P_ttestPott_ValenAftInter]=ttest(betaInter_valenAft);
[H_ttestPott_ValenAftSlope, P_ttestPott_ValenAftSlope]=ttest(betaSlope_valenAft);


% Using Pothoff analysis, Valence does not enter in the model, the slope is
% not significant!

[H_MagnitudeBetaAloneBef, P_MagnitudeBetaAloneAft]=ttest(MagnitudeBetaBef);  
[H_MagnitudeBetaAloneAft, P_MagnitudeBetaAloneAft]=ttest(MagnitudeBetaAft);
[H_pairedTMagnitudeAlone, P_pairedTMagnitudeAlone]=ttest(MagnitudeBetaBef,MagnitudeBetaAft);

%& if magnitude alone  in the linear model, it is significant in Before, after, and no differece between them in paired Ttest


%% ANOVA simple  This is completely Fixed effects though

% RT=RT2AFC; val=diff_vall2AFC;


RTpar=[];
valpar=[];
Valence=[];
PID=[];
Magnitude=[];

for i=1:length(val)
    
   RTpar(end+1:end+length(RT{1}),1)=RT{i}(:,1);
   valpar(end+1:end+length(RT{1}),1)=val{i}(:,1);

   Valence(end+1:end+length(RT{1}),1)=val{i}(:,1)<0; % pos or neg
   PID(end+1:end+length(RT{1}),1)= i;
   
   RTpar(end-length(RT{1})+1:end,2)=RT{i}(:,2);
   valpar(end-length(RT{1})+1:end,2)=val{i}(:,2);
                           %(25, 50, 75, 100)
   Valence(end-length(RT{1})+1:end,2)=val{i}(:,2)<0; % pos or neg
   PID(end-length(RT{1})+1:end,2)= i;
end
Magnitude(:,1)=abs(valpar(:,1));                          %(25, 50, 75, 100)
Magnitude(:,2)=abs(valpar(:,2)); 


for i=1:2
    
    [P_anova{i},T_anova{i},STATS_anova{i},TERMS_anova{i}]=anovan(RTpar(:,i), {Magnitude(:,i), Valence(:,i)}, 'varnames', {'Mag', 'Valence'});
    [P_anovaInt{i},T_anovaInt{i},STATS_anovaInt{i},TERMS_anovaInt{i}]=anovan(RTpar(:,i), {Magnitude(:,i), Valence(:,i)}, 'varnames', {'Mag', 'Valence'}, 'model', 'interaction');
end




%% LME Model for each

% tried to make RManova as well, but, 1st, each condition is repeated many
% times and a different number of times for each condition, completely
% unbalanced, moreover some patients have 'missing data', not all patients 
% give answers for each of the 8 possible conditins, so as missing data

% First, trying to analyze the effect before and after

for i=1:2
tbl_lme= table(RTpar(:,i) , Magnitude(:,i) , Valence(:,i) , PID(:,i), 'VariableNames', {'RT', 'Magnitude', 'Valence' , 'PID'} );
lme_model{i}= fitlme(tbl_lme,  'RT~ 1+Magnitude+Valence + (Valence+Magnitude|PID)');
lme_model2{i}= fitlme(tbl_lme, 'RT~ 1+Magnitude+Valence + (1|PID)+(Magnitude-1|PID)+(Valence-1|PID)');
lme_model3{i}= fitlme(tbl_lme, 'RT~ 1+Magnitude+Valence + (Magnitude|PID) + (Valence-1|PID)');
lme_model4{i}= fitlme(tbl_lme, 'RT~ 1+Magnitude*Valence + (Valence*Magnitude|PID)');
lme_model5{i}= fitlme(tbl_lme, 'RT~ 1+Magnitude*Valence + (1|PID)+(Magnitude-1|PID)+(Valence-1|PID)');
lme_model6{i}= fitlme(tbl_lme, 'RT~ 1+Magnitude*Valence + (Magnitude|PID) + (Valence-1|PID)');
lme_model7{i}= fitlme(tbl_lme, 'RT~ 1+Magnitude*Valence + (1|PID)+(Magnitude-1|PID)+(Valence-1|PID)+(Valence:Magnitude-1|PID)');
% lme_model8{i}= fitlme(tbl_lme, 'RT~ 1+Magnitude+Magnitude:Valence + (Valence*Magnitude|PID)');
% lme_model9{i}= fitlme(tbl_lme, 'RT~ 1+Magnitude+Magnitude:Valence + (1|PID)+(Magnitude-1|PID)+(Valence-1|PID)+(Valence:Magnitude-1|PID)');
% lme_model10{i}= fitlme(tbl_lme,'RT~ 1+Magnitude+Magnitude:Valence + (1|PID)+(Magnitude-1|PID)+(Valence-1|PID)');
% lme_model11{i}= fitlme(tbl_lme,'RT~ 1+Magnitude+Magnitude:Valence + (Magnitude|PID)+(Valence-1|PID)');
% 

end


RTtot = [RTpar(:,1); RTpar(:,2)];
valtot= [valpar(:,1);valpar(:,2)];
PIDtot= [PID(:,1);PID(:,2)];
Valencetot= [Valence(:,1);Valence(:,2)];
Magnitudetot= [Magnitude(:,1);Magnitude(:,2)];
Time=[zeros(size(Magnitude(:,1))); ones(size(Magnitude(:,2)))];
[nor,lambda]=boxcox(RTtot);


tbl_lme= table(RTtot , Magnitudetot, Valencetot , PIDtot, Time, 'VariableNames', {'RT', 'Magnitude', 'Valence' , 'PID', 'Time' } );

% here pooling together the befroe and after as no test showed difference
% maybe centering intercept to lowest measured RT, but still not very
% meaningful

lme_model20{12}= fitlme(tbl_lme, 'RT~ 1+Magnitude+Valence + (1|PID)');
lme_model20{1}= fitlme(tbl_lme, 'RT~ 1+Magnitude+Valence + (Valence+Magnitude|PID)');
lme_model20{2}= fitlme(tbl_lme, 'RT~ 1+Magnitude+Valence + (1|PID)+(Magnitude-1|PID)+(Valence-1|PID)');
lme_model20{3}= fitlme(tbl_lme, 'RT~ 1+Magnitude+Valence + (Magnitude|PID) + (Valence-1|PID)');
lme_model20{4}= fitlme(tbl_lme, 'RT~ 1+Magnitude*Valence + (Valence*Magnitude|PID)');
lme_model20{5}= fitlme(tbl_lme, 'RT~ 1+Magnitude*Valence + (1|PID)+(Magnitude-1|PID)+(Valence-1|PID)');
lme_model20{6}= fitlme(tbl_lme, 'RT~ 1+Magnitude*Valence + (Magnitude|PID) + (Valence-1|PID)');
lme_model20{7}= fitlme(tbl_lme, 'RT~ 1+Magnitude*Valence + (1|PID)+(Magnitude-1|PID)+(Valence-1|PID)+(Valence:Magnitude-1|PID)');
% lme_model20{8}= fitlme(tbl_lme, 'RT~ 1+Magnitude+Magnitude:Valence + (Valence*Magnitude|PID)');
% lme_model20{9}= fitlme(tbl_lme, 'RT~ 1+Magnitude+Magnitude:Valence + (1|PID)+(Magnitude-1|PID)+(Valence-1|PID)+(Valence:Magnitude-1|PID)');
% lme_model20{10}= fitlme(tbl_lme,'RT~ 1+Magnitude+Magnitude:Valence + (1|PID)+(Magnitude-1|PID)+(Valence-1|PID)');
% lme_model20{11}= fitlme(tbl_lme,'RT~ 1+Magnitude+Magnitude:Valence + (Magnitude|PID)+(Valence-1|PID)');



lme_model30= fitglme(tbl_lme, 'RT~ 1+Magnitude+Valence + (1|PID)', 'Distribution', 'Gamma', 'Link', 'reciprocal');
figure;
plotResiduals(lme_model30,'fitted','ResidualType','Pearson');
figure;
plotResiduals(lme_model30,'histogram','ResidualType','Pearson')
lme_model31= fitglme(tbl_lme, 'RT~ 1+Magnitude+Valence + (1|PID)', 'Distribution', 'Normal', 'Link', 'reciprocal');
figure;
plotResiduals(lme_model31,'fitted','ResidualType','Pearson');
figure;
plotResiduals(lme_model31,'histogram','ResidualType','Pearson')
% % Plot the Pearson residuals versus lagged residuals, to check for correlation among the residuals. The conditional independence assumption in GLME implies that the conditional Pearson residuals are approximately uncorrelated




% Now all data but separating the before and after
% 
% lme_model21{1}= fitlme(tbl_lme, 'RT~ 1+Magnitude+Valence + Time +(Magnitude|PID:Time)+(Valence-1|PID:Time) ');
% lme_model21{2}= fitlme(tbl_lme, 'RT~ 1+Magnitude+Valence + Time +(Magnitude+Valence|PID:Time)');
% lme_model21{3}= fitlme(tbl_lme, 'RT~ 1+Magnitude+Valence + Time +(1|PID:Time)+(Magnitude-1|PID:Time)+(Valence-1|PID:Time)');
% lme_model21{4}= fitlme(tbl_lme, 'RT~ 1+Magnitude*Valence + Time +(Magnitude*Valence|PID:Time)');
% lme_model21{1}= fitlme(tbl_lme, 'RT~ 1+Magnitude*Valence + Time +(Magnitude|PID:Time)+(Valence-1|PID:Time) ');
% lme_model21{2}= fitlme(tbl_lme, 'RT~ 1+Magnitude*Valence + Time +(Magnitude+Valence|PID:Time)');
% lme_model21{3}= fitlme(tbl_lme, 'RT~ 1+Magnitude*Valence + Time +(1|PID:Time)+(Magnitude-1|PID:Time)+(Valence-1|PID:Time)');
% lme_model21{6}= fitlme(tbl_lme, 'RT~ 1+Magnitude + Time +(1|PID:Time)+(Magnitude-1|PID:Time)');
% lme_model21{7}= fitlme(tbl_lme, 'RT~ 1+Magnitude + Time +(Magnitude|PID:Time)');
% lme_model21{8}= fitlme(tbl_lme, 'RT~ 1+Magnitude +(Magnitude|PID:Time)');
% lme_model21{8}= fitlme(tbl_lme, 'RT~ 1+Magnitude +(Magnitude|PID/Time)');
% lme_model21{9}= fitlme(tbl_lme, 'RT~ 1+Magnitude +(1|PID:Time)+(Magnitude-1|PID:Time)');


% in the paper about theory of reactuon times they encourage to use random
% intercepts in RT analysis!!! so use this as reference. Moreover, with 
% age and diseases the RT change and increase. 


% idea, make a inverse-gamma lienar model , subject wise ,sot compare
% slopes, and the n chekc if residuals are ACF or PACF.



%% plotting RT to test if inverse gaussian fits better than normal

rt=RTtot;
x  = 0:50:1000; % reaction time bins (ms)
N  = hist(rt,x);  % count reactions in bin
N  = 100*N./sum(N); % percentage
 
figure(1) % graphics
h  = bar(x,N); % histogram/bar plot
% Making figure nicer to look at
set(h,'FaceColor',[.7 .7 .7]); % setting the color of the bars to gray
box off; % remove the "box" around the graph
axis square;
xlabel('Reaction time (ms)'); % and setting labels is essential for anyone to understand graphs!
ylabel('Probability (%)');
xlim([0 3000]);



% Inverse reaction time
rtinv  = 1./rt; % inverse reaction time / promptness (ms-1)
 
n    = numel(x); % number of bins in reaction time plot
x    = linspace(1/2000,0.01,n); % promptness bins
N    = hist(rtinv,x);
N    = 100*N./sum(N);
 
figure(2)
h = bar(x*1000,N);
hold on
set(h,'FaceColor',[.7 .7 .7]);
box off
axis square;
xlabel('Promptness (s^{-1})');
ylabel('Probability (%)');
title('Reciprocal time axis');
set(gca,'YTick',0:5:100,'XTick',0:1:8);
xlim([0 8]);
 
% Does this look like a Gaussian?
% Let's plot a Gaussian curve with mean and standard deviation from the
% promptness data
mu  = mean(rtinv);
sd  = std(rtinv);
y  = normpdf(x,mu,sd);
y  = y./sum(y)*100;
plot(x*1000,y,'ks-','LineWidth',2,'MarkerFaceColor','w');

[h,p]  = kstest(zscore(rt)); % for reaction time
if h
  str    = ['The null hypothesis that the reaction time distribution is Gaussian distributed is rejected, with P = ' num2str(p)];
else
  str    = ['The null hypothesis that the reaction time distribution is Gaussian distributed is NOT rejected, with P = ' num2str(p)];
end
disp(str); % display the results in command window
[h,p]  = kstest(zscore(rtinv)); % for inverse reaction time
if h
  str    = ['The null hypothesis that the inverse reaction time distribution is Gaussian distributed is rejected, with P = ' num2str(p)];
else
  str    = ['The null hypothesis that the inverse reaction time distribution is Gaussian distributed is NOT rejected, with P = ' num2str(p)];
end
disp(str); % display the results in command window

% Cumulative probability
% Normalized scale and nicer shape
x = sort(1000*rtinv);
n = numel(rtinv); % number of data points
y = 100*(1:n)./n; % cumulative probability for every data point
figure(3)
plot(x,y,'k.');
hold on
 
% Now, plot it cumulative probability in quantiles
% this is easier to compare between different distributions
p    = [1 2 5 10:10:90 95 98 99]/100; % some arbitrary probabilities
q    = quantile(rtinv,p); % instead of hist, let's use quantiles
 
h = plot(q*1000,p*100,'ko','LineWidth',2,'MarkerFaceColor','r');
hold on
xlabel('Promptness (s^{-1})');
ylabel('Cumulative probability (%)');
title('Cumulative probability plot');
box off
axis square;
set(gca,'YTick',0:10:100,'XTick',0:1:8);
xlim([0 8]);
legend(h,'Quantiles','Location','SE');
 
% % optional, save figure
% print('-depsc','-painter',[mfilename '3']);

cdf     = q;
myerf       = 2*cdf - 1;
myerfinv    = sqrt(2)*erfinv(myerf);
chi         = myerfinv;

 % Probit
figure(4)
% raw data
x = -1./sort((rt)); % multiply by -1 to mirror abscissa
n = numel(rtinv); % number of data points
% y = pa_probit((1:n)./n); % cumulative probability for every data point converted to probit scale
% plot(x,y,'k.');
% hold on
%  
% % quantiles
% p    = [1 2 5 10:10:90 95 98 99]/100;
% probit  = pa_probit(p);
% q    = quantile(rt,p);
% q    = -1./q;
% xtick  = sort(-1./(150+[0 pa_oct2bw(50,-1:5)])); % some arbitrary xticks
%  
% plot(q,probit,'ko','Color','k','MarkerFaceColor','r','LineWidth',2);
% hold on
% set(gca,'XTick',xtick,'XTickLabel',-1./xtick);
% xlim([min(xtick) max(xtick)]);
% set(gca,'YTick',probit,'YTickLabel',p*100);
% ylim([pa_probit(0.1/100) pa_probit(99.9/100)]);
% axis square;
% box off
% xlabel('Reaction time (ms)');
% ylabel('Cumulative probability');
% title('Probit ordinate');
%  
% % this should be a straight line
% x = q;
% y = probit;
% b = regstats(y,x);
% h = pa_regline(b.beta,'k-');
% set(h,'Color','r','LineWidth',2);
%  
%  
% % optional, save figure
% print('-depsc','-painter',[mfilename '4']);






%% RManova - 2 way
% h=1;
% for kk=[1:7,9:11]
%     
%    RTrm(h,1:60)= RT{kk}(:,1)';
%    h=h+1;  
% end
% 
% tab=table([1:7,9:11]', 'VariableNames', {'PID' });
% 
% %  h=1;
% for hh=1:60
%     var_name=['RT' num2str(hh)]
%     tbl= table( RTrm(:,hh), 'VariableNames', {var_name});
%     tab=[tab, tbl]
% %     h=h+1;
% end

% 
% level_mag=unique(Magnitude);
% level_val=unique(Valence);
% 
% for i=[6]
%    RTrm=RT{i}(:,1);
%    valrm=val{i}(:,1);
%    Valence=val{i}(:,1)>0; % pos or neg
%    Magnitude=abs(val{i}(:,1));
%    PID= i;
%    
%    indValPosMag25= find((Magnitude==25) & Valence)
%    indValPosMag50= find((Magnitude==50) & Valence)
%    indValPosMag75= find((Magnitude==75) & Valence)
%    indValPosMag100= find((Magnitude==100) & Valence)
%    indValNegMag25= find((Magnitude==25) & (Valence==0))
%    indValNegMag50= find((Magnitude==50) & (Valence==0))
%    indValNegMag75= find((Magnitude==75) & (Valence==0))
%    indValNegMag100= find((Magnitude==100) & (Valence==0))
%    
%    hh=[0,0,0,0,1,1,1,1];
%    h=1;
%    for jj=[25:25:100, 25:25:100]
%        eval([' ind= find((Magnitude==' num2str(jj) ' ) & (Valence==' num2str(hh(h)) ';))' ])
%    
%        RTrmfinal(i,h)=nanmean(RTrm(ind));
%        h=h+1;
%    end
% 
% end
% 
% h=1;
% for kk=[1:7,9:11]
%     
%    RTrm(h,1:60)= RT{kk}(:,1)';
%    h=h+1;  
% end
% 
% tab=table([1:7,9:11]', 'VariableNames', {'PID' });



%% 2 way anova repeated measures Fixed effects though

% for rating response time, and 2 factors, saliency and valence
% maybe too unbalanced for matlab... possible with R


% Haynes paper on salience and value By using both appetitive and aversive cues, this task allows for
% linearly independent (i.e., uncorrelated) levels of value and salience
% (Fig. 1), they do multiple regression on Rating trials of salince and
% value and test if betas are positive individually, then for group just
% random effect with t-test. 
% for comparisson of value vs salience, as orthogonal, if both positive,
% then rmANOVA with factors value and valence, should be an interaction.
% 
% tbl= table(RTtot, PID, Magnitude, Valence, 'VariableNames',{'RTtot','PID','Magnitude','Valence'});
% fitrm(tbl)


% mag=unique(Magnitude); valen=unique(Valence);
% h=1;
% for i=1:length(mag)
%    
%     ind1= find(Magnitude==mag(i));
%     for j=1:length(valen)
%         ind2=find(Valence==valen(i));
%         ind3{h}=intersect(ind1, ind2);
%         h=h+1;
%     end
% end
% 
%    
%     
%   
%    RTtotal=RTtot(ind);
%    valtot(end+1:end+length(RT{1}))=val{i}(:,1);
%    Magnitude=abs(valtot);                          %(25, 50, 75, 100)
%    Valence(end+1:end+length(RT{1}))=val{i}(:,1)>0; % pos or neg
%    PID(end+1:end+length(RT{1}))= i;
% 
% 
% 
%     ind1= find(Magnitude==mag(i));
%     for j=1:length(valen)
%         ind2=find(Valence==valen(i));
%         ind3{h}=intersect(ind1, ind2);
%         h=h+1;
%     end









% Random effect analysis of second level analysis. Here, I analyse the ZPF
% of all patients, and test if bigger thanzero. Due to CLT I suppose they
% follow normal distribution. This is similar to fmri second level with 
% random effects. This is more restrictive than fixed effects, but it
% allows for exptrapolation. 


% http://saxelab.mit.edu/resources/papers/in_press/federetal.pdf

% I use random effects t oavoid “fixed-effects fallacy”(Clark, 1973), r
% the unfounded inference that conclusions about items sampled from a
% population generalize to the entire item population. 

% %. In treating  subjects as a “ fixed ” variable, the variability of an effect across subjects
% is not taken into account, which is especially problematic if there is
% substantial variability in the effect size across participants. Consider
% an instance in which only one out of
% fi
% ve subjects shows a very large
% effect for some condition, and the other subjects show no effect.
% Averaging the effect across subjects, the entire group will appear to
% exhibit a medium-sized effect
%


