% function [H,P,CI, stats]=response_time_multilevel

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
   
        [rho2AFCabs{i}, rhoRat{i}, pval2AFCabs{i}, pvalRat{i}, ZPF{i}, pvalZPF{i}, ZPFRat{i}, pvalZPFrat{i},lm_2afc{i}, lm_rat{i}, lm_interact{i}, lm_MagValen{i}, tbetas{i}, pbetas{i}, RT{i}, val{i},  pottlm{i}, compFval{i},interPval{i}] = response_time_ana(f{i},type);
  
end


h=1; 
for i=[1:7, 9:11]     % skip patient 39 because he did not eat food, not used in any test of this type

    % 2AFC correlations
    
    for j=1:length(ZPF{i})
       eval(['t1' num2str(j) '(h) = ZPF{i}(j);'])
       eval(['t2' num2str(j) '(h) = rho2AFCabs{i}(j,1);'])
       eval(['t3' num2str(j) '(h) = rho2AFCabs{i}(j,2);'])
    end
    
    
    
    for j=1:length(tbetas{i})
        
     % 2AFC Regression
       eval(['t4' num2str(j) '(h) = tbetas{i}(j) ;'])
    end
    
    for j=1:2
     
     % 2AFC Regression
       eval(['t5' num2str(j) '(h) = lm_2afc{i}{j}.Coefficients.tStat(2);'])
       eval(['t6' num2str(j) '(h) = lm_2afc{i}{j}.Coefficients.Estimate(2);'])
     
     % Rating Regression  
       eval(['t7' num2str(j) '(h) = lm_rat{i}{j}.Coefficients.tStat(2);'])
       eval(['t8' num2str(j) '(h) = lm_rat{i}{j}.Coefficients.tStat(3);'])
       eval(['t9' num2str(j) '(h) = lm_rat{i}{j}.Coefficients.Estimate(2);'])
       eval(['t10' num2str(j) '(h) = lm_rat{i}{j}.Coefficients.Estimate(3);'])
       eval(['t14' num2str(j) '(h) = lm_interact{i}{j}.Coefficients.Estimate(4);'])
       eval(['t15' num2str(j) '(h) = lm_interact{i}{j}.Coefficients.tStat(4);']) 
     % Magnitude and Valence Regression  
       eval(['t16' num2str(j) '(h) = lm_MagValen{i}{j}.Coefficients.Estimate(2);'])
       eval(['t17' num2str(j) '(h) = lm_MagValen{i}{j}.Coefficients.tStat(2);']) 
       eval(['t18' num2str(j) '(h) = lm_MagValen{i}{j}.Coefficients.Estimate(3);'])
       eval(['t19' num2str(j) '(h) = lm_MagValen{i}{j}.Coefficients.tStat(3);']) 
       eval(['t20' num2str(j) '(h) = lm_MagValen{i}{j+2}.Coefficients.Estimate(4);'])
       eval(['t21' num2str(j) '(h) = lm_MagValen{i}{j+2}.Coefficients.tStat(4);']) 
     % Rating correlations, 1 difference then each correlation
      
       eval(['t11' num2str(j) '(h) = ZPFRat{i}(j);'])
       eval(['t12' num2str(j) '(h) = rhoRat{i}(j,1);'])
       eval(['t13' num2str(j) '(h) = rhoRat{i}(j,2);'])
       
%        eval(['t16' num2str(j) '(h) = lm_sum{i}{j}.Coefficients.Estimate(2);'])
%        eval(['t17' num2str(j) '(h) = lm_sum{i}{j}.Coefficients.tStat(2);'])
%      sum only in the model if difference not included, moreoer, quite
%      collinear in collintest, and correlted 0.6, way way more than I
%      expected
    end
    h=h+1;
end


% 2AFC Correlation tests

for i=1:length(ZPF{1}) 
    if eval(['kstest(t1' num2str(i) ') ==0'])     
        eval(['[H_2corrAB{i},P_2corrAB{i},CI_2corrAB{i}, stats_2corrAB{i}]=ttest(t1' num2str(i) ');']);
    end

       
        eval(['[P_2corrB{i},H_2corrB{i}] =signrank(t2' num2str(i) ');']);
        eval(['[P_2corrA{i},H_2corrA{i}] =signrank(t3' num2str(i) ');']);
        eval(['[P_2corrAB2{i},H_2corrAB2{i}] =signrank(t1' num2str(i) ');']);
end


% Rating correlation tests

for i=1:2
   
    eval(['[P_RcorrAB{i}, H_RcorrAB{i}]=signrank(t11' num2str(i) ');']);
    eval(['[P_RcorrB{i}, H_RcorrB{i}]  =signrank(t12' num2str(i) ');']);
    eval(['[P_RcorrA{i}, H_RcorrA{i}]  =signrank(t13' num2str(i) ');']);
end



% Regression tests

for i=1:length(tbetas{1})
    

    if eval(['kstest(t4' num2str(i) ') ==0'])      
       eval(['[H_betasAB{i},P_betasAB{i},CI_betasAB{i}, stats_betasAB{i}]=ttest(t4' num2str(i) ');']);
    end
end

for i=1:2    
    % Difference value in 2AFC
    eval(['[P_2tstat{i}, H_2tstat{i}]=signrank(t5' num2str(i) ');']);
    eval(['[P_2beta{i}, H_2beta{i}]=signrank(t6' num2str(i) ');']);
    
    % Value and Salience in ratinf
    eval(['[P_RtstatVal{i}, H_RtstatVal{i}]=signrank(t7' num2str(i) ');']);
    eval(['[P_RtstatSal{i}, H_RtstatSal{i}]=signrank(t8' num2str(i) ');']);
    eval(['[P_RbetaVal{i}, H_RbetaVal{i}]=signrank(t9' num2str(i) ');']);
    eval(['[P_RbetaSal{i}, H_RbetaSal{i}]=signrank(t10' num2str(i) ');']);
    eval(['[P_RbetaInteraction{i}, H_RbetaInteraction{i}]=signrank(t14' num2str(i) ');']);
    eval(['[P_RtstatInteraction{i}, H_RtstatInteraction{i}]=signrank(t15' num2str(i) ');']);
    
    % Magnitude and valence in rating
    eval(['[P_betaMag{i}, H_betaMag{i}]=signrank(t16' num2str(i) ');']);
    eval(['[P_tstatMag{i}, H_tstatMag{i}]=signrank(t17' num2str(i) ');']);
    eval(['[P_betaValence{i}, H_betaValence{i}]=signrank(t18' num2str(i) ');']);
    eval(['[P_tstatValence{i}, H_tstatValence{i}]=signrank(t19' num2str(i) ');']);
    eval(['[P_betaMagValInter{i}, H_betaMagValInter{i}]=signrank(t20' num2str(i) ');']);
    eval(['[P_tstatMagValInter{i}, H_tstatMagValInter{i}]=signrank(t21' num2str(i) ');']);
%     eval(['[P_betaSum{i}, H_betaSum{i}]=signrank(t16' num2str(i) ');']);
%     eval(['[P_tstatSum{i}, H_tstatSum{i}]=signrank(t17' num2str(i) ');']);
end



%% ANOVA simple  This is completely Fixed effects though

RTpar=[];
valpar=[];
Valence=[];
PID=[];
Magnitude=[];

for i=1:length(val)
    
   RTpar(end+1:end+length(RT{1}),1)=RT{i}(:,1);
   valpar(end+1:end+length(RT{1}),1)=val{i}(:,1);

   Valence(end+1:end+length(RT{1}),1)=val{i}(:,1)>0; % pos or neg
   PID(end+1:end+length(RT{1}),1)= i;
   
   RTpar(end-length(RT{1})+1:end,2)=RT{i}(:,2);
   valpar(end-length(RT{1})+1:end,2)=val{i}(:,2);
                           %(25, 50, 75, 100)
   Valence(end-length(RT{1})+1:end,2)=val{i}(:,2)>0; % pos or neg
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
%give answeres for each of the 8 possible conditins, so as missing data

for i=1:2
tbl_lme= table(RTpar(:,i) , Magnitude(:,i) , Valence(:,i) , PID(:,i), 'VariableNames', {'RT', 'Magnitude', 'Valence' , 'PID'} );
lme_model{i}= fitlme(tbl_lme,  'RT~ 1+Magnitude+Valence + (Valence+Magnitude|PID)');
lme_model2{i}= fitlme(tbl_lme, 'RT~ 1+Magnitude+Valence + (1|PID)+(Magnitude+Valence-1|PID)');
lme_model3{i}= fitlme(tbl_lme, 'RT~ 1+Magnitude*Valence + (Valence+Magnitude|PID)');
lme_model4{i}= fitlme(tbl_lme, 'RT~ 1+Magnitude*Valence + (1|PID)+(Magnitude+Valence-1|PID)');
lme_model5{i}= fitlme(tbl_lme, 'RT~ 1+Magnitude*Valence + (Valence*Magnitude|PID)');
lme_model6{i}= fitlme(tbl_lme, 'RT~ 1+Magnitude*Valence + (1|PID)+(Magnitude*Valence-1|PID)');
lme_model7{i}= fitlme(tbl_lme, 'RT~ Magnitude*Valence + (1|PID)+(Magnitude*Valence-1|PID)');
end


RTtot = [RTpar(:,1); RTpar(:,2)];
valtot= [valpar(:,1);valpar(:,2)];
PIDtot= [PID(:,1);PID(:,2)];
Valencetot= [Valence(:,1);Valence(:,2)];
Magnitudetot= [Magnitude(:,1);Magnitude(:,2)];
Time=[zeros(size(Magnitude(:,1))); ones(size(Magnitude(:,2)))];

tbl_lme= table(RTtot , Magnitudetot, Valencetot , PIDtot, Time, 'VariableNames', {'RT', 'Magnitude', 'Valence' , 'PID', 'Time' } );
lme_model10{1}= fitlme(tbl_lme, 'RT~ 1+Magnitude+Valence + Time +(1|PID:Time)+(Magnitude-1|PID:Time)+(Valence-1|PID:Time) ');
lme_model10{2}= fitlme(tbl_lme, 'RT~ 1+Magnitude+Valence + Time +(Magnitude+Valence|PID:Time)');
lme_model10{3}= fitlme(tbl_lme, 'RT~ 1+Magnitude+Valence + Time +(1|PID:Time)+(Magnitude+Valence-1|PID:Time)');
lme_model10{4}= fitlme(tbl_lme, 'RT~ 1+Magnitude*Valence + Time +(1|PID:Time)+(Magnitude*Valence-1|PID:Time)');




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


