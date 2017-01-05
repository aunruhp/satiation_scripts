function [rho2AFCabs, rhoRat, pval2AFCabs, pvalRat, ZPF, pvalZPF, ZPFRat, pvalZPFrat,lm_2afc, lm_rat, lm_ratInter, lm_MagValen, lm_salienceonly, tbetas, pbetas, response_time_rating ,val_rat , pottlm, compFval, interPval,response_time_2AFC, val_2AFC_abs_out ] = response_time_ana (DirName, type)

% In this analysis I test if the higher the value difference between the 2
% products, the faster the decision, because according to the
% drift-diffusion model it should be so as the threshold is reached earlier.
% With this I want to analyse if the INITIATION of the response is faster,
% of course, then a specific likert values has to be chosen, and the higher
% the difference the more button presses needed, thus only the initiation is
% checked. It should be negatively correlated to absolute difference.
%
%
%
% Output: 
% rho... are the values from Spearman correlation for Response Times and
% different values. 
% pval.. are the pvals dor these correlations
% a lot of stuff joe!



dbstop if error


%% Load the data
if nargin ==0
    [satcue, ranking, events,~,DirName] = load_behavioural_data;
    
elseif nargin>0
     [satcue, ranking, events] = load_behavioural_data(DirName);
end

if nargin <2
    type=0;  % type 0 means Pearsons matrix of correlation, if 1 means Spearmans
end


% check if any event is a time out, so event == 31

if any(events(:,2)==31)
    warning(['patient' DirName 'has' sum(events(:,2)==31) 'time outs, check it and avoid them for your analysis as outliers'])
end



%% Extract the real reaction times

[RT] = extract_reaction_times(events);

% satcue{3}(:,13)  was used at the beginning of this analysis, 
% saved in runsat_output, but those where the end of the
% response. Thus, not really what I wanted
% Now I calculate directly from the events the distance between the
% likert onset and the first response
   



%% Start analysis, correlations of RT with values

h=1;
for jj=[3,6]  % 3 is before pause, and 6 after it
    
    
    %% RT to use
    response_time_2AFC(:,h)= RT{jj};             % RT for 2AFC vector, 190 entries
    response_time_rating(:,h)= RT{jj-1};         % RT for rating vector, 60 entries
    
    
    
    %% Choice trials analysis (2AFC)
    
    
    
    val_A_rat=nan(length(satcue{jj}(:,14)),1);             % value of image on A given in previous rating trials (-300 to 300)
    val_B_rat=nan(length(satcue{jj}(:,14)),1);             % value of image on B given in previous rating trials (-300 to 300)
    val_left_rat=  nan(length(satcue{jj}),1);              % value of image on the left given in previous rating trials (-300 to 300)
    val_right_rat= nan(length(satcue{jj}),1);              % value of image on the right given in previous rating trials (-300 to 300)
    val_A_2AFC=nan(length(satcue{jj}(:,14)),1);            % value of image on A given in previous whole 2AFC trials (0 to 1900)
    val_B_2AFC=nan(length(satcue{jj}(:,14)),1);            % value of image on B given in previous whole 2AFC trials (0 to 1900)
    val_A_ranking=nan(length(satcue{jj}(:,14)),1);         % rank of image on A, calculated after all the choice trials (0 to 20) 
    val_B_ranking=nan(length(satcue{jj}(:,14)),1);         % rank of image on B, calculated after all the choice trials (0 to 20) 
    val_A_2AFCtotal=nan(length(satcue{jj}(:,14)),1);            % value of image on A given in previous whole 2AFC trials (0 to 1900)
    val_B_2AFCtotal=nan(length(satcue{jj}(:,14)),1);            % value of image on B given in previous whole 2AFC trials (0 to 1900)

    for i=1:20
        
        ind=find(cell2mat(satcue{jj}(:,2))==i);            % when this stimulus was on left
        ind2=find(cell2mat(satcue{jj}(:,3))==i);           % when this stimulus was on right
        lr_to_ab=cell2mat(satcue{jj}(ind,8));              % conversion factor from left_rigth to A_B
        lr_to_ab2=cell2mat(satcue{jj}(ind2,8));            % conversion factor from left_rigth to A_B
        indA = cat(1,((lr_to_ab) .* ind), (~lr_to_ab2).* ind2);  indA(indA==0)= [];  %indices for A
        indB = cat(1,((~lr_to_ab) .* ind), (lr_to_ab2).* ind2);  indB(indB==0)= [];  %indices for B
        
        
        valRat= ranking{jj-1}(i,4);                        % rating value (-300 to 300)  
        valRank= ranking{jj}(i,3);                         % ranking (0 to 20)
        val2AFC= ranking{jj}(i,4);                         % choice rating (0-1900)
        val2AFCtotal= ranking{jj}(i,5);                         % choice rating (-1900-1900)
        
                
        val_A_rat(indA)= valRat;
        val_B_rat(indB)= valRat;
%         val_left_rat(ind)= valRat;
%         val_right_rat(ind2)= valRat;
        val_A_ranking(indA)= valRank;
        val_B_ranking(indB)= valRank;
        val_A_2AFC(indA)= val2AFC;
        val_B_2AFC(indB)= val2AFC;
        val_A_2AFCtotal(indA)= val2AFCtotal;
        val_B_2AFCtotal(indB)= val2AFCtotal;
       
        
    end

     
    
    %% Correlations of reponse time and value difference and absolute difference
    
    % Output correlation 
    
    % likert 2AFC of EACH TRIAL (-100 to 100) and response_time for choice,
    % TRIAL BY TRIAL ANALYSIS betwwen output and RT
    val_2AFC_out(:,h)= cell2mat(satcue{jj}(:,14));
    [rho2AFC(1,h), pval2AFC(1,h)]=corr(response_time_2AFC(:,h), val_2AFC_out(:,h), 'type', 'Spearman');  %any correlation between
    
    
    % absolute value of likert 2AFC of EACH TRIAL (-100 to 100) and response_time choice trial by trial
    % TRIAL BY TRIAL ANALYSIS betwwen output and RT  
    val_2AFC_abs_out(:,h)= abs(val_2AFC_out(:,h));
    [rho2AFCabs(1,h), pval2AFCabs(1,h)]=corr(response_time_2AFC(:,h), val_2AFC_abs_out(:,h), 'type', 'Spearman');
    
    
%     
%     % Input correlation 
%     
%     %form here on value of all ratings, or of all rankings, not the real
%     %output of the decision
%     
%     % difference of input likert rating A-B (-300-300) and response_time for choice (0-300)
%     val_diff_rat(:,h)= val_A_rat - val_B_rat;
%     [rho2AFC(2,h), pval2AFC(2,h)]=corr(response_time_2AFC(:,h), val_diff_rat(:,h), 'type', 'Spearman');
%     
    % absolute difference of likert rating (A-B) (0-600) and response_time for choice
    val_abs_rat(:,h)= abs(val_A_rat - val_B_rat);
    [rho2AFCabs(2,h), pval2AFCabs(2,h)]=corr(response_time_2AFC(:,h), val_abs_rat(:,h), 'type', 'Spearman');
    [rho2AFCabs(3,h), pval2AFCabs(3,h)]=corr(val_2AFC_abs_out(:,h), val_abs_rat(:,h), 'type', 'Spearman');
%     
%     % difference of likert 2AFC_rating and response_time for choice  (0-1900)
%     val_diff_2AFC(:,h)= val_A_2AFC - val_B_2AFC;
%     [rho2AFC(3,h), pval2AFC(3,h)]=corr(response_time_2AFC(:,h), val_diff_2AFC(:,h), 'type', 'Spearman');
%     
%     % absolute difference of likert 2AFC_rating and response_time for choice (0-1900)
%     val_abs_2AFC(:,h)= abs(val_diff_2AFC(:,h));
%     [rho2AFCabs(3,h), pval2AFCabs(3,h)]=corr(response_time_2AFC(:,h), val_abs_2AFC(:,h), 'type', 'Spearman');
%     
%     
%     % difference of likert ranking and response_time for choice
%     val_diff_rank(:,h)= val_A_ranking - val_B_ranking;
%     [rho2AFC(4,h), pval2AFC(4,h)]=corr(response_time_2AFC(:,h), val_diff_rank(:,h), 'type', 'Spearman');
%     
%     % absolute difference of ranking rating and response_time for choice
%     val_abs_rank(:,h)= abs(val_diff_rank(:,h));
%     [rho2AFCabs(4,h), pval2AFCabs(4,h)]=corr(response_time_2AFC(:,h), val_abs_rank(:,h), 'type', 'Spearman');
%     
%    % difference of likert 2AFC_rating and response_time for choice  (-1900-1900)
%     val_diff_total(:,h)= val_A_2AFCtotal - val_B_2AFCtotal;
%     [rho2AFC(5,h), pval2AFC(5,h)]=corr(response_time_2AFC(:,h), val_diff_total(:,h), 'type', 'Spearman');
%     
%      % absolute difference of likert 2AFC_rating and response_time for choice (-1900-1900)
%     val_abs_total(:,h)= abs(val_diff_total(:,h));
%     [rho2AFCabs(5,h), pval2AFCabs(5,h)]=corr(response_time_2AFC(:,h), val_abs_total(:,h), 'type', 'Spearman');
%     
%     
% %     
% %     % Freaking out bro
%     
% %     % difference of likert rating (left - right) and response_time for choice
% %     val_diff_lr= val_left_rat - val_right_rat;
% %     [rho2AFC(6,h), pval2AFC(6,h)]=corr(response_time_2AFC, val_diff_lr, 'type', 'Spearman');
% %     % difference of likert rating (left - right) and response_time for choice
% %     val_diff_lr_abs= abs(val_diff_lr);
% %     [rho2AFCabs(6,h), pval2AFCabs(6,h)]=corr(response_time_2AFC, val_diff_lr_abs, 'type', 'Spearman');
%     
%     
%     % this analysis correlates the sum of the values with the RT
%     % difference of input likert rating A-B (-300-300) and response_time for choice (0-300)
%     val_sum_rat(:,h)= val_A_rat + val_B_rat;
%     [rho2AFC(7,h), pval2AFC(7,h)]=corr(response_time_2AFC(:,h), val_sum_rat(:,h), 'type', 'Spearman');
%     
%     % this analysis correlates the absolute sum of the values with the RT
%     % difference of input likert rating A-B (-300-300) and response_time for choice (0-300)
    val_sum_rat_abs(:,h)= abs(val_A_rat + val_B_rat);
    [rho2AFCabs(7,h), pval2AFCabs(7,h)]=corr(response_time_2AFC(:,h), val_sum_rat_abs(:,h), 'type', 'Spearman');
    
%     % this analysis correlates the absolute sum of the saleincy, thus the absolute sum of the distance to zero with the RT
%     % difference of input likert rating A-B (-300-300) and response_time for choice (0-300)
%     val_sum_rat_sal(:,h)= abs(val_A_rat) + abs(val_B_rat);
%     [rho2AFCabs(8,h), pval2AFCabs(8,h)]=corr(response_time_2AFC(:,h), val_sum_rat_sal(:,h), 'type', 'Spearman');
%     
%      
%     % this analysis correlates the absolute difference of the saleincy, thus the absolute sum of the distance to zero with the RT
    % difference of input likert rating A-B (-300-300) and response_time for choice (0-300)
%     val_sum_rat_sal_abs(:,h)= abs(val_A_rat) - abs(val_B_rat);
%     [rho2AFCabs(9,h), pval2AFCabs(9,h)]=corr(response_time_2AFC(:,h), val_sum_rat_sal_abs(:,h), 'type', 'Spearman');

%     %absolute of absolute difference saliency
%     val_sum_rat_sal_abs2(:,h)= abs(abs(val_A_rat) - abs(val_B_rat));
%     [rho2AFCabs(10,h), pval2AFCabs(10,h)]=corr(response_time_2AFC(:,h), val_sum_rat_sal_abs2(:,h), 'type', 'Spearman');

    
    

    
    
    
    %% Rating trials analysis
    
    val_rat_rat(:,h)=nan(length(satcue{jj-1}(:,14)),1);
    val_rat_2AFC(:,h)=nan(length(satcue{jj-1}(:,14)),1);
    val_rat_rank(:,h)=nan(length(satcue{jj-1}(:,14)),1);
    val_rat_total(:,h)=nan(length(satcue{jj-1}(:,14)),1);
    
    
    for i=1:20
        
        ind3=find(cell2mat(satcue{jj-1}(:,2))==i);
        valRat=  ranking{jj-1}(i,4); % ranking
        val2AFC= ranking{jj}(i,4);
        valRank= ranking{jj}(i,3);
        valTotal= ranking{jj}(i,5);
        val_rat_rat(ind3,h)= valRat;
        val_rat_2AFC(ind3,h)= val2AFC;
        val_rat_rank(ind3,h)= valRank;
        val_rat_total(ind3,h)= valTotal;
        
    end
    
    
%OUTPUT

  %VALUE
    % Evaluation likert rating and response_time for rating, the value of each trial is used (-100 to 100)
    val_rat(:,h)= cell2mat(satcue{jj-1}(:,14));
    [rhoRat(1,h), pvalRat(1,h)]=corr(response_time_rating(:,h), val_rat(:,h), 'type', 'Spearman');
  %SALIENCE  
    % Absolute likert rating and response_time for rating, the value of each trial is used (-100 to 100)
    val_rat_abs(:,h)= abs(cell2mat(satcue{jj-1}(:,14)));
    [rhoRat(2,h), pvalRat(2,h)]=corr(response_time_rating(:,h), val_rat_abs(:,h), 'type', 'Spearman');
%     
% % INPUT
%     %VALUE
%     % Evaluation likert rating and response_time for rating, the value of all the trials is used (-300 to 300)
%     [rhoRat(3,h), pvalRat(3,h)]=corr(response_time_rating(:,h), val_rat_rat(:,h), 'type', 'Spearman');
%     
%     % 2AFC Likert rating and response_time for rating, , the value of all the trials is used (0 to 1900)    
%     [rhoRat(4,h), pvalRat(4,h)]=corr(response_time_rating(:,h), val_rat_2AFC(:,h), 'type', 'Spearman');
%     
%     % Ranking and response_time for rating, the value of all the trials is used (0 to 20)
%     [rhoRat(5,h), pvalRat(5,h)]=corr(response_time_rating(:,h), val_rat_rank(:,h), 'type', 'Spearman');
%     
%      % Evaluation likert rating and response_time for rating, the value of all the trials is used (-1900 to 1900)
%     [rhoRat(6,h), pvalRat(6,h)]=corr(response_time_rating(:,h), val_rat_total(:,h), 'type', 'Spearman');
%     
%     %SALIENCE
%     % Evaluation likert rating and response_time for rating, the value of all the trials is used (-300 to 300)
%     [rhoRat(7,h), pvalRat(7,h)]=corr(response_time_rating(:,h), abs(val_rat_rat(:,h)), 'type', 'Spearman');
%     
%     % 2AFC Likert rating and response_time for rating, , the value of all the trials is used (0 to 1900)    
%     [rhoRat(8,h), pvalRat(8,h)]=corr(response_time_rating(:,h), abs(val_rat_2AFC(:,h)), 'type', 'Spearman');
%     
%     % Ranking and response_time for rating, the value of all the trials is used (0 to 20)
%     [rhoRat(9,h), pvalRat(9,h)]=corr(response_time_rating(:,h), abs(val_rat_rank(:,h)), 'type', 'Spearman');
%     
%     % 2AFC Likert rating and response_time for rating, , the value of all the trials is used (0 to 1900)    
%     [rhoRat(10,h), pvalRat(10,h)]=corr(response_time_rating(:,h), abs(val_rat_total(:,h)), 'type', 'Spearman');
%     
    %% Subject-wise regression of RT
  

    
%    %% to delete
%       
%     tbl1=table((response_time_2AFC(:,h)), (val_sum_rat_abs(:,h)), (val_2AFC_abs_out(:,h)),'VariableNames',{'RT','sum','diff'});
%     tbl2=table((response_time_rating(:,h)), (val_rat(:,h)), (val_rat_abs(:,h)),'VariableNames',{'RT','val','sal'});
%     
%     lm_rat{h} =fitglm(tbl2,'RT ~ val+ sal', 'link', 'reciprocal','Distribution', 'gamma'); %, 'Distribution', 'gamma'
%     lm_ratInter{h} =fitglm(tbl2,'RT ~ val* sal', 'link', 'reciprocal','Distribution', 'gamma');    
%     lm_2afc{h}=fitglm(tbl1,'RT ~ diff', 'link', 'reciprocal','Distribution', 'gamma');
% 
%     % Sum is not in the model when combined with difference, so not used
%       
% 
%     valin=val_rat(:,h);  
%     ind=find(valin>=0);
%     valence=zeros(size(valin));
%     valence(ind)=1;
%     magnitude= abs(valin);
%     tbl3=table((response_time_rating(:,h)), (magnitude), (valence),'VariableNames',{'RT','magnit','valen'});
%     lm_MagValen{h} =fitglm(tbl3,'RT ~ magnit + valen', 'link', 'reciprocal','Distribution', 'gamma');
%     lm_MagValen{h+2} =fitglm(tbl3,'RT ~ magnit * valen', 'link', 'reciprocal','Distribution', 'gamma');
% %     lm_Mag{h} =fitglm(tbl3,'RT ~ magnit', 'link', 'reciprocal','Distribution', 'gamma');
%     rt=response_time_rating(:,h);
%     rtinv  = 1./rt; 
%     [h1,p]  = kstest(zscore(rt)); % for reaction time
%     if p<0.05/10
%         str    = [DirName, ' Rating The null hypothesis that the reaction time distribution is Gaussian distributed is rejected, with P = ' num2str(p)];
%     else
%         str    = [DirName, ' Rating The null hypothesis that the reaction time distribution is Gaussian distributed is NOT rejected, with P = ' num2str(p)];
%     end
%     disp(str); % display the results in command window
%     [h2,p]  = kstest(zscore(rtinv)); % for inverse reaction time
%     if p<0.05/10
%         str    = [DirName, ' Rating The null hypothesis that the inverse reaction time distribution is Gaussian distributed is rejected, with P = ' num2str(p)];
%     else
%         str    = [DirName, ' Rating The null hypothesis that the inverse reaction time distribution is Gaussian distributed is NOT rejected, with P = ' num2str(p)];
%     end
%     disp(str); % display the results in command window
%     [~, lambda]=boxcox(rt)
%     
%     rt=response_time_2AFC(:,h);
%     rtinv  = 1./rt; 
%     [h1,p]  = kstest(zscore(rt)); % for reaction time
%     if p<0.05/10
%         str    = [DirName, ' 2AFC The null hypothesis that the reaction time distribution is Gaussian distributed is rejected, with P = ' num2str(p)];
%     else
%         str    = [DirName, ' 2AFC The null hypothesis that the reaction time distribution is Gaussian distributed is NOT rejected, with P = ' num2str(p)];
%     end
%     disp(str); % display the results in command window
%     [h2,p]  = kstest(zscore(rtinv)); % for inverse reaction time
%     if p<0.05/10
%         str    = [DirName, ' 2AFC The null hypothesis that the inverse reaction time distribution is Gaussian distributed is rejected, with P = ' num2str(p)];
%     else
%         str    = [DirName, ' 2AFC The null hypothesis that the inverse reaction time distribution is Gaussian distributed is NOT rejected, with P = ' num2str(p)];
%     end
%     disp(str); % display the results in command window
%      
%     
%     [~, lambda]=boxcox(rt)
    
       
    tbl1=table(zscore(response_time_2AFC(:,h)), zscore(val_sum_rat_abs(:,h)), zscore(val_2AFC_abs_out(:,h)),'VariableNames',{'RT','sum','diff'});
    tbl2=table(zscore(response_time_rating(:,h)), zscore(val_rat(:,h)), zscore(val_rat_abs(:,h)),'VariableNames',{'RT','val','sal'});
    
    
    lm_rat{h} =fitlm(tbl2,'RT ~ val+ sal');
    lm_ratInter{h} =fitlm(tbl2,'RT ~ val* sal');   
    lm_salienceonly{h} =fitlm(tbl2,'RT ~ sal');    
    lm_2afc{h}=fitlm(tbl1,'RT ~ diff');

    % Sum is not in the model when combined with difference, so not used
      

    valin=val_rat(:,h);  
    ind=find(valin>=0);
    valence=zeros(size(valin));
    valence(ind)=1;
    magnitude= abs(valin);
    
    
    % standarized coefficients when a dummy variable is present shoud be
    % divided by 2*std, not only one std. Dummy stay the same, paper about
    % standarization with dummy variables. 
    % http://www.stat.columbia.edu/~gelman/research/published/standardizing7.pdf
    RT2AFC_corr_for_dummy=(response_time_rating(:,h)- mean(response_time_rating(:,h)))/ (2*std(response_time_rating(:,h)));
    magnitude_corr_for_dummy=(magnitude - mean(magnitude))/(2*std(magnitude));
    
    tbl3=table(RT2AFC_corr_for_dummy, magnitude_corr_for_dummy , valence,'VariableNames',{'RT','magnit','valen'});
    lm_MagValen{h} =fitlm(tbl3,'RT ~ magnit * valen');
   
    
    h=h+1;
    
    
    if sum(magnitude==0)>=1
        error('There is one answer that equals zero!! check your regrtession for RT and delete these values while comparing pos and neg domain')
    end 
end



% Wuensch paper analysis, ZPF test for comparisson correlations (Pearson),
% not sure if correct with Spearmann

% 2AFC
  
  [ZPF(1),pvalZPF(1)]=ZPF_test(response_time_2AFC, val_2AFC_abs_out, type);
%   [ZPF(2),pvalZPF(2)]=ZPF_test(response_time_2AFC, val_abs_rat, type);
%   [ZPF(3),pvalZPF(3)]=ZPF_test(response_time_2AFC, val_abs_2AFC, type);
%   [ZPF(4),pvalZPF(4)]=ZPF_test(response_time_2AFC, val_abs_rank, type);
%   [ZPF(5),pvalZPF(5)]=ZPF_test(response_time_2AFC, val_abs_total, type);
%   
  
% Rating

  [ZPFRat(1), pvalZPFrat(1)]= ZPF_test(response_time_rating, val_rat, type);
  [ZPFRat(2), pvalZPFrat(2)]= ZPF_test(response_time_rating, val_rat_abs, type);
%   [ZPFRat(3), pvalZPFrat(3)]= ZPF_test(response_time_rating, val_rat_rat, type);
%   [ZPFRat(4), pvalZPFrat(4)]= ZPF_test(response_time_rating, val_rat_2AFC, type);
%   [ZPFRat(5), pvalZPFrat(5)]= ZPF_test(response_time_rating, val_rat_rank, type);
%   [ZPFRat(6), pvalZPFrat(6)]= ZPF_test(response_time_rating, abs(val_rat_rat), type);
%   [ZPFRat(7), pvalZPFrat(7)]= ZPF_test(response_time_rating, abs(val_rat_2AFC), type);
%   [ZPFRat(8), pvalZPFrat(8)]= ZPF_test(response_time_rating, abs(val_rat_rank), type);
%   [ZPFRat(9), pvalZPFrat(9)]= ZPF_test(response_time_rating, abs(val_rat_2AFC), type);
%   [ZPFRat(10), pvalZPFrat(10)]= ZPF_test(response_time_rating, abs(val_rat_rank), type);
    



% Wuensch paper analysis, modified ttest, mainly a
% Welsh t-test with dof=n1+n2-2m-2; where m is the number of estimations, 
% the intercept is not taken into account in m. Again in the paper.


  dof= 190 + 190 - 2 -2;
  [tbetas(1), pbetas(1)]=ttestbetas(lm_2afc, dof, 2);

  dof= 60 + 60 - (2*2) -2;
  [tbetas(2), pbetas(2)]=ttestbetas(lm_rat, dof, 2);
  
  dof= 60 + 60 - (2*2) -2;
  [tbetas(3), pbetas(3)]=ttestbetas(lm_rat, dof, 3);


% Wuensch paper analysis,  Pothoff anlysis, hierarchic regression, allows
% for comparisson of intrcepts and slopes of both regressions. It is more
% effective than the ttest because it does not use the betas but the data,
% so no rounding. Moreover more flexible as it allows to shift the
% intercept to check effects ad different levels. 
  
% Problem, I use each variable in the Rating separately, despite both are
% included in the model salinecy and value. Maybe wrong? Check it out!
  
  [pottlm{1}, compFval{1}, interPval{1}]=Potthoff_analysis(response_time_2AFC, val_2AFC_abs_out);
  [pottlm{2}, compFval{2}, interPval{2}]=Potthoff_analysis(response_time_rating, val_rat);
  [pottlm{3}, compFval{3}, interPval{3}]=Potthoff_analysis(response_time_rating, val_rat_abs);
  [pottlm{4}, compFval{4}, interPval{4}]=Potthoff_MagnitudeValence(response_time_rating, val_rat);

  % the fourth analysis is an analysis of the relatinship of magnitude vs
  % valence using Pottoff, so magnitude as continuous and valence as dummy
  % varaiable, including data both pre and post, but without analysing if
  % pre post are different in this case. 

  
  %% 

  % Vincentized group RT distributions analysis
  
  % link.springer.com/article/10.3758%2Fs13423-011-0053-5
  % this consists in constructing a quantilized distribution of RTs for
  % each condition, for example two conditions, positive and negative
  % valence. Then, for each individual, bin the RT and find the mean of the
  % each bin. Then go for group level and find the mean of each bin time. 
  
  
  
  %% simeon recommended using GLZ models, so a mixed model with an gamma link (inverse function)
  % many papers suggest this, actually, in DDM, as the accumulator of
  % evidence receives an innovation with normal distribution, the latencies
  % follow an inverse normal distribution (LATER model of RT), so what
  % Simeon said is quite useful. Some other people use Exgaussian
  % distribution (a mixture between normal and exponential) and other a
  % lognormal distribution for RT.
  % https://www.ncbi.nlm.nih.gov/pubmed/16754533
  
  % some assumptions of gamma distribution are that the ration
  % mean/variance is constant, thus the higher the value, the higher the
  % variance in this part of the value.  http://civil.colorado.edu/~balajir/CVEN6833/lectures/GammaGLM-01.pdf
  % moreover it should be used only por distributions between  bigger than zero and
  % +infinite, never negative values and never zero.  
  
  % R example of gamma glm http://civil.colorado.edu/~balajir/CVEN6833/lectures/GammaGLM-01.pdf
  % The residual deviance for the GLM-Gamma  mo del in section 8.1 was
  % 599.09 on 798 degrees of freedom with an AIC  of 2479.4. The GLM-Normal
  % mo del has deviance of 2156.6 and the AIC  as 3069.6 
  
  % so, for example run the LMM you made and the same with gamma function,
  % compare AIC, BIC devaince... and keep the best one.
  
  % more literature on 
  % link.springer.com/article/10.3758%2Fs13423-011-0053-5
  % http://www.ane.pl/pdf/4320.pdf
  % http://www.sciencedirect.com/science/article/pii/0022249665900143
  % http://cogprints.org/6603/1/recin-psychrev.pdf
  
  % the Recinormal distribution, so the inverse of the nornal is what the
  % RT could follow
  
  % nature neuroscience papaer on inverse gamma of RT http://www.nature.com/neuro/journal/v3/n8/full/nn0800_827.html
  % and a website about this LATER theory of RT http://www.cudos.ac.uk/later.htm
  
end 




