function decision_bias

% In principle, according to rangle, if an object has positive value, and is
% paid more attention than other one, its relative value increases,
% moreover, if is negative, more negative. So in my case, A shown 2x longer
% than B, so I expect that if A has positive value, it will have a tendency
% to be selected and if negative, it will have a tendency not to be selected
% compared to B, just because of the time computing.
%
% In our case, it is not clear if the time on screen is going to affect the
% decision, better not, as only the combinations (not the permutations)
% where displayed. So, no perfect control as if a pair of images is shown 
% once in one order, it is not repeated in the oppossite order. 
%
% Output: 
% pvals of every test performed
%


[satcue, ranking] = load_behavioural_data;



%% Create a matrix that converts the original left-right layout in A-B, to compare the influence of time A (2 seconds) and B (1 second) and/or order of stimulation


% This loop creates a matrix with the position when a specific stimulus happend.
% Converts columns from left-right to A-B. It also changes the sign of the
% value when the stimlus displayed on left, as negative means  actually
% positive for that stimlus ID.


for i=1:20
    
    indL=find(cell2mat(satcue{3}(:,2))==i);        % when in the 2AFC was this image diplayed on the left
    indR=find(cell2mat(satcue{3}(:,3))==i);        % when in the 2AFC was this image diplayed on the rigth
    lr_to_ab=cell2mat(satcue{3}(indL,8));          % conversion factor from left-right to A-B for columns on left   %rmember if ==1 left=A/right=B, if ==zero otherwise
    lr_to_ab2=cell2mat(satcue{3}(indR,8));         % conversion factor from left-right to A-B for columns on right
    val_L= - cell2mat(satcue{3}(indL,14));      % values, sign is reversed, because as it is on the left, if negative it is positive actually
    val_R=  cell2mat(satcue{3}(indR,14));      % values, not reversed as when on right, if positive it is positive for this stim, and if negative as well 
    ext_ind=[indL' zeros(1,length(indR))]';        % extension of each column for concatenatin
    ext_ind2=[zeros(1,length(indL)) indR']';
    mat=[ext_ind ext_ind2 vertcat(lr_to_ab, lr_to_ab2) vertcat(val_L, val_R)];   % concatenation in a matrix
    mat((mat(1:length(indL),3)==0),[1,2])                      =  mat((mat(1:length(indL),3)==0),[2,1]);                          %  Reverse the position of the columns when the conv factor is 0 
    mat(find(mat(length(indL)+1:end,3)==0)+length(indL),[1,2])  =  mat(find(mat(length(indL)+1:end,3)==0)+length(indL),[2,1]);      %  Reverse the position of the columns when the conv factor is 0
    
    values_A=mat(mat(:,1)~=0,4);            % when in first column is an index, then the stimulus was on A
    values_B=mat(mat(:,2)~=0,4);            % when in second column is an index, then the stimulus was on B
    

    
%     if ranking{2}(i,4)> 0                                             % This follows Rangel's model, in which if positive value, the longer the time payed attention, the more value coded
%         [pval_pos(i),h1]=ranksum(values_A, values_B, 'tail', 'right');   % therefore, if the value is positive, the probablity of choosing it when A is higher than in B, as more time and more value is coded, thus A>B as H1
%         if h1>0
%             disp(['ID ' , num2str(i), ' has a significant attentional bias'])
%         end
%     end
%     if ranking{2}(i,4)< 0                                             % This follows Rangel's model, in which if negative value, the longer the time payed attention, the less value coded
%         [pval_neg(i),h2]=ranksum(values_A, values_B, 'tail', 'left');    % therefore, if the value is negative, the probablity of choosing it when A is lower than in B, as more time and more rejection is coded, , thus A<B as H1
%         if h2>0
%             disp(['ID ' , num2str(i), ' has a significant attentional bias'])
%         end
%     end
%     
%     if ranking{3}(i,4)> 950                                            % The same but for the classification of pos/neg I use the final value during 2AFC choices, if > half the top value, it is positive
%         [pval_pos_2AFC(i),h3]=ranksum(values_A, values_B, 'tail', 'right');   % therefore, if the value is positive, the probablity of choosing it when A is higher than in B, as more time and more value is coded, thus A>B as H1
%         if h3>0
%             disp(['ID ' , num2str(i), ' has a significant attentional bias'])
%         end
%     end
%     if ranking{3}(i,4)< 950                                            % The same but for the classification of pos/neg I use the final value during 2AFC choices, if < half the top value, it is negative
%         [pval_neg_2AFC(i),h4]=ranksum(values_A, values_B, 'tail', 'left');    % therefore, if the value is negative, the probablity of choosing it when A is lower than in B, as more time and more rejection is coded, , thus A<B as H1
%         if h4>0
%             disp(['ID ' , num2str(i), ' has a significant attentional bias'])
%         end
%     end
%                                                  % This is just to make sure I am not loosing an important positive test in the opposite direction                                              
%    [pval_2tail(i),h5]=ranksum(values_A, values_B, 'tail', 'both');
%    if h5>0
%        disp(['ID ' , num2str(i), ' has a significant attentional bias'])
%    end
%    
%     
%    
%     % These are some alternative anlaysis, taken individually the value
%     % given in the likert as outcome to check if positive outcomes are more
%     % positive when in A than in B, so only locally, trial by trial anlysis. 
%     % I am not sure if this analysis makes any sense, but probably not
%     
%     
%     %only positive values
%     posmat   = mat(mat(:,4)>0,[1:4]);
%     values_A_pos = posmat(posmat(:,1)~=0,4);
%     values_B_pos = posmat(posmat(:,2)~=0,4);
%     if ~isempty(values_A_pos) && ~isempty(values_B_pos)
%         [pval_pos_alternative(i),h1]=ranksum(values_A_pos, values_B_pos, 'tail', 'right'); %higher in A if positive as longer time;
%     end
%     
%     %only negative values
%     negmat       = mat(mat(:,4)<0,[1:4]);
%     values_A_neg = negmat(negmat(:,1)~=0,4);
%     values_B_neg = negmat(negmat(:,2)~=0,4);
%     if ~isempty(values_A_neg) && ~isempty(values_B_neg)
%         [pval_neg_alternative(i),h2]=ranksum(values_A_neg, values_B_neg, 'tail', 'left'); %lower in A as negative and longer in time
%     end
%     if h1+h2>0
%         disp(['ID ' , num2str(i), ' has a significat attentional bias'])
%     end
%  
%         
%     %% Any tendency towards what was shown left or right during 2AFC likert?  Just in case. 
%     
%     [pval_lr(i), h6] = ranksum (val_L, val_R);
%     if h6>0
%         disp(['ID ' , num2str(i), ' has a significat bias due to the position on screen when shown at choice time (left vs right)'])
%     end
%     
%     
%     
%     %% Analysis of the effect of A-B and L-R on the value given to each stimulus. 
%     % Anova of the explanation of varaice of value by temporal stimulation
%     % or position during choice, and if any interaction
%     % a 2-way ANOVA with groups A-B and L-R
%     
%     % Be careful, I have seen sometimes where all val_RA have -100 and all
%     % val_LB = 100, this looks weird to me... so maybe check this also in several ways
%     
% 
%         
%     indLA=(lr_to_ab) .* indL;  indLA(indLA==0)=[];   % indeces where stimulus was shown as A and on the left
%     indLB=(~lr_to_ab) .* indL; indLB(indLB==0)=[];   % indeces where stimulus was shown as B and on the left
%     indRA=(~lr_to_ab2).* indR; indRA(indRA==0)=[];   % indeces where stimulus was shown as A and on the right
%     indRB=(lr_to_ab2).* indR;  indRB(indRB==0)=[];   % indeces where stimulus was shown as B and on the right
%     val_LA= - cell2mat(satcue{3}(indLA,14));         % values of when stimulus was shown as A and on the left, signed reversed so that it is relative to this stimulus
%     val_LB= - cell2mat(satcue{3}(indLB,14));         % values of when stimulus was shown as B and on the left, signed reversed so that it is relative to this stimulus
%     val_RA= cell2mat(satcue{3}(indRA,14));           % values of when stimulus was shown as A and on the right
%     val_RB= cell2mat(satcue{3}(indRB,14));           % values of when stimulus was shown as B and on the right
%     val=[val_LA; val_LB; val_RA; val_RB];            % all the values i none vector for ANOVA, first left, second right. 
%     g1=[ones(length(val_L),1); zeros(length(val_R),1)];              % group, with 1=left, 0=right.  
%     g2=[ones(length(val_LA),1); zeros(length(val_LB),1); ones(length(val_RA),1); zeros(length(val_RB),1)]; % for 1=A, 0=B. Intercalated
%     
%     [p,t, stats, terms] = anovan(val, {g1,g2}, 'varnames', {'L-R', 'A-B'}, 'display', 'off');
%     [p2,t2, stats2, terms2] = anovan(val, {g1,g2}, 'varnames',{'L-R', 'A-B'}, 'model', 'interaction', 'display', 'off');
%     
%     
%     if any(p<0.05)
%         disp(['ID ' , num2str(i), ' has a significat bias in single stimulus anova, stimulus' num2str(i) ])
%     end
%     
%     if any(p2<0.05)
%         disp(['ID ' , num2str(i), ' has a significat bias in single stimulus anova, stimulus' num2str(i) ])
%     end
%     
% 
%     %% Factorial analysis for A-B vs L-R to check any dependency
%     % For this I do Pearson's Chi-squared test for independency of factors
%     % in a contingency table, H0 is independency, so that the position does not change the frequency of selection of this stimulus 
%     % H1 is dependency, but
%     % but this is only for the number of times every options was chosen, so the
%     % freqeuncy/probability of choosing one or other (length of each value)
%     % Moreover Fisher exact test is performed, a non-parametric, as the
%     % samples are lower than 20, and not every cell has at least 5, but
%     % close... 
%     
%     
%         
%     fr_LA=length(find(cell2mat(satcue{3}(indLA,14))<0)) + length(find(cell2mat(satcue{3}(indRB,14))<0));
%     fr_RB=length(find(cell2mat(satcue{3}(indRB,14))>0)) + length(find(cell2mat(satcue{3}(indLA,14))>0));
%     fr_RA=length(find(cell2mat(satcue{3}(indRA,14))>0)) + length(find(cell2mat(satcue{3}(indLB,14))>0));
%     fr_LB=length(find(cell2mat(satcue{3}(indLB,14))<0)) + length(find(cell2mat(satcue{3}(indRA,14))<0));
%     cont_tab=[fr_LA, fr_LB; fr_RA, fr_RB];
%     
%     [p_chi, x2_chi] = chisquarecont(cont_tab);                 % parametric test for contingency tables, to test in any nonrandom association, chi-squared
%     [h_fisher,p_fisher,stats_fisher] = fisherexacttest(cont_tab);   % non parametric test for contingency tables, to test in any nonrandom association, fisher exact test
%     
%     if (p_chi)<0.05
%         disp(['ID ' , num2str(i), ' has a significat bias in single stimulus X2, stimulus' num2str(i) ])
%     end
%     
%     if (p_fisher)<0.05
%         disp(['ID ' , num2str(i), ' has a significat bias in single stimulus Fisher, stimulus' num2str(i) ])
%     end
    
end



%% Here I analyse not stimulus by stimulus, but general tendencies.

% left vs right in general, are values when on left different than on right?
% in this case, positive values are when chosen right, negative when
% chosen left. I test if, when chosen left, a tendency to give them
% more or less value than whrn chosen the one on the right
% also if the products on left are chosen more often than on the right
% this is different to the anlysis made in th stimulus, there the value
% was of, given the stimulus, which value it got in R vs L. Here is,
% given chosen something on the left, des it have different values than
% when chosen one of the right, so, if a tendency to give products
% on one side more value. If positive, mean, for example, that on right
% closer to zero when selected than on the left. So more value to left.

left_value  = abs(cell2mat(satcue{3}((cell2mat(satcue{3}(:,14))<0), 14)));  % absolute value of negative ones
right_value = cell2mat(satcue{3}((cell2mat(satcue{3}(:,14))>0), 14));
[pval_lr_gen, h7]= ranksum(left_value, right_value);                       % if any difference is found, the patient has tendency to give more value to products on one side or another when chosen?

if h7>0
    disp('This patient has a significat bias due to the position on screen when shown at choice time (left vs right)')
% elseif abs(length(left_value) - length(right_value)) / 190 > 0.05
%     disp('This patient has a significat bias due to the position on screen when shown at choice time (left vs right), left or right were selected more than 0.05 times more often')
end

%% Binomial test, I made this but n general the binomial test is for small samples, when more samples,
% there are better option such as likelihood ratio test or wald-test

% Binomial test to test if left chosen more than right, analysing freq not
% value given, trials 190, successes = number of
% chosen left, expected probability = 0.5

left_chosen= length(left_value);
if left_chosen <= 190/2
    p_bino_l = sum(binopdf([1:left_chosen], 190, 0.5));
    if p_bino_l< 0.05
    disp('This patient has a significat bias due to the position on screen when shown at choice time (left vs right), two-tail binomial test')
    end
elseif left_chosen > 190/2
    p_bino_r = sum(binopdf([left_chosen:190], 190, 0.5));
    if p_bino_r<0.05
    disp('This patient has a significat bias due to the position on screen when shown at choice time (left vs right), two-tail binomial test')
    end
end


%% likelihood-ratio tests

% relative likelihood-ratio test, to complement binomial test.
% by Neyman-Pearson lemma, the likelihood-ratio test is the most poweful
% test, so the highest sensitivity, at a specific alpha, if all paramters
% of the model are known. Thus if this test is negative, none will be
% positive. 

uMLE= mle(left_chosen, 'distribution', 'binomial', 'ntrials', 190); 

uLogL=  sum(log(binopdf(left_chosen, 190, uMLE))); 
rLogL=  sum(log(binopdf(left_chosen, 190, 0.5)));

[h_lratioLR, p_lratioLR, RatioLR, CriticalValueLR] = lratiotest(uLogL, rLogL,1);  % make the test likelihood ratio with the negative loglikelihood and 1 degree of freedom. 
  

if h_lratioLR>0
    disp('This patient has a significat bias in likelihood ratio test')
end



% in this case the unrestricted model is not unrestriucted, as not less parameters but hust different. 
% 
% % manual version 
% succ= left_chosen; %number of successes
% trials=190;
% dof=1; % for simple hypothesis testing of parameters 1 dof is used. When complex models are used (composite hypotheis testing) the difference in parameters is used. 
% p0=  0.5;
% 
% stats_lratio= 2* (succ * log(uMLE/p0) + (trials-succ)* log((1-uMLE)/ (1- p0)));
% pValue_lratio = 1-chi2cdf(stats_lratio,dof);




%% Real val A vs B (if value given == 0, no choice so not analysed)


% Same as previous analysis but for A vs B when chosen


% this loop extracts the information relative to A, B and if on R or on L 

indRB=[]; indLA=[];  indLB=[]; indRA=[];
fr_LBgen= 0; fr_RBgen= 0; fr_RAgen= 0; fr_LAgen= 0;

for jj= 1:190
   if cell2mat(satcue{3}(jj,8))==1                       % if conversion factor == 1, so L=A, R=B
        if cell2mat(satcue{3}(jj,14))<0                  % if negative, thus left and A (LA)
            valA(jj)  = - cell2mat(satcue{3}(jj,14));    % absolute value for A
            valLA(jj) = - cell2mat(satcue{3}(jj,14));    % absolute value for LA
            indLA=[indLA; jj];                           % indices for LA
            fr_LAgen= fr_LAgen + 1;                      % how many times chosen LA
        elseif cell2mat(satcue{3}(jj,14))>0              % if positive, thus right and B (RB) 
            valB(jj) = cell2mat(satcue{3}(jj,14));
            valRB(jj) = cell2mat(satcue{3}(jj,14));
            indRB=[indRB; jj];
            fr_RBgen= fr_RBgen + 1;
        end
   elseif cell2mat(satcue{3}(jj,8))==0                  % if conversion factor == 1, so L=B, R=A
        if cell2mat(satcue{3}(jj,14))>0                 % if positive, thus right and A (RA)     
            valA(jj) = cell2mat(satcue{3}(jj,14));
            valRA(jj) = cell2mat(satcue{3}(jj,14));
            indRA=[indRA; jj];
            fr_RAgen= fr_RAgen + 1;
        elseif cell2mat(satcue{3}(jj,14)) < 0           % if negative, thus left and B (LB)
            valB(jj) = - cell2mat(satcue{3}(jj,14));
            valLB(jj) = - cell2mat(satcue{3}(jj,14));
            indLB=[indLB; jj];
            fr_LBgen= fr_LBgen + 1;
        end
   end
end
valA(valA==0)=[]; valRA(valRA==0)=[]; valLA(valLA==0)=[];  % delete zeros
valB(valB==0)=[]; valLB(valLB==0)=[]; valRB(valRB==0)=[];


% Rank sum analysis for difference, if positive tendency t ogive more value
% to A than to B

[pval_AB_gen, h8]  = ranksum(valA,valB);


if h8>0
    disp('This patient has a significat bias ranksum (A vs B)')
% elseif abs(length(valA) - length(valB)) / 190 > 0.05
%     disp('This patient has a significat bias due to the position on screen when shown at choice time (A vs B), left or right were selected more than 0.05 times more often')
end



%crossed-test to check interactions
[pval_inter, h9]  = ranksum([valLA,valRB], [valLB, valRA]);


if h9>0
    disp('This patient has a significat interaction A-B,L-R, rank sum test')
end




% Binomial test to test if A chosen more than B, analysing freq not
% value given, trials 190, successes = number of
% chosen left, expected probability = 0.5

A_chosen= length(valA);
if A_chosen <= 190/2
    p_bino_l = sum(binopdf([1:A_chosen], 190, 0.5));
    if p_bino_l< 0.05 
        disp('This patient has a significat bias due to the position on screen when shown (A vs B), one-tail binomial test')
    end
elseif A_chosen > 190/2
    p_bino_r = sum(binopdf([A_chosen:190], 190, 0.5));
    if  p_bino_r<0.05
        disp('This patient has a significat bias due to the position on screen when shown (A vs B), one-tail binomial test')
    end
end



%% likelihood-ratio tests

% relative likelihood-ratio test, to complement binomial test.
% by Neyman-Pearson lemma, the likelihood-ratio test is the most poweful
% test, so the highest sensitivity, at a specific alpha, if all paramters
% of the model are known. Thus if this test is negative, none will be
% positive. 

uMLE= mle(A_chosen, 'distribution', 'binomial', 'ntrials', 190); 

uLogL=  sum(log(binopdf(A_chosen, 190, uMLE))); 
rLogL=  sum(log(binopdf(A_chosen, 190, 0.5)));

[h_lratioAB, p_lratioAB, RatioAB, CriticalValueAB] = lratiotest(uLogL, rLogL,1);  % make the test likelihood ratio with the negative loglikelihood and 1 degree of freedom. 
                                                                          % in this case the unrestricted model is not unrestriucted, as not less parameters but hust different. 
if h_lratioAB>0
    disp('This patient has a significat bias in likelihood ratio test')
end

% % manual version 
% succ= A_chosen; %number of successes
% trials=190;
% dof=1; % for simple hypothesis testing of parameters 1 dof is used. When complex models are used (composite hypotheis testing) the difference in parameters is used. 
% p0=  0.5;
% 
% stats_lratio= 2* (succ * log(uMLE/p0) + (trials-succ)* log((1-uMLE)/ (1- p0)));
% pValue_lratio = 1-chi2cdf(stats_lratio,dof);
% 




%% General analysis of the effect of A-B and L-R on the value given to each stimulus.
% Anova of the explanation of varaice of value by temporal stimulation
% or position during choice, and if any interaction


g1gen([indRA; indRB])=1;   %indeces for left == 1 in g1
g1gen([indLA; indLB])=2;   %indeces for right == 2 in g1
g2gen([indLA; indRA])=1;   %indeces for A == 1 in g2
g2gen([indLB; indRB])=2;   %indeces for B == 2 in g2
val= abs(cell2mat(satcue{3}(:,14))); %just the values given

[p3,t3, stats3, terms3] = anovan(val, {g1gen,g2gen}, 'varnames', {'L-R', 'A-B'});
[p4,t4, stats4, terms4] = anovan(val, {g1gen,g2gen}, 'varnames', {'L-R', 'A-B'}, 'model', 'interaction');


if any(p3<0.05)
    disp(' This patient has a significat bias in whole stimulus anova' )
end

if any(p4<0.05)
    disp( 'This patient has a significat interaction in whole stimulus anova' )
end


% 
% [p5,t5, stats5, terms5] = kruskalwallis(val, {g1gen,g2gen}, 'varnames', {'L-R', 'A-B'});
% [p6,t6, stats6, terms6] = kruskalwallis(val, {g1gen,g2gen}, 'varnames', {'L-R', 'A-B'}, 'model', 'interaction');
%    


% Internet he ANOVA procedure with fixed factors and equal sample sizes works quite well even when the assumption of normality is violated, unless one or more of the distributions are highly skewed or the variances are very different. 


%% Problem, the distribution of the values is not normal, actually not eevan close as 0 value is approx the mean in all patients but
% the 0 and around are usually the least frquent ones. Morover, no
% continuous, this decreases the potency of the test for interaction, even
% though it can be quite good for main effects, thus here I perfomr the
% aligned-rank anova, which transfomrs the data to a rank-normalized
% distribution before the normal anova for interactions. 
% it controls for FA, and has more power to detect interaction even if both
% main effects are present.
% Because ranking is a nonlinear transformation, in these data it removes the significant interaction. Aligned ranking fixes this problem.
% http://depts.washington.edu/aimgroup/proj/art/,  aligned rank tests (ART) in Higgins (2003) Introduction to Modern Nonparametric Statistics, Duxbury Press. Higgins gives the data example:
% http://faculty.washington.edu/wobbrock/pubs/chi-11.06.pdf
% Example: If you have two factors (X1 and X2), and the response (Y), you will run three ANOVAs, each using the same input model (X1, X2, X1×X2), but using a different response variable, one for each aligned-and-ranked Y. That is, one ANOVA will use the response for which Y was aligned-and-ranked for X1. The second ANOVA will use the response for which Y was aligned-and-ranked for X2. The third ANOVA will use the response for which Y was aligned-and-ranked for X1×X2. When interpreting the results in each ANOVA's output, only look at the main effect or interaction for which Y was aligned and ranked. So you would extract one result from each of three ANOVAs, for three total results. 


% Arttool citation
%     Kay M and Wobbrock J (2016). ARTool: Aligned Rank Transform for Nonparametric Factorial ANOVAs. R package version 0.10.4, https://github.com/mjskay/ARTool.
% 
%     Wobbrock J, Findlater L, Gergle D and Higgins J (2011). “The Aligned Rank Transform for Nonparametric Factorial Analyses Using Only ANOVA Procedures.” In Proceedings of the ACM Conference on Human Factors in Computing Systems (CHI '11), pp. 143–146. http://depts.washington.edu/aimgroup/proj/art/.




%% Factorial analysis for A-B vs L-R to check any dependency on frequency of choice
% For this I do Pearson's Chi-squared test for independency of factors
% in a contingency table, H0 is independency, so that the position does not change the frequency of selection of this stimulus

cont_tab=[fr_LAgen, fr_LBgen; fr_RAgen, fr_RBgen];
[p_chigen, x2_chigen] = chisquarecont(cont_tab);                 % parametric test for contingency tables, to test in any nonrandom association, chi-squared
[h_fishergen,p_fishergen,stats_fishergen] = fisherexacttest(cont_tab);   % non parametric test for contingency tables, to test in any nonrandom association, fisher exact test

% [pval,wald_stat]=Barnardextest(fr_LAgen, fr_LBgen, fr_RAgen, fr_RBgen); %this is a new test (1950's or so, which is brand new for math!) and does a NON-parametric 2x2 contingency tables test, but more sensitive than fiser exact test and less restrctions on experimental conditions, quite good
    
if (p_chigen)<0.05
    disp(['ID ' , num2str(i), ' has a significat bias in p_chigen' ])
end

if (p_fishergen)<0.05
    disp(['ID ' , num2str(i), ' has a significat bias in fisherexacttest' ])
end

% if any(pval)<0.05
%     disp(['ID ' , num2str(i), ' has a significat bias in Barnardextest'])
% end
% The test was first published by George Alfred Barnard (1945) (link to the original paper in Nature). 
% Mehta and Senchaudhuri (2003) explain why Barnard’s test can be more powerful than Fisher’s under certain conditions:



end