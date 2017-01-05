function [rankBias_number_Rat, rankBias_number_2afc, rankBias_number_Tot,meanBias_number_Rat, meanBias_number_2afc,meanBias_number_Tot]=attention_bias(DirName)

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


dbstop if error


%% Load the data
if nargin ==0
    [satcue, ranking] = load_behavioural_data;
elseif nargin>0
     [satcue, ranking] = load_behavioural_data(DirName);
end

%% Create a matrix that converts the original left-right layout in A-B, to compare the influence of time A (2 seconds) and B (1 second) and/or order of stimulation


% This loop creates a matrix with the position when a specific stimulus happend.
% Converts columns from left-right to A-B. It also changes the sign of the
% value when the stimlus displayed on left, as negative means  actually
% positive for that stimlus ID.

rankBias_number_Rat={''};
rankBias_number_2afc={''};
rankBias_number_Tot={''};

meanBias_number_Rat={''};
meanBias_number_2afc={''};
meanBias_number_Tot={''};


h=1;
for jj=[3,6]

rankBias_Rat=0;
rankBias_2afc=0;
rankBias_Tot=0;

meanBias_Rat=0;
meanBias_2afc=0;
meanBias_Tot=0;

for i=1:20
    
    indL=find(cell2mat(satcue{jj}(:,2))==i);        % when in the 2AFC was this image diplayed on the left
    indR=find(cell2mat(satcue{jj}(:,3))==i);        % when in the 2AFC was this image diplayed on the rigth
    lr_to_ab=cell2mat(satcue{jj}(indL,8));          % conversion factor from left-right to A-B for columns on left   %rmember if ==1 left=A/right=B, if ==zero otherwise
    lr_to_ab2=cell2mat(satcue{jj}(indR,8));         % conversion factor from left-right to A-B for columns on right
    val_L= - cell2mat(satcue{jj}(indL,14));      % values, sign is reversed, because as it is on the left, if negative it is positive actually
    val_R=  cell2mat(satcue{jj}(indR,14));      % values, not reversed as when on right, if positive it is positive for this stim, and if negative as well 
    ext_ind=[indL' zeros(1,length(indR))]';        % extension of each column for concatenatin
    ext_ind2=[zeros(1,length(indL)) indR']';
    mat=[ext_ind ext_ind2 vertcat(lr_to_ab, lr_to_ab2) vertcat(val_L, val_R)];   % concatenation in a matrix
    mat((mat(1:length(indL),3)==0),[1,2])                      =  mat((mat(1:length(indL),3)==0),[2,1]);                          %  Reverse the position of the columns when the conv factor is 0 
    mat(find(mat(length(indL)+1:end,3)==0)+length(indL),[1,2])  =  mat(find(mat(length(indL)+1:end,3)==0)+length(indL),[2,1]);      %  Reverse the position of the columns when the conv factor is 0
    
    values_A=mat(mat(:,1)~=0,4);            % when in first column is an index, then the stimulus was on A
    values_B=mat(mat(:,2)~=0,4);            % when in second column is an index, then the stimulus was on B
    
 

    
    if ranking{jj-1}(i,4)> 0                                             % This follows Rangel's model, in which if positive value, the longer the time payed attention, the more value coded
        [pval_pos(i),h1]=ranksum(values_A, values_B, 'tail', 'right');   % therefore, if the value is positive, the probablity of choosing it when A is higher than in B, as more time and more value is coded, thus A>B as H1
        if h1>0
            disp(['ID ', DirName(5:7), ' stimulus '  , num2str(i), ' has a significant attentional bias Rating'])
            rankBias_Rat=rankBias_Rat+1;
        end
    end
    if ranking{jj-1}(i,4)< 0                                             % This follows Rangel's model, in which if negative value, the longer the time payed attention, the less value coded
        [pval_neg(i),h2]=ranksum(values_A, values_B, 'tail', 'left');    % therefore, if the value is negative, the probablity of choosing it when A is lower than in B, as more time and more rejection is coded, , thus A<B as H1
        if h2>0
            disp(['ID ', DirName(5:7), ' stimulus ' , num2str(i), ' has a significant attentional bias Rating'])
            rankBias_Rat=rankBias_Rat+1;
        end
    end
    
    if ranking{jj}(i,4)> 950                                            % The same but for the classification of pos/neg I use the final value during 2AFC choices, if > half the top value, it is positive
        [pval_pos_2AFC(i),h3]=ranksum(values_A, values_B, 'tail', 'right');   % therefore, if the value is positive, the probablity of choosing it when A is higher than in B, as more time and more value is coded, thus A>B as H1
        if h3>0
            disp(['ID ', DirName(5:7), ' stimulus '  , num2str(i), ' has a significant attentional bias 2AFC'])
            rankBias_2afc=rankBias_2afc+1;
        end
    end
    if ranking{jj}(i,4)< 950                                            % The same but for the classification of pos/neg I use the final value during 2AFC choices, if < half the top value, it is negative
        [pval_neg_2AFC(i),h4]=ranksum(values_A, values_B, 'tail', 'left');    % therefore, if the value is negative, the probablity of choosing it when A is lower than in B, as more time and more rejection is coded, , thus A<B as H1
        if h4>0
            disp(['ID ', DirName(5:7), ' stimulus ' , num2str(i), ' has a significant attentional bias 2AFC'])
            rankBias_2afc=rankBias_2afc+1;
        end
    end
    
    if ranking{jj}(i,5)> 0                                          % The same but for the classification of pos/neg I use the final value during 2AFC choices, if > half the top value, it is positive
        [pval_pos_2AFC(i),h3]=ranksum(values_A, values_B, 'tail', 'right');   % therefore, if the value is positive, the probablity of choosing it when A is higher than in B, as more time and more value is coded, thus A>B as H1
        if h3>0
            disp(['ID ', DirName(5:7), ' stimulus '  , num2str(i), ' has a significant attentional bias Total'])
            rankBias_Tot=rankBias_Tot+1;
        end
    end
    if ranking{jj}(i,5)< 0                                           % The same but for the classification of pos/neg I use the final value during 2AFC choices, if < half the top value, it is negative
        [pval_neg_2AFC(i),h4]=ranksum(values_A, values_B, 'tail', 'left');    % therefore, if the value is negative, the probablity of choosing it when A is lower than in B, as more time and more rejection is coded, , thus A<B as H1
        if h4>0
            disp(['ID ', DirName(5:7), ' stimulus ' , num2str(i),' has a significant attentional bias Total'])
            rankBias_Tot=rankBias_Tot+1;
        end
    end
    
    
    
    %% Here no ranksum test, just if it is higher or not
    
    
    if ranking{jj-1}(i,4)> 0                                             % This follows Rangel's model, in which if positive value, the longer the time payed attention, the more value coded
        if (mean(values_A)> mean(values_B))   % therefore, if the value is negative, the probablity of choosing it when A is lower than in B, as more time and more rejection is coded, , thus A<B as H1
            meanBias_Rat=meanBias_Rat+1;
        end
    end
    if ranking{jj-1}(i,4)< 0                                             % This follows Rangel's model, in which if negative value, the longer the time payed attention, the less value coded
        if (mean(values_A)< mean(values_B))   % therefore, if the value is negative, the probablity of choosing it when A is lower than in B, as more time and more rejection is coded, , thus A<B as H1
            meanBias_Rat=meanBias_Rat+1;
        end
        
    end
    
    if ranking{jj}(i,4)> 950                                            % The same but for the classification of pos/neg I use the final value during 2AFC choices, if > half the top value, it is positive
        if (mean(values_A)> mean(values_B))   % therefore, if the value is negative, the probablity of choosing it when A is lower than in B, as more time and more rejection is coded, , thus A<B as H1
            meanBias_2afc=meanBias_2afc+1;
        end
    end
    if ranking{jj}(i,4)< 950                                            % The same but for the classification of pos/neg I use the final value during 2AFC choices, if < half the top value, it is negative
        if (mean(values_A)< mean(values_B))   % therefore, if the value is negative, the probablity of choosing it when A is lower than in B, as more time and more rejection is coded, , thus A<B as H1
            meanBias_2afc=meanBias_2afc+1;
        end
    end
    
    if ranking{jj}(i,5)> 0                                          % The same but for the classification of pos/neg I use the final value during 2AFC choices, if > half the top value, it is positive
        if (mean(values_A)> mean(values_B))   % therefore, if the value is negative, the probablity of choosing it when A is lower than in B, as more time and more rejection is coded, , thus A<B as H1
            meanBias_Tot=meanBias_Tot+1;
        end
    end
    if ranking{jj}(i,5)< 0                                           % The same but for the classification of pos/neg I use the final value during 2AFC choices, if < half the top value, it is negative
        if (mean(values_A)< mean(values_B))   % therefore, if the value is negative, the probablity of choosing it when A is lower than in B, as more time and more rejection is coded, , thus A<B as H1
            meanBias_Tot=meanBias_Tot+1;
        end
    end
    
end
   
rankBias_number_Rat{h}=rankBias_Rat;
rankBias_number_2afc{h}=rankBias_2afc;
rankBias_number_Tot{h}=rankBias_Tot;
 
meanBias_number_Rat{h}=meanBias_Rat;
meanBias_number_2afc{h}=meanBias_2afc;
meanBias_number_Tot{h}=meanBias_Tot;
 


h=h+1;
   


%                                                  % This is just to make sure I am not loosing an important positive test in the opposite direction                                              
%    [pval_2tail(i),h5]=ranksum(values_A, values_B, 'tail', 'both');
%    if h5>0
%        disp(['ID ' , num2str(i), ' has a significant attentional bias'])
%    end
%    
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
% %  
% %         
% %     %% Any tendency towards what was shown left or right during 2AFC likert?  Just in case. 
% %     
%     [pval_lr(i), h6] = ranksum (val_L, val_R);
%     if h6>0
%         disp(['ID ' , num2str(i), ' has a significat bias due to the position on screen when shown at choice time (left vs right)'])
%     end
%     
%     
% %     
% %     %% Analysis of the effect of A-B and L-R on the value given to each stimulus. 
% %     % Anova of the explanation of varaice of value by temporal stimulation
% %     % or position during choice, and if any interaction
% %     % a 2-way ANOVA with groups A-B and L-R
% %     
% %     % Be careful, I have seen sometimes where all val_RA have -100 and all
% %     % val_LB = 100, this looks weird to me... so maybe check this also in several ways
% %     
% % 
% %         
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
% %     %% Factorial analysis for A-B vs L-R to check any dependency
% %     % For this I do Pearson's Chi-squared test for independency of factors
% %     % in a contingency table, H0 is independency, so that the position does not change the frequency of selection of this stimulus 
% %     % H1 is dependency, but
% %     % but this is only for the number of times every options was chosen, so the
% %     % freqeuncy/probability of choosing one or other (length of each value)
% %     % Moreover Fisher exact test is performed, a non-parametric, as the
% %     % samples are lower than 20, and not every cell has at least 5, but
% %     % close... 
% %     
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





end