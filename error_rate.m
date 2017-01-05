function error_rate

% in this analysis I will check if there is a correlation between the error
% rate (intransitivity) and the difference in value in 2AFC. In principle,
% the higher the difference, the less intransitivity, according to
% Diff-Drift models, as random effects are less effective. 
% the deal is: should I do linear regression and compare slopes or better
% correlation and compare its coefficients? 
% mmain problem, there are very few intrnasitivities, thus, maybe not very
% useful this analysis. 



[satcue, ranking] = load_behavioural_data;



val_rat_rat= nan(length(satcue{2}(:,14)),2);
val_rat_2AFC=nan(length(satcue{2}(:,14)),2);
val_rat_rank=nan(length(satcue{2}(:,14)),2);

for i=1:20
    ind=find(cell2mat(satcue{2}(:,2))==i);
    valRat=  ranking{2}(i,4); % ranking
    val2AFC= ranking{3}(i,4);
    valRank= ranking{3}(i,3);
    
    val_rat_rat(ind,1)= valRat;
    val_rat_2AFC(ind,1)= val2AFC;
    val_rat_rank(ind,1)= valRank;
    
    ind2=find(cell2mat(satcue{5}(:,2))==i);
    valRat=  ranking{5}(i,4); % ranking
    val2AFC= ranking{6}(i,4);
    valRank= ranking{6}(i,3);
    val_rat_rat(ind2,2)= valRat;
    val_rat_2AFC(ind2,2)= val2AFC;
    val_rat_rank(ind2,2)= valRank;
    
    
end        




h=1;
for jj=[3]  % 3 is before pause, and 6 after it
    
    
   
    
    %% Choice trials analysis (2AFC)
    
    
    
    val_A_rat=nan(length(satcue{jj}(:,14)),1);             % value of image on A given in previous rating trials (-300 to 300)
    val_B_rat=nan(length(satcue{jj}(:,14)),1);             % value of image on B given in previous rating trials (-300 to 300)
    val_left_rat=  nan(length(satcue{jj}),1);              % value of image on the left given in previous rating trials (-300 to 300)
    val_right_rat= nan(length(satcue{jj}),1);              % value of image on the right given in previous rating trials (-300 to 300)
    val_A_2AFC=nan(length(satcue{jj}(:,14)),1);            % value of image on A given in previous whole 2AFC trials (0 to 1900)
    val_B_2AFC=nan(length(satcue{jj}(:,14)),1);            % value of image on B given in previous whole 2AFC trials (0 to 1900)
    val_A_ranking=nan(length(satcue{jj}(:,14)),1);         % rank of image on A, calculated after all the choice trials (0 to 20) 
    val_B_ranking=nan(length(satcue{jj}(:,14)),1);         % rank of image on B, calculated after all the choice trials (0 to 20) 
    
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
        
                
        val_A_rat(indA)= valRat;
        val_B_rat(indB)= valRat;
        val_left_rat(ind)= valRat;
        val_right_rat(ind2)= valRat;
        val_A_ranking(indA)= valRank;
        val_B_ranking(indB)= valRank;
        val_A_2AFC(indA)= val2AFC;
        val_B_2AFC(indB)= val2AFC;
        
        
    end
end

val_diff= val_right_rat - val_left_rat; %if r>l, result is positive, if r<l, negative
val_abs= abs(val_diff);
[rho2AFCdiff, pval2AFCdiff]=corr(cell2mat(satcue{3}(:,14)), val_diff, 'type', 'Spearman');  % the higher the value difference, the higher should be the value given, so if R>L, then both positive
[rho2AFCabs, pval2AFCabs]=corr(abs(cell2mat(satcue{3}(:,14))), val_abs, 'type', 'Spearman');

figure; plot (val_abs, abs(cell2mat(satcue{3}(:,14))),'o')



%% intransitivities analysis 
% 
% [pval_intrans, intransitivity_ratio, diff_intran]=intransitivities_analysis;

%% partial_error analysis

% this is a very weak anlysis. As patients have very few
% error, compardd to monkeys, and Diff-Drift predicts reduction
% of error with higher difference in value
% for this I use this analysis, which maps the response (in 9 levels)
% to the difference in value (9 levels from 0 to 600)

% level_val=[];
% step= 600/9;
% for i=1:length(val_abs)
%     if     val_abs(i)< step*1
%             level_val(i)=1;
%     elseif val_abs(i)< step*2
%             level_val(i)=2;
%     elseif val_abs(i)< step*3
%             level_val(i)=3;
%     elseif val_abs(i)< step*4
%             level_val(i)=4;
%     elseif val_abs(i)< step*5
%             level_val(i)=5;
%     elseif val_abs(i)< step*6
%             level_val(i)=6;
%     elseif val_abs(i)< step*7
%             level_val(i)=7;
%     elseif val_abs(i)< step*8
%             level_val(i)=8;
%     else
%             level_val(i)=9;
%     end
% end
% 
% 
% level_val2=[];
% resp_val=cell2mat(satcue{3}(:,14));
% for i=1:length(resp_val)
%     if     resp_val(i)== -100
%             level_val2(i)=1;
%     elseif resp_val(i)== -75
%             level_val2(i)=2;
%     elseif resp_val(i)== -50
%             level_val2(i)=3;
%     elseif resp_val(i)== -25
%             level_val2(i)=4;
%     elseif resp_val(i)== 0
%             level_val2(i)=5;
%     elseif resp_val(i)== 25
%             level_val2(i)=6;
%     elseif resp_val(i)== 50
%             level_val2(i)=7;
%     elseif resp_val(i)== 75
%             level_val2(i)=8;
%     elseif resp_val(i)== 100
%             level_val2(i)=9;
%     end
% end


level_val=[];
step= 600/5;
for i=1:length(val_abs)
    if     val_abs(i)< step*1
            level_val(i)=1;
    elseif val_abs(i)< step*2
            level_val(i)=2;
    elseif val_abs(i)< step*3
            level_val(i)=3;
    elseif val_abs(i)< step*4
            level_val(i)=4;
    else %if val_abs(i)< step*5
            level_val(i)=5;
%     elseif val_abs(i)< step*6
%             level_val(i)=6;
%     elseif val_abs(i)< step*7
%             level_val(i)=7;
%     elseif val_abs(i)< step*8
%             level_val(i)=8;
%     else
%             level_val(i)=9;
    end
end


level_val2=[];
resp_val=cell2mat(satcue{3}(:,14));
for i=1:length(resp_val)
    if     resp_val(i)== -100
            level_val2(i)=5;
    elseif resp_val(i)== -75
            level_val2(i)=4;
    elseif resp_val(i)== -50
            level_val2(i)=3;
    elseif resp_val(i)== -25
            level_val2(i)=2;
    elseif resp_val(i)== 0
            level_val2(i)=1;
    elseif resp_val(i)== 25
            level_val2(i)=2;
    elseif resp_val(i)== 50
            level_val2(i)=3;
    elseif resp_val(i)== 75
            level_val2(i)=4;
    elseif resp_val(i)== 100
            level_val2(i)=5;
    end
end


partial_error= level_val - level_val2;
partial_error_abs= abs(partial_error);
partial_error_corr = partial_error_abs >= 2;

coeff_error=[];
freq_error=[];
for kk=1:5
coeff_error(kk)= sum(partial_error_abs(level_val==kk));
freq_error(kk)= sum(partial_error_corr(level_val==kk));
end

%% Find the identities of left and right images during 2AFC

% ident_left= cell2mat(satcue{3}(:,2));
% ident_right= cell2mat(satcue{3}(:,3));
% 
% %% Create looping variables and 


% j=1;
% for i= ident_left'
%     rank1(j,1)= ranking{2}(i,4); %for each left image, check ID and include the ranking obtained during evaluation
%     j=j+1;
% end
% 
% 
% j=1;
% for i= ident_right'
%     rank2(j,1)= ranking{2}(i,4);
%     j=j+1;
% end
% 

% rank_difference=rank1>rank2;
% value=cell2mat(satcue{3}(:,14));
% diff_intran=[];
% intrans=0;
% 
% h=1;
% for i=1:length(rank_difference)
%     if rank_difference(i)==1
%         if value(i)>0
%             disp(['intransitivity found in ' num2str(i)]);
%             intrans=intrans +1;
%             diff_intran(h) = abs(rank1(i)-rank2(i));
%             h=h+1;
%         end
%     end
%     if rank_difference(i)==0
%         if value(i)<0
%             disp(['intransitivity found in ' num2str(i)]);
%             intrans=intrans +1;
%             diff_intran(h)= abs(rank1(i)-rank2(i));
%             h=h+1;
%         end
%     end
% end
% 
% intransitivity_ratio=intrans/length(ident_left);
% pval= sum(binopdf(intrans:length(ident_left), length(ident_left), 0.05));


end



