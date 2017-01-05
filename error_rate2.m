function [out_coeff, out_corrfreq, out_noncorr_freq, out_corrCoeff, corrDiffError, pvalDiffError, corrInOut, pvalcorrInOut, corrDiffErrorAbs, pvalDiffErrorAbs, corrInOutAbs, pvalcorrInOutAbs,corrRTError, pvalcorrRTError,corrRTErrorFreq, pvalcorrRTErrorFreq]= error_rate2(DirName)

% in this analysis I will check if there is a correlation between the error
% rate (intransitivity) and the difference in value in 2AFC. In principle,
% the higher the difference, the less intransitivity, according to
% Diff-Drift models, as random effects are less effective.
% the deal is: should I do linear regression and compare slopes or better
% correlation and compare its coefficients?
% mmain problem, there are very few intrnasitivities, thus, maybe not very
% useful this analysis.



if nargin ==0
    [satcue, ranking, events] = load_behavioural_data;
elseif nargin>0
    [satcue, ranking, events] = load_behavioural_data(DirName);
end


[RT] = extract_reaction_times(events);



h=1;
for jj=[3,6]  % 3 is before pause, and 6 after it
    

    val_L_rat=nan(length(satcue{jj}(:,14)),1);             % value of image on A given in previous rating trials (-300 to 300)
    val_R_rat=nan(length(satcue{jj}(:,14)),1);             % value of image on B given in previous rating trials (-300 to 300)
    val_L_2AFC=nan(length(satcue{jj}(:,14)),1);            % value of image on A given in previous whole 2AFC trials (0 to 1900)
    val_R_2AFC=nan(length(satcue{jj}(:,14)),1);            % value of image on B given in previous whole 2AFC trials (0 to 1900)
    val_L_ranking=nan(length(satcue{jj}(:,14)),1);         % rank of image on A, calculated after all the choice trials (0 to 20)
    val_R_ranking=nan(length(satcue{jj}(:,14)),1);         % rank of image on B, calculated after all the choice trials (0 to 20)
    val_L_total=nan(length(satcue{jj}(:,14)),1);         % rank of image on A, calculated after all the choice trials (0 to 20)
    val_R_total=nan(length(satcue{jj}(:,14)),1);         % rank of image on B, calculated after all the choice trials (0 to 20)
    
    
    for i=1:20
        
        indL=find(cell2mat(satcue{jj}(:,2))==i);            % when this stimulus was on left
        indR=find(cell2mat(satcue{jj}(:,3))==i);           % when this stimulus was on right

        valRat= ranking{jj-1}(i,4);                        % rating value (-300 to 300)
        valRank= ranking{jj}(i,3);                         % ranking (0 to 20)
        val2AFC= ranking{jj}(i,4);                         % choice rating (0-1900)
        valTot= ranking{jj}(i,5);                         % choice rating (-1900-1900)
        
        
        val_L_rat(indL)= valRat;
        val_R_rat(indR)= valRat;
        val_L_ranking(indL)= valRank;
        val_R_ranking(indR)= valRank;
        val_L_2AFC(indL)= val2AFC;
        val_R_2AFC(indR)= val2AFC;
        val_L_total(indL)= valTot;
        val_R_total(indR)= valTot;

    end
    
    
    % Correlation between value given and value input difference?

    val_diff= val_R_rat - val_L_rat ; %if r>l, result is positive, if r<l, negative
    val_abs= abs(val_diff);

    
    %% partial_error analysis
    
    % this is a very weak anlysis. As patients have very few
    % error, compardd to monkeys, and Diff-Drift predicts reduction
    % of error with higher difference in value
    % for this I use this analysis, which maps the response (in 9 levels)
    % to the difference in value (9 levels from 0 to 600)
    

    % Expected value
    
    expected_val=[];

    step=1200/9; shift=1200/2;
    
    for i=1:length(val_abs)
        if  val_diff(i)< step*1 -shift
            expected_val(i)=1;
        elseif val_diff(i)< step*2 -shift
            expected_val(i)=2;
        elseif val_diff(i)< step*3 -shift
            expected_val(i)=3;
        elseif val_diff(i)< step*4 -shift
            expected_val(i)=4;
        elseif val_diff(i)< step*5 -shift
            expected_val(i)=5;
        elseif val_diff(i)< step*6 -shift
            expected_val(i)=6;
        elseif val_diff(i)< step*7 -shift
            expected_val(i)=7;
        elseif val_diff(i)< step*8 -shift
            expected_val(i)=8;
        elseif   val_diff(i)<= step*9 -shift
            expected_val(i)=9;
        else 
            expected_val(i)=0;
        end
    end
    
    if any(expected_val==0); 
        error('some expected values were not classified');
    end 
    
    % Real value
    
    real_val=[];
    resp_val=cell2mat(satcue{jj}(:,14));
    for i=1:length(resp_val)
        if     resp_val(i)== -100
            real_val(i)=1;
        elseif resp_val(i)== -75
            real_val(i)=2;
        elseif resp_val(i)== -50
            real_val(i)=3;
        elseif resp_val(i)== -25
            real_val(i)=4;
        elseif resp_val(i)== 0
            real_val(i)=5;
        elseif resp_val(i)== 25
            real_val(i)=6;
        elseif resp_val(i)== 50
            real_val(i)=7;
        elseif resp_val(i)== 75
            real_val(i)=8;
        elseif resp_val(i)== 100
            real_val(i)=9;

        end
    end
    
    if any(resp_val==0); 
        error('some answer values were not classified')
    end 
    
    % Resdiuals!
    
    partial_error= expected_val - real_val;
    partial_error_abs= abs(partial_error);       %
    partial_error_freq = partial_error_abs >= 1; %
    partial_error_freqcorr = partial_error_abs >= 2; %
    partial_corrCoef= partial_error_abs; partial_corrCoef(partial_corrCoef==1)=0;     % not included as error if only one step of difference
    
    
    coeff_error=[];
    freq_corrError=[];
    non_corrFreq=[];
    error_total=cell('');
    pos_error=cell('');
    neg_error=cell('');
    corrCoef=[];
    
    for kk=1:9
        coeff_error(kk)= sum(partial_error_abs(expected_val==kk));
        freq_corrError(kk) = sum(partial_error_freqcorr(expected_val==kk));
        non_corrFreq(kk)= sum(partial_error_freq(expected_val==kk));
        corrCoef(kk)= sum(partial_corrCoef(expected_val==kk));
        error_total{kk}  = partial_error(expected_val==kk);
        neg_error{kk}= error_total{kk}(error_total{kk}<0);
        pos_error{kk}= error_total{kk}(error_total{kk}>0);
        
    end
    
    
%     for i=1:5
%         
%         final_coeff(6-i)= (coeff_error(i) + coeff_error(10-i)) / (length(find(expected_val==i))+ length(find(expected_val==10-i)));
%         fincorral_freq(6-i) = (freq_corrError (i) + freq_corrError (10-i)) / (length(find(expected_val==i))+ length(find(expected_val==10-i)));
%         final_noncorr_freq(6-i) = (non_corrFreq (i) + non_corrFreq (10-i)) / (length(find(expected_val==i))+ length(find(expected_val==10-i)));
%         fincorrCoef(6-i)  = (corrCoef  (i) + corrCoef  (10-i)) / (length(find(expected_val==i))+ length(find(expected_val==10-i)));
%     end

    for i=1:9

        final_coeff(i)= coeff_error(i)/ length(find(expected_val==i));
        fincorral_freq(i) = freq_corrError (i)  / length(find(expected_val==i));
        final_noncorr_freq(i) = non_corrFreq (i) / length(find(expected_val==i));
        fincorrCoef(i)  = corrCoef (i) / length(find(expected_val==i));
    end


    
    out_coeff{h}=final_coeff;
    out_corrfreq{h}= fincorral_freq;
    out_noncorr_freq{h}=final_noncorr_freq;
    out_corrCoeff{h}= fincorrCoef;
    
    [corrDiffError(h), pvalDiffError(h)]=corr(val_diff(:), partial_error_abs(:), 'type', 'Spearman');
    [corrDiffErrorAbs(h), pvalDiffErrorAbs(h)]=corr(val_abs(:), partial_error_abs(:), 'type', 'Spearman');
    [corrInOut(h), pvalcorrInOut(h)]=corr(cell2mat(satcue{jj}(:,14)), val_diff, 'type', 'Spearman');  % the higher the value difference, the higher should be the value given, so if R>L, then both positive
    [corrInOutAbs(h), pvalcorrInOutAbs(h)]=corr(abs(cell2mat(satcue{jj}(:,14))), val_abs, 'type', 'Spearman');
    [corrRTError(h,1), pvalcorrRTError(h,1)]=corr(partial_error_abs(:), RT{jj}, 'type', 'Spearman');  % the higher the value difference, the higher should be the value given, so if R>L, then both positive
    [corrRTErrorFreq(h,1), pvalcorrRTErrorFreq(h,1)]=corr(partial_error_freq(:), RT{jj}, 'type', 'Spearman');
    [corrRTError(h,2), pvalcorrRTError(h,2)]=corr(partial_corrCoef(:), RT{jj}, 'type', 'Spearman');  % the higher the value difference, the higher should be the value given, so if R>L, then both positive
    [corrRTErrorFreq(h,2), pvalcorrRTErrorFreq(h,2)]=corr(partial_error_freqcorr(:), RT{jj}, 'type', 'Spearman');
    
    
    h=h+1;
    
end


end


