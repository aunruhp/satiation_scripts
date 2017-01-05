

f={};
try
    cd('/media/Projects/Alex/Satiation Analysis');      %/Volumes/MTL/CS4 output')   %orginal form michael '/Volumes/Barnacle/Data/Satiation_data/Satiation_session_data'
catch
    cd('/media/Projects/Alex/');      %/Volumes/MTL/CS4 output')
end


h=1;
DirName= dir;
for a = 1:length(DirName)
    if regexp(DirName(a).name, 'Sat_+[0-9][0-9][0-9]')
        f{end+1,1} = [DirName(a).name];
        
        
        %
        % if nargin ==0
        %     [satcue, ranking, events] = load_behavioural_data;
        % elseif nargin>0
        [satcue, ranking, events] = load_behavioural_data(DirName(a).name);
        % end
        
        
        [RT] = extract_reaction_times(events);
        
        
        
        for jj=[3,6]
            val_rat=nan(length(satcue{jj-1}(:,14)),1);             % value of image on A given in previous rating trials (-300 to 300)
            val_2AFC=nan(length(satcue{jj-1}(:,14)),1);            % value of image on A given in previous whole 2AFC trials (0 to 1900)
            val_ranking=nan(length(satcue{jj-1}(:,14)),1);         % rank of image on A, calculated after all the choice trials (0 to 20)
            val_total=nan(length(satcue{jj-1}(:,14)),1);         % rank of image on A, calculated after all the choice trials (0 to 20)
            
            
            
            for i=1:20
                
                ind=find(cell2mat(satcue{jj-1}(:,2))==i);            % when this stimulus was on left
                
                
                valRat= ranking{jj-1}(i,4);                        % rating value (-300 to 300)
                valRank= ranking{jj}(i,3);                         % ranking (0 to 20)
                val2AFC= ranking{jj}(i,4);                         % choice rating (0-1900)
                valTot= ranking{jj}(i,5);                         % choice rating (-1900-1900)
                
                
                val_rat(ind)= valRat;
                val_L_ranking(ind)= valRank;
                val_L_2AFC(ind)= val2AFC;
                val_L_total(ind)= valTot;
                
                product_value(i,:) = cell2mat(satcue{jj-1}(ind,14));% this is to reckon the variance of each product, then all variances sorted by expectesd value as the rest. More random-like as no pool of different products (which is more fixed-like)
            end
            
            
            % Correlation between value given and value input difference?
            
            
            val_expected=val_rat;
            %% partial_error analysis
            
            % this is a very weak anlysis. As patients have very few
            % error, compardd to monkeys, and Diff-Drift predicts reduction
            % of error with higher difference in value
            % for this I use this analysis, which maps the response (in 9 levels)
            % to the difference in value (9 levels from 0 to 600)
            
            
            % Expected value
            
            expected_val=[];
            
            step=600/9; shift=600/2;
            
            for i=1:length(val_expected)
                if  val_expected(i)< step*1 -shift
                    expected_val(i)=1;
                elseif val_expected(i)< step*2 -shift
                    expected_val(i)=2;
                elseif val_expected(i)< step*3 -shift
                    expected_val(i)=3;
                elseif val_expected(i)< step*4 -shift
                    expected_val(i)=4;
                elseif val_expected(i)< step*5 -shift
                    expected_val(i)=5;
                elseif val_expected(i)< step*6 -shift
                    expected_val(i)=6;
                elseif val_expected(i)< step*7 -shift
                    expected_val(i)=7;
                elseif val_expected(i)< step*8 -shift
                    expected_val(i)=8;
                elseif val_expected(i)<= step*9 -shift
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
            resp_val=cell2mat(satcue{jj-1}(:,14));
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
            
            
            for i=1:9
                
                final_coeff(i)= coeff_error(i)/ length(find(expected_val==i));
                fincorral_freq(i) = freq_corrError (i)  / length(find(expected_val==i));
                final_noncorr_freq(i) = non_corrFreq (i) / length(find(expected_val==i));
                fincorrCoef(i)  = corrCoef  (i) / length(find(expected_val==i));
            end
            
            
            out_coeff{h,jj/3}=final_coeff;
            out_corrfreq{h,jj/3}= fincorral_freq;
            out_noncorr_freq{h,jj/3}=final_noncorr_freq;
            out_corrCoeff{h,jj/3}= fincorrCoef;
            
            
            manhattan_distance_comp=nan(9,60);
            manhattan_distance_level=nan(9,1);
            euclidian_level=nan(9,1);
            for i=1:9
                manhattan_distance_comp(i,1:length(find(expected_val==i))) = (val_expected(expected_val==i))/3 - (resp_val(expected_val==i));  % this means, I distribute the values based on their expected values/ divided by 3 , and then I
                manhattan_distance_level(i) = mean(abs((val_expected(expected_val==i))/3 - (resp_val(expected_val==i))));
                euclidian_level(i) = sqrt(sum(((val_expected(expected_val==i))/3 - (resp_val(expected_val==i))).^2))/(length(find(expected_val==i)));
              
            end

            out_manhattan_distance{h,jj/3}= manhattan_distance_level';
            out_manhattan_distance_components{h,jj/3}= manhattan_distance_comp';
            out_euclidean_distance{h,jj/3}= euclidian_level';
            
           
            
            % variance/std product-wise analysis
            
            product_variance = std(product_value, 0, 2);
            product_sum = sum(product_value,  2);
                        
            
            expected_val_product=[];  % In this case this is going to be a sorting vector, with 20 indeces form 1 to 9 to sort them
            
            step=600/9; shift=600/2;
            
            for i=1:length(product_variance)
                if  product_sum(i)< step*1 -shift
                    expected_val_product(i)=1;
                elseif product_sum(i)< step*2 -shift
                    expected_val_product(i)=2;
                elseif product_sum(i)< step*3 -shift
                    expected_val_product(i)=3;
                elseif product_sum(i)< step*4 -shift
                    expected_val_product(i)=4;
                elseif product_sum(i)< step*5 -shift
                    expected_val_product(i)=5;
                elseif product_sum(i)< step*6 -shift
                    expected_val_product(i)=6;
                elseif product_sum(i)< step*7 -shift
                    expected_val_product(i)=7;
                elseif product_sum(i)< step*8 -shift
                    expected_val_product(i)=8;
                elseif product_sum(i)<= step*9 -shift
                    expected_val_product(i)=9;
                else
                    expected_val_product(i)=0;
                end
            end
            
            if any(expected_val_product==0);
                error('some expected values were not classified');
            end

            
            for kk=1:9
                variance_mean(kk) = mean(product_variance(expected_val_product==kk));

            end
            
            out_variance_level{h,jj/3} = variance_mean;
            
        end
        h=h+1;
    end
    
end

for jj=[1:7,9:11]
    
    A1(jj,1:9)=out_coeff{jj,1}; C1(jj,1:9)=out_corrfreq{jj,1}; B1(jj,1:9)=out_noncorr_freq{jj,1}; D1(jj,1:9)=out_corrCoeff{jj,1}; E1(jj,1:9)=out_manhattan_distance{jj,1}; F1(jj,1:9)=out_euclidean_distance{jj,1}; G1(jj,1:9)=out_variance_level{jj,1};
    A2(jj,1:9)=out_coeff{jj,2}; C2(jj,1:9)=out_corrfreq{jj,2}; B2(jj,1:9)=out_noncorr_freq{jj,2}; D2(jj,1:9)=out_corrCoeff{jj,2}; E2(jj,1:9)=out_manhattan_distance{jj,2}; F2(jj,1:9)=out_euclidean_distance{jj,2}; G2(jj,1:9)=out_variance_level{jj,2};
end

% beforeA=nansum(A1,1)/10; beforeB=nansum(B1,1)/10; beforeC=nansum(C1,1)/10; beforeD=nansum(D1,1)/10;
% afterA=nansum(A2,1)/10; afterB=nansum(B2,1)/10; afterC=nansum(C2,1)/10; afterD=nansum(D2,1)/10;

beforeA=nanmean(A1,1); beforeC=nanmean(C1,1); beforeB=nanmean(B1,1); beforeD=nanmean(D1,1); beforeE=nanmean(E1,1); beforeF=nanmean(F1,1); beforeG=nanmean(G1,1); 
afterA=nanmean(A2,1); afterC=nanmean(C2,1); afterB=nanmean(B2,1); afterD=nanmean(D2,1); afterE=nanmean(E2,1); afterF=nanmean(F2,1); afterG=nanmean(G2,1);       




figure;
bar([beforeA;afterA]);
title('Partial Error Coefficient'); ylabel('Mean Partial Error Coefficient')

figure;
bar([beforeB;afterB]);
title('Partial Error Frequency');ylabel('Mean Partial  Error Frequency' )

figure;
bar([beforeC;afterC]);
title('Partial Error Frequency Corrected');ylabel('Mean Partial Error Frequency' )

figure;
bar([beforeD;afterD]);
title('Partial Error Coefficient Corrected');ylabel('Mean Partial Error Coefficient' )

figure;
bar([beforeE; afterE]);
title('Mean Absolute Error (Manhattan distance)');ylabel('Mean Manhattan distance')

figure;
bar([beforeF; afterF]);
title('Mean Euclidean distance');ylabel('Mean euclidean distance')

figure;
bar([beforeG; afterG]);
title('Mean Variance Error');ylabel('Mean Standard Deviation')


%% Inference

% in this part you have to select which error signal is used, partial
% error, absolute distance, variance...



% comparisson pre post of the errors, individually, so, signrank test for
% error distribution before and after

i=1;
for h=[1:6,9:11]  
    [pval_ErrorDistrib(i)] = signrank(A1(h,:), A2(h,:));
    i=i+1;
end


% correlation and regression analysis

% Patients 37 and 39 excluded in this analysis



% subject-wise

i=1;

for h=[1:6,9:11]
    A1s=A1(h,:);  valtocorr= (-100:25:100); dummytodelete = [0,0,0,0, nan ,1,1,1,1]; dummytocorr = [0,0,0,0,1,1,1,1]; 


    sValCorrA1(i)= corr(A1s(~isnan(A1s))', valtocorr(~isnan(A1s))');
    sSalCorrA1(i)= corr(A1s(~isnan(A1s))', abs(valtocorr(~isnan(A1s)))');
    
    A2s=A2(h,:);
    sValCorrA2(i)= corr(A2s(~isnan(A2s))', valtocorr(~isnan(A2s))');
    sSalCorrA2(i)= corr(A2s(~isnan(A2s))', abs(valtocorr(~isnan(A2s)))');
    
    [p1{i}]=fitlm([zscore(valtocorr(~isnan(A1s)))', zscore(abs(valtocorr(~isnan(A1s))))'], zscore(A1s(~isnan(A1s)))');
    [p2{i}]=fitlm([zscore(valtocorr(~isnan(A2s)))', zscore(abs(valtocorr(~isnan(A2s))))'], zscore(A2s(~isnan(A2s)))');
    
    betasVal1(i)= p1{i}.Coefficients.Estimate(2);
    betasSal1(i)= p1{i}.Coefficients.Estimate(3);
    betasVal2(i)= p2{i}.Coefficients.Estimate(2);
    betasSal2(i)= p2{i}.Coefficients.Estimate(3);
    
    
    
    A1sdummy=A1s(~isnan(dummytodelete)); A2sdummy=A2s(~isnan(dummytodelete));  % exclude the 0 values to compare positive and negative domain
    
    tab1 = table([dummytocorr(~isnan(A1sdummy))'],[zscore(abs(valtocorr(~isnan(A1sdummy))))'], zscore(A1sdummy(~isnan(A1sdummy)))', 'VariableNames', {'valence','sal', 'A1'});

    tab2 = table([dummytocorr(~isnan(A2sdummy))'],[zscore(abs(valtocorr(~isnan(A2sdummy))))'], zscore(A2sdummy(~isnan(A2sdummy)))', 'VariableNames', {'valence','sal', 'A2'});


    [p3{i}]=fitlm(tab1, 'A1~valence*sal');
    [p4{i}]=fitlm(tab2, 'A2~valence*sal');
    
    betasVal3(i)= p3{i}.Coefficients.Estimate(2);
    betasSal3(i)= p3{i}.Coefficients.Estimate(3);
    betasSalVal3(i)= p3{i}.Coefficients.Estimate(4);
    betasVal4(i)= p4{i}.Coefficients.Estimate(2);
    betasSal4(i)= p4{i}.Coefficients.Estimate(3);
    betasSalVal4(i)= p4{i}.Coefficients.Estimate(4);
    
    
    % standarized coefficients when a dummy variable is present shoud be
    % divided by 2*std, not only one std. Dummy stay the same, paper about
    % standarization with dummy variables. 
    % http://www.stat.columbia.edu/~gelman/research/published/standardizing7.pdf
    
    
    std_coeff_Sal_dummy1=(abs(valtocorr(~isnan(A1sdummy))) - mean(abs(valtocorr(~isnan(A1sdummy)))))/(2*std(abs(valtocorr(~isnan(A1sdummy)))));
    std_coeff_Sal_dummy2=(abs(valtocorr(~isnan(A2sdummy))) - mean(abs(valtocorr(~isnan(A2sdummy)))))/(2*std(abs(valtocorr(~isnan(A2sdummy)))));
    std_coeff_A1_dummy1=(A1sdummy(~isnan(A1sdummy))- mean(A1sdummy(~isnan(A1sdummy)))) / (2* std(A1sdummy(~isnan(A1sdummy))));
    std_coeff_A2_dummy2=(A2sdummy(~isnan(A2sdummy))- mean(A2sdummy(~isnan(A2sdummy)))) / (2* std(A2sdummy(~isnan(A2sdummy))));
    
    tab3 = table(dummytocorr(~isnan(A1sdummy))', std_coeff_Sal_dummy1', std_coeff_A1_dummy1', 'VariableNames', {'valence','sal', 'A1'});

    tab4 = table(dummytocorr(~isnan(A2sdummy))', std_coeff_Sal_dummy2', std_coeff_A2_dummy2', 'VariableNames', {'valence','sal', 'A2'});

    [p5{i}]=fitlm(tab3, 'A1~valence*sal');
    [p6{i}]=fitlm(tab4, 'A2~valence*sal');


    betasVal5(i)= p5{i}.Coefficients.Estimate(2);
    betasSal5(i)= p5{i}.Coefficients.Estimate(3);
    betasSalVal5(i)= p5{i}.Coefficients.Estimate(4);
    betasVal6(i)= p6{i}.Coefficients.Estimate(2);
    betasSal6(i)= p6{i}.Coefficients.Estimate(3);
    betasSalVal6(i)= p6{i}.Coefficients.Estimate(4);

    i=i+1;
end

disp(' on sample ttests to compare value betas, if significant')
[~,pval_Val1]=ttest(betasVal1),[~,pval_Val2]=ttest(betasVal2)
disp(' on sample ttests to compare valence betas, if significant')
[~,pval_Valence3]=ttest(betasVal3),[~,pval_Valence4]=ttest(betasVal4),[~,pval_Valence5]=ttest(betasVal5),[~,pval_Valence6]=ttest(betasVal6)
disp(' on sample ttests to compare salience betas, if significant')
[~,pval_Sal1]=ttest(betasSal1),[~,pval_Sal2]=ttest(betasSal2),[~,pval_Mag3]=ttest(betasSal3),[~,pval_Mag4]=ttest(betasSal4),[~,pval_Mag5]=ttest(betasSal5),[~,pval_Mag6]=ttest(betasSal6)
disp(' on sample ttests to compare interactions magnitude valence betas, if significant')
[~,pval_MagValence3]=ttest(betasSalVal3),[~,pval_MagValence4]=ttest(betasSalVal4),[~,pval_MagValence5]=ttest(betasSalVal5),[~,pval_MagValence6]=ttest(betasSalVal6)


% test before and after

disp('paired ttests to compare before vs after , if significant')
[~,pval_Sal_whole]=ttest (betasSal1,    betasSal2)
[~,pval_Val_whole]=ttest (betasVal1,    betasVal2)
[~,pval_Mag_whole]=ttest (betasSal5,    betasSal6)
[~,pval_MagValence_whole]=ttest (betasSalVal5, betasSalVal6)


% individual ttest pre-post (welsch correction) with Bonferoni correction 

% Wuensch paper analysis, modified ttest, mainly a
% Welsh t-test with dof=n1+n2-2m-2; where m is the number of estimations, 
% the intercept is not taken into account in m. Again in the paper.


dof= 9 + 9 - 2 -2; % 9 values in each regression, 2 estimations

i=1;
for h=1:length(p1)
    
    [tbetasVal(i), pbetasVal(i)]=ttestbetas_for_ErrorAna(p1{h}, p2{h}, dof, 2); %the last 2 is the position of the estimate
    [tbetasSal(i), pbetasSal(i)]=ttestbetas_for_ErrorAna(p1{h}, p2{h}, dof, 3);
    [tbetasMag(i), pbetasMag(i)]=ttestbetas_for_ErrorAna(p5{h}, p6{h}, dof-1, 3); % 1 dof less as the interaction is calculated
    [tbetasMagValence(i), pbetasMagValence(i)]=ttestbetas_for_ErrorAna(p5{h}, p6{h}, dof-1, 4);  % 1 dof less as the interaction is calculated
    i=i+1;
    
end
   
disp(' subject-wise welsch ttest, pre vs post, bonferoni corrected, if 1, a significant difference')
pbetasVal<0.05/length(p1),pbetasSal<0.05/length(p1),pbetasMag<0.05/length(p1), pbetasMagValence<0.05/length(p1)

% no Potthoff analysis is possible here, as multivariable case. But if I consider this as only Salience, then it is possible 

