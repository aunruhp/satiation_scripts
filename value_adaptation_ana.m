function value_adaptation_ana

    

        
        % this analysis thies to test if the value of the products affects
        % subsewuent values, such as value adaptation in context dependent Decision
        % making. Thus I make a window of rating/rankings, and I expect to see that
        
        
        dbstop if error
        [path]= uigetdir('/media/Projects/Alex/Satiation Analysis', 'Select a patient folder');
        cd (path)
        %load behavioural results
        pathfolders=dir;
        for i=1:length(pathfolders)
            if regexp(pathfolders(i).name, 'Files')
                cd (pathfolders(i).name)
                filesfolder=dir;
                for j=1:length(filesfolder)
                    if regexp(filesfolder(j).name, 'runsat_output')
                        load (filesfolder(j).name)
                    end
                end
            end
            
        end
        
        
        %cd (path)
        
        %correct the errors in the sign of rating
        satcue{1,2}(:,14) = mat2cell(cellfun(@(x) -x,satcue{1,2}(:,14)),ones(1,size(satcue{1,2},1)),1);
        satcue{1,5}(:,14) = mat2cell(cellfun(@(x) -x,satcue{1,5}(:,14)),ones(1,size(satcue{1,5},1)),1);
        % Fix ranking (inverted Likert in rating trials):
        ranking{1,2}(:,3) = 3-ranking{1,2}(:,3);
        ranking{1,5}(:,3) = 3-ranking{1,5}(:,3);
        ranking{1,2}(:,4) = -ranking{1,2}(:,4);
        ranking{1,5}(:,4) = -ranking{1,5}(:,4);
        
        
        %% analysis using rating
        
        
        
        %paramters
        for j=3
            
            window_size=j;
            starting_point=window_size + 1;
            
            if starting_point <= window_size
                error('window size has to be smaller than the starting point of the loop!!')
            end
            
            
            shift_starting=starting_point - 1;
            window_bef=zeros(size(satcue{1,2}(:,14),1)-starting_point,1);
            rating_bef=zeros(size(satcue{1,2}(:,14),1)-starting_point,1);
            window_aft=zeros(size(satcue{1,5}(:,14),1)-starting_point,1);
            rating_aft=zeros(size(satcue{1,5}(:,14),1)-starting_point,1);
            
            for i=starting_point:size(satcue{1,2}(:,14),1) %skip first 3 because no previous window to compare
                
                window_bef(i-shift_starting)= mean(cell2mat(satcue{2}(i-window_size:i-1,14)));
                rating_bef(i-shift_starting)= cell2mat(satcue{2}(i,14));
                
            end
            [rho_bef(j),pvalrho_bef(j)]=corr(window_bef, rating_bef, 'type', 'Spearman');
            [tau_bef(j),pvaltau_bef(j)]=corr(window_bef, rating_bef, 'type', 'Kendall');
            
            
            
            
            for i=starting_point:size(satcue{1,5}(:,14),1) %skip first 3 because no previous window to compare
                
                window_aft(i-shift_starting)= mean(cell2mat(satcue{5}(i-window_size:i-1,14)));
                rating_aft(i-shift_starting)= cell2mat(satcue{5}(i,14));
                
            end
            [rho_aft(j),pvalrho_aft(j)]=corr(window_aft, rating_aft, 'type', 'Spearman');
            [tau_aft(j),pvaltau_aft(j)]=corr(window_aft, rating_aft, 'type', 'Kendall');
        end
        
        figure;
        plot((1:length(pvalrho_bef)),pvalrho_bef, 'k-', (1:length(pvalrho_aft)), pvalrho_aft , 'b-', (1:length(pvalrho_aft)), ones(length(pvalrho_aft),1)*0.05 , 'r-', ones(length(pvalrho_aft),1)*0.1 , 'r-');
        % figure;
        % plot((1:length(pvaltau_bef)),pvaltau_bef, 'k-', (1:length(pvaltau_aft)), pvaltau_aft, 'b-');
        
        
        
        
        %% trying stuff
        for j=3
            window_size=j;
            window_bef2=zeros(size(satcue{1,2}(:,14),1)-1,1);
            rating_bef2=zeros(size(satcue{1,2}(:,14),1)-1,1);
            window_aft2=zeros(size(satcue{1,5}(:,14),1)-1,1);
            rating_aft2=zeros(size(satcue{1,5}(:,14),1)-1,1);
            
            for i=1:size(satcue{1,2}(:,14),1) %skip first 3 because no previous window to compare
                
                if  i == 1
                    
                elseif i <= window_size
                    window_bef2(i-1)= mean(cell2mat(satcue{2}(1:i-1,14)));
                    rating_bef2(i-1)= cell2mat(satcue{2}(i,14));
                    
                elseif i > window_size
                    window_bef2(i-1)= mean(cell2mat(satcue{2}(i-window_size:i-1,14)));
                    rating_bef2(i-1)= cell2mat(satcue{2}(i,14));
                end
                
                
                
                
                if  i == 1
                    
                elseif i <= window_size
                    window_aft2(i-1)= mean(cell2mat(satcue{5}(1:i-1,14)));
                    rating_aft2(i-1)= cell2mat(satcue{5}(i,14));
                    
                elseif i > window_size
                    window_aft2(i-1)= mean(cell2mat(satcue{5}(i-window_size:i-1,14)));
                    rating_aft2(i-1)= cell2mat(satcue{5}(i,14));
                end
                
            end
            
            [rho_bef2(j),pvalrho_bef2(j)]=corr(window_bef2, rating_bef2, 'type', 'Spearman');
            [tau_bef2(j),pvaltau_bef2(j)]=corr(window_bef2, rating_bef2, 'type', 'Kendall');
            
            
            [rho_aft2(j),pvalrho_aft2(j)]=corr(window_aft2, rating_aft2, 'type', 'Spearman');
            [tau_aft2(j),pvaltau_aft2(j)]=corr(window_aft2, rating_aft2, 'type', 'Kendall');
            
            %         [rho_bef3(j),pvalrho_bef3(j)]=corr(window_bef2/max(window_bef2), rating_bef2/max(rating_bef2), 'type', 'Spearman');
            %         [tau_bef3(j),pvaltau_bef3(j)]=corr(window_bef2/max(window_bef2), rating_bef2/max(rating_bef2), 'type', 'Kendall');
            %
            %
            %         [rho_aft3(j),pvalrho_aft3(j)]=corr(window_aft2/max(window_aft2), rating_aft2/max(rating_bef2), 'type', 'Spearman');
            %         [tau_aft3(j),pvaltau_aft3(j)]=corr(window_aft2/max(window_aft2), rating_aft2/max(rating_bef2), 'type', 'Kendall');
            %
        end
        
        figure;
        plot((1:length(pvalrho_bef2)),pvalrho_bef2, 'k-', (1:length(pvalrho_aft2)), pvalrho_aft2 , 'b-',(1:length(pvalrho_aft)), ones(length(pvalrho_aft),1)*0.05 , 'r-', ones(length(pvalrho_aft),1)*0.1 , 'r-');
        % figure;
        % plot((1:length(pvalrho_bef3)),pvalrho_bef3, 'k-', (1:length(pvalrho_aft3)), pvalrho_aft3 , 'b-');
        
        
        
        for j=3
            window_size=j;
            window_bef4=zeros(size(satcue{1,2}(:,14),1)-1,1);
            rating_bef4=zeros(size(satcue{1,2}(:,14),1)-1,1);
            window_aft4=zeros(size(satcue{1,5}(:,14),1)-1,1);
            rating_aft4=zeros(size(satcue{1,5}(:,14),1)-1,1);
            
            for i=1:size(satcue{1,2}(:,14),1) %skip first 3 because no previous window to compare
                
                if  i == 1
                    window_bef4(i)= 0;
                    rating_bef4(i)= cell2mat(satcue{2}(i,14));
                elseif i <= window_size
                    window_bef4(i)= (cell2mat(satcue{2}(1:i-1,14)))'*(1:i-1)'/sum(1:i-1);
                    rating_bef4(i)= cell2mat(satcue{2}(i,14));
                    
                elseif i > window_size
                    window_bef4(i)= (cell2mat(satcue{2}(i-window_size:i-1,14)))'*(1:window_size)'/sum(1:window_size);
                    rating_bef4(i)= cell2mat(satcue{2}(i,14));
                end
                
                
                
                
                if  i == 1
                    window_aft4(i)= 0;
                    rating_aft4(i)= cell2mat(satcue{5}(i,14));
                elseif i <= window_size
                    window_aft4(i)= (cell2mat(satcue{5}(1:i-1,14)))'*(1:i-1)'/sum(1:i-1);
                    rating_aft4(i)= cell2mat(satcue{5}(i,14));
                    
                elseif i > window_size
                    window_aft4(i)= (cell2mat(satcue{5}(i-window_size:i-1,14)))'*(1:window_size)'/sum(1:window_size);
                    rating_aft4(i)= cell2mat(satcue{5}(i,14));
                end
                
            end
            
            [rho_bef4(j),pvalrho_bef4(j)]=corr(window_bef4, rating_bef4, 'type', 'Spearman');
            [tau_bef4(j),pvaltau_bef4(j)]=corr(window_bef4, rating_bef4, 'type', 'Kendall');
            
            
            [rho_aft4(j),pvalrho_aft4(j)]=corr(window_aft4, rating_aft4, 'type', 'Spearman');
            [tau_aft4(j),pvaltau_aft4(j)]=corr(window_aft4, rating_aft4, 'type', 'Kendall');
            
            %         [rho_bef3(j),pvalrho_bef3(j)]=corr(window_bef2/max(window_bef2), rating_bef2/max(rating_bef2), 'type', 'Spearman');
            %         [tau_bef3(j),pvaltau_bef3(j)]=corr(window_bef2/max(window_bef2), rating_bef2/max(rating_bef2), 'type', 'Kendall');
            %
            %
            %         [rho_aft3(j),pvalrho_aft3(j)]=corr(window_aft2/max(window_aft2), rating_aft2/max(rating_bef2), 'type', 'Spearman');
            %         [tau_aft3(j),pvaltau_aft3(j)]=corr(window_aft2/max(window_aft2), rating_aft2/max(rating_bef2), 'type', 'Kendall');
            %
        end
        
        figure;
        plot((1:length(pvalrho_bef4)),pvalrho_bef4, 'k-', (1:length(pvalrho_aft4)), pvalrho_aft4 , 'b-', (1:length(pvalrho_aft)), ones(length(pvalrho_aft),1)*0.05 , 'r-', ones(length(pvalrho_aft),1)*0.1 , 'r-');
        
        
        
        
        jj=3;
        val_rat_rat=nan(length(satcue{jj-1}(:,14)),1);
        val_rat_2AFC=nan(length(satcue{jj-1}(:,14)),1);
        val_rat_rank=nan(length(satcue{jj-1}(:,14)),1);
        
        for i=1:20
            ind3=find(cell2mat(satcue{2}(:,2))==i);
            valRat=  ranking{jj-1}(i,4); % ranking
            val2AFC= ranking{jj}(i,4);
            valRank= ranking{jj}(i,3);
            val_rat_rat(ind3)= valRat;
            val_rat_2AFC(ind3)= val2AFC;
            val_rat_rank(ind3)= valRank;
            
        end
        
        
        
        
        
        
        
        
  % Arima models

        autocorr(cell2mat(satcue{2}(:,14)), 20, 0, 1.645);  %accept 90% CI as it is one tailed, it  should be NEGATIVE
        figure;
        parcorr(cell2mat(satcue{2}(:,14)),20, 0, 1.645)
        [cs,csc, csf]= adftest(cell2mat(satcue{2}(:,14)), 'lags', (0:5));      % test for stationarty. If both test indicate stationar, I could use this model for an arima, if non-stationary --> needs transformations of either mean (differencing) or variance (logarithmic transform)
        [cs0,csc0, csf0]= kpsstest(cell2mat(satcue{2}(:,14)), 'lags', (0:5)); 
        
        satdiff= diff(cell2mat(satcue{2}(:,14)));                              % Be aware that differeciating a random iid dataset could create artifact-autocorrelations, mainly negative ones, which is what you are looking for, so, if a new autocorr appears in differencing, could be an artifact
        [cs1,csc1, csf1]= adftest(satdiff, 'lags', (0:5));                     % Even though it could be an artifact (more even if negative correlation), I want to test for higher-order adaptations, such as 
        [cs2,csc2, csf2]= kpsstest(satdiff, 'lags', (0:5));
        figure;
        autocorr(satdiff, 20, 0, 1.645)
        figure;
        parcorr(satdiff,20, 0, 1.645)
        
        satdiff2= diff(cell2mat(satcue{2}(:,14)),2);
        [cs3,csc3, csf3]= adftest(cell2mat(satcue{2}(:,14)), 'lags', (0:5));
        [cs4,csc4, csf4]= kpsstest(cell2mat(satcue{2}(:,14)), 'lags', (0:5));
        figure;
        autocorr(satdiff2,20, 0, 1.645)
        figure;
        parcorr(satdiff2,20, 0, 1.645)
        
        %it looks like differentiating increases the time-dependecy, so the
        % immediate differencies have more autocorrelation. This makes sense as
        % it is the difference what should be changed
        
        
        % provin iid of data
        [h_jb, p_jb] = jbtest(cell2mat(satcue{2}(:,14)));  % test for normality of the values, this would show simmetry, so no skewness towards any value, but also kurtosis==3, but I am nto sure it need be gaussian to prove iid 
        
        
        y= rating_bef4;
        X=[ones(size(y)), window_bef4, (val_rat_rat./3)];  % I scale down the val_rat_rat to match the scale of the response, so male it the mean from -100 to 100 and not -300 to 300, cahnging the scale gives proper betas, but does NOT cahnge the significance of the t-test for the coeffifcient nor te R-squared
        [a1,b1,c1,d1,e1] = regress(y,X);
        [a2,b2,c2,d2,e2] = stepwisefit(X(2:3), y);         % Actually,
        
        
        mdl1= arima(1,0,1);
        mdl2= arima(2,0,2);
        mdl3= arima(1,0,0);
        mdl4= arima(0,0,1);
        mdl5= arima(0,1,1); % I take on on MA as PAFC tails off and AFC cuts-off after lag one
        mdl6= arima(0,2,1); %second order difference
        [est1]= estimate(mdl1, cell2mat(satcue{2}(:,14)));
        [est2]= estimate(mdl2, cell2mat(satcue{2}(:,14)));
        [est3]= estimate(mdl3, cell2mat(satcue{2}(:,14)));
        [est4]= estimate(mdl4, cell2mat(satcue{2}(:,14)));
        [est5]= estimate(mdl5, cell2mat(satcue{2}(:,14)));
        [est6]= estimate(mdl6, cell2mat(satcue{2}(:,14)));
        
        
        % once every model is constructed, test for validity. Necessary
        % that residual show no autocorr, ARCH test negative (for
        % heteroscedasticity, othrewise logtransfomration), exclude autocorrelation with ljung-Box Q
        % test for residuals, check if normal with Kolgomorov-Smirnov test.
        % mayby I could use jbtest for normality or even Bia-Ng tests for 
        % normlaity/skewness (2005, journal of business and economical statistics)
        % Of course only if t-tests coefficiente of the model are correct
        % those who survive, compare with BIC or AIC. (bayesian inference
        % methods), the lowest BIC/AIC the better. 

        
        
end




% % 
% % % this analysis thies to test if the value of the products affects
% % % subsewuent values, such as value adaptation in context dependent Decision
% % % making. Thus I make a window of rating/rankings, and I expect to see that
% % 
% % 
% % dbstop if error
% % [path]= uigetdir('/media/Projects/Alex/Satiation Analysis', 'Select a patient folder');
% % cd (path)
% % %load behavioural results 
% % pathfolders=dir;
% % for i=1:length(pathfolders)
% %     if regexp(pathfolders(i).name, 'Files')
% %        cd (pathfolders(i).name)
% %        filesfolder=dir;
% %        for j=1:length(filesfolder)
% %          if regexp(filesfolder(j).name, 'runsat_output')
% %              load (filesfolder(j).name)
% %          end
% %        end
% %     end
% %     
% % end
% % 
% % %cd (path)
% % 
% % %correct the errors in the sign of rating
% %     satcue{1,2}(:,14) = mat2cell(cellfun(@(x) -x,satcue{1,2}(:,14)),ones(1,size(satcue{1,2},1)),1);
% %     satcue{1,5}(:,14) = mat2cell(cellfun(@(x) -x,satcue{1,5}(:,14)),ones(1,size(satcue{1,5},1)),1);
% %     % Fix ranking (inverted Likert in rating trials):
% %     ranking{1,2}(:,3) = 3-ranking{1,2}(:,3);
% %     ranking{1,5}(:,3) = 3-ranking{1,5}(:,3);
% %     ranking{1,2}(:,4) = -ranking{1,2}(:,4);
% %     ranking{1,5}(:,4) = -ranking{1,5}(:,4);
% %     
% % 
% % %% Analysis using rating 
% % figure;
% % autocorr(cell2mat(satcue{2}(:,14)),20, 0, 1.645)  % plot autocorr with 90% CI, as it should be negatively correlated, I can can accept 90% as if it where one-tail test
% % 
% % [adf_h, adf_p]= adftest(cell2mat(satcue{2}(:,14)), 'lags', [0:5]);   %this is the unit root test, if negative, means non-stationarity, if positive, then means either non stationarty or a trend
% % %parameters
% % for j=3
% %     
% %     window_size=j;
% %     starting_point=window_size + 1;
% % 
% %     if starting_point <= window_size
% %         error('window size has to be smaller than the starting point of the loop!!')
% %     end
% % 
% %     
% %     shift_starting=starting_point - 1;
% %     window_bef=zeros(size(satcue{1,2}(:,14),1)-starting_point,1);
% %     rating_bef=zeros(size(satcue{1,2}(:,14),1)-starting_point,1);
% %     window_aft=zeros(size(satcue{1,5}(:,14),1)-starting_point,1);
% %     rating_aft=zeros(size(satcue{1,5}(:,14),1)-starting_point,1);
% % 
% %     for i=starting_point:size(satcue{1,2}(:,14),1) %skip first 3 because no previous window to compare
% % 
% %         window_bef(i-shift_starting)= mean(cell2mat(satcue{2}(i-window_size:i-1,14)));
% %         rating_bef(i-shift_starting)= cell2mat(satcue{2}(i,14));
% % 
% %     end 
% %         [rho_bef(j),pvalrho_bef(j)]=corr(window_bef, rating_bef, 'type', 'Spearman');
% %         [tau_bef(j),pvaltau_bef(j)]=corr(window_bef, rating_bef, 'type', 'Kendall');
% %         [r_bef(j),pvalr_bef(j)]=corr(window_bef, rating_bef, 'type', 'Pearson');
% % 
% % 
% %  
% % 
% %     for i=starting_point:size(satcue{1,5}(:,14),1) %skip first 3 because no previous window to compare
% % 
% %         window_aft(i-shift_starting)= mean(cell2mat(satcue{5}(i-window_size:i-1,14)));
% %         rating_aft(i-shift_starting)= cell2mat(satcue{5}(i,14));
% % 
% %     end 
% %         [rho_aft(j),pvalrho_aft(j)]=corr(window_aft, rating_aft, 'type', 'Spearman');
% %         [tau_aft(j),pvaltau_aft(j)]=corr(window_aft, rating_aft, 'type', 'Kendall');
% %         [r_aft(j),pvalr_aft(j)]    =corr(window_aft, rating_aft, 'type', 'Pearson');
% % end
% % % figure;
% % % plot((1:length(pvalrho_bef)),pvalrho_bef, 'k-', (1:length(pvaltau_bef)), pvaltau_bef, 'b-', (1:length(pvalr_bef)), pvalr_bef, '-r');
% % % figure;
% % % plot((1:length(pvalrho_aft)),pvalrho_aft, 'k-', (1:length(pvaltau_aft)), pvaltau_aft, 'b-', (1:length(pvalr_aft)), pvalr_aft, '-r');
% % 
% % 
% % % ideas, acf function, to check autocorrelation of the signal over time,
% % % autoregressive models, or moving-average models, that make regression
% % % models for the time-series, depending on previous points. 
% % % also remember that a regression model usually needs iid error (independent)
% % % if the noise depends on time, then autoression allows to control for this
% % % and have better regressions. 
% % 
% % 
% % % 
% % %     for i=1:size(satcue{1,2}(:,14),1) %skip first 3 because no previous window to compare
% % % 
% % %         if i<window_size
% % %         window_bef(i-shift_starting)= mean(cell2mat(satcue{2}(i-window_size:i-1,14)));
% % %         rating_bef(i-shift_starting)= cell2mat(satcue{2}(i,14));
% % %         elseif i==window_size
% % %         window_bef(i-shift_starting)= mean(cell2mat(satcue{2}(i-window_size:i-1,14)));
% % %         rating_bef(i-shift_starting)= cell2mat(satcue{2}(i,14));
% % %         elseif i>window_size       
% % %         window_bef(i-shift_starting)= mean(cell2mat(satcue{2}(i-window_size:i-1,14)));
% % %         rating_bef(i-shift_starting)= cell2mat(satcue{2}(i,14));
% % %         end
% % %         
% % %     end 
% % %         [rho_bef(j),pvalrho_bef(j)]=corr(window_bef, rating_bef, 'type', 'Spearman');
% % %         [tau_bef(j),pvaltau_bef(j)]=corr(window_bef, rating_bef, 'type', 'Kendall');
% % %         [r_bef(j),pvalr_bef(j)]=corr(window_bef, rating_bef, 'type', 'Pearson');
% % 
% % 
% % 
% % %regression analysis
% %     jj=3;
% %     val_rat_rat=nan(length(satcue{jj-1}(:,14)),1);
% %     val_rat_2AFC=nan(length(satcue{jj-1}(:,14)),1);
% %     val_rat_rank=nan(length(satcue{jj-1}(:,14)),1);
% %     
% %     for i=1:20
% %         ind3=find(cell2mat(satcue{2}(:,2))==i);
% %         valRat=  ranking{jj-1}(i,4); % ranking
% %         val2AFC= ranking{jj}(i,4);
% %         valRank= ranking{jj}(i,3);
% %         val_rat_rat(ind3)= valRat;
% %         val_rat_2AFC(ind3)= val2AFC;
% %         val_rat_rank(ind3)= valRank;
% %         
% %     end
% %     
