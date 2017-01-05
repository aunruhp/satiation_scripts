function value_adaptation_ana2




% this analysis thies to test if the value of the products affects
% subsewuent values, such as value adaptation in context dependent Decision
% making. Thus I make a window of rating/rankings, and I expect to see that



[satcue, ranking] = load_behavioural_data;



% extract information about covariates, ie, mean values of stimulus over time to correlate.



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

%scaling values
val_rat_rat= val_rat_rat./3;             % so it goes from -100 to +100
val_rat_2AFC= (val_rat_2AFC ./9.5) -100; % so it goes from -100 to +100


%% select which value use
val= val_rat_2AFC;
% val= val_rat_rat;

a1=cell(1,5);b1=cell(1,5);c1=cell(1,5);d1=cell(1,5);e1=cell(1,5);
a2=cell(1,5);b2=cell(1,5);c2=cell(1,5);d2=cell(1,5);e2=cell(1,5);
a3=cell(1,5);b3=cell(1,5);c3=cell(1,5);d3=cell(1,5);e3=cell(1,5);
a4=cell(1,5);b4=cell(1,5);c4=cell(1,5);d4=cell(1,5);e4=cell(1,5);

for j=1:5
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
    
    
    
    % linear regression models
    
    
    y= rating_bef4;
    X=[ones(size(y)), val(:,1), window_bef4];  % I scale down the val_rat_rat to match the scale of the response, so male it the mean from -100 to 100 and not -300 to 300, cahnging the scale gives proper betas, but does NOT cahnge the significance of the t-test for the coeffifcient nor te R-squared
    y1= rating_aft4;
    X1=[ones(size(y1)), val(:,2), window_aft4];
    
    [a1{j},b1{j},c1{j},d1{j},e1{j}] = regress(y,X);
    [a2{j},b2{j},c2{j},d2{j},e2{j}] = stepwisefit(X(:,2:3), y);
    [a3{j},b3{j},c3{j},d3{j},e3{j}] = regress(y1,X1);
    [a4{j},b4{j},c4{j},d4{j},e4{j}] = stepwisefit(X1(:,2:3), y1);
    
    
    
end

figure;
plot((1:length(pvalrho_bef4)),pvalrho_bef4, 'k-', (1:length(pvalrho_aft4)), pvalrho_aft4 , 'b-', (1:length(pvalrho_aft4)), ones(length(pvalrho_aft4),1)*0.05 , 'r-', (1:length(pvalrho_aft4)), ones(length(pvalrho_aft4),1)*0.1 , 'r-');



%%  %% Regression analysis



for jj=[2,5]
   
    
    loc_var=floor(jj/2);  % floor(jj/2) means if jj=2, 1, if jj=5, then 2
    Y= cell2mat(satcue{jj}(:,14));
    X=[ones(size(val(:,loc_var))), val(:,loc_var)];
    
    % check how the lin reg behaves, check residuals
    
    % prove normal gaussian data
    [h_jb, p_jb] = jbtest(Y);  % test for normality of the values, this would show simmetry, so no skewness towards any value, but also kurtosis==3, but I am nto sure it need to be gaussian to prove iid
    [h_ks, p_ks] = kstest(Y);
    
    figure;
    plot(1:length(Y), Y, 'b-', 1:length(X(:,2)), X(:,2), 'k-')
    [a,b,c,d,e] = regress(Y, X );
    h_arch =archtest(c, 'lags', (1:10)) %should show homodesticity
    h_lqb = lbqtest(c, 'lags', (1:10)) % should show no autocorrelations in residuals
    figure; plot(c);
    figure; autocorr(c ,20, 0, 1.645); %no autocorr
    figure; parcorr(c ,20, 0, 1.645); %no autocorr
    figure; qqplot(c);  %should show normality
%     figure; plot(1:length(Y), Y, 'b-', length(X(:,2)), X(:,2), 'k-')
%     
    
    
    
    % Arima models   %remember, variance and dependent variable must be
    % homoscedastic, but independent varaiable can be heteroscedastic.
    figure;
    autocorr(Y, 20, 0, 1.645);  %accept 90% CI as it is one tailed, it  should be NEGATIVE
    figure;
    parcorr(Y,20, 0, 1.645)
    [asdf,csc, csf]= adftest(Y, 'lags', (0:5))      % test for unit-root vs stationarty. If both test indicate stationary, I could use this to model an arima, if non-stationary --> needs transformations of either mean (differencing) or variance (logarithmic transform)
    [cs0,csc0, csf0]= kpsstest(Y, 'lags', (0:5))  % test if it is unit-root vs trend-stationary. Cobined with unit-root gives more info, as adf test could give unit root, adn then find it was a trend-statioraity, so a drift term
    
%     satdiff= diff(Y);                              % Be aware that differeciating a random iid dataset could create artifact-autocorrelations, mainly negative ones, which is what you are looking for, so, if a new autocorr appears in differencing, could be an artifact
%     [cs1,csc1, csf1]= adftest(satdiff, 'lags', (0:5));                     % Even though it could be an artifact (more even if negative correlation), I want to test for higher-order adaptations, such as
%     [cs2,csc2, csf2]= kpsstest(satdiff, 'lags', (0:5));
%     figure;
%     autocorr(satdiff, 20, 0, 1.645)
%     figure;
%     parcorr(satdiff,20, 0, 1.645)
%     
%     satdiff2= diff(Y,2);
%     [cs3,csc3, csf3]= adftest(satdiff2, 'lags', (0:5));  % unit root test, check if for a AR model, so with AR(x) the coefficient of the AR terms  are ALL less than one significantly, thus statinary as tend t odecrease error effect over time, if rho==1, the coefficent rho is the coefficient of AR, then it is a random walk model, thus non-stationary. If enegative, no non-stationarity in any AR(x) model, so it test for all AR models
%     [cs4,csc4, csf4]= kpsstest(satdiff2, 'lags', (0:5));
%     figure;
%     autocorr(satdiff2,20, 0, 1.645)
%     figure;
%     parcorr(satdiff2,20, 0, 1.645)
%     
%     [cs5,csc5, csf5]= adftest(X(:,2), 'lags', (0:5));
%     [cs6,csc6, csf6]= kpsstest(X(:,2), 'lags', (0:5));
    %it looks like differentiating increases the time-dependecy, so the
    % immediate differencies have more autocorrelation. This makes sense as
    % it is the difference what should be changed
    
    
    
    
    
    % Arima models.
    % this models, assume that the inovation proces, so the inputs (the
    % stimulus in out case), are iid, with a normal(0, sigma²), thus
    % the values of the products (the mean values), should follow a
    % normal dist or t-student dist, with 0 mean. This is not the case
    % so arima model, without covariates, breaks the assumptions, thus
    % not reliable
    
    try
        mdl1= arima(1,0,1);  % becareful with this one, white noise can be modelled as an ARMA (1,0,1), with coefficients -0.5 for AR parameter and +0.5 for MA parameter. Or close to it. As our PAFC/AFC not clear this could fit a wrong model here.
        mdl2= arima(2,0,2);
        mdl3= arima(3,0,3);
        mdl4= arima(1,0,0);
        mdl5= arima(2,0,0);
        mdl6= arima(0,0,1);
        mdl7= arima(0,0,2);
        mdl8= arima(0,0,3); % I take on on MA as PAFC tails off and AFC cuts-off after lag one
        mdl9= arima(0,1,1); %second order difference
        mdl10=arima(0,2,1);
%         [est1, covar1, log1, info1]= estimate(mdl1, Y);
%         [est2, covar2, log2, info2]= estimate(mdl2, Y);
%         [est3, covar3, log3, info3]= estimate(mdl3, Y);
        [est4, covar4, log4, info4]= estimate(mdl4, Y);
        [est5, covar5, log5, info5]= estimate(mdl5, Y);
        [est6, covar6, log6, info6]= estimate(mdl6, Y);
        [est7, covar7, log7, info7]= estimate(mdl7, Y);
        [est8, covar8, log8, info8]= estimate(mdl8, Y);
%         [est9, covar9, log9, info9]= estimate(mdl9, Y);  %check variance explained. If it increases with differencing order, means over-differencing
%         [est10, covar10, log10, info10]= estimate(mdl10, Y);  %check variance explained. If it increases with differencing order, means over-differencing
%         
        
        % once every model is constructed, test for validity. Necessary
        % that residual show no autocorr check ACF anc PACF should be all
        % negative,  ARCH test negative (for
        % heteroscedasticity, othrewise logtransfomration, random noise has constant variance
        % if not consant more parameters AR/MA to correct), exclude autocorrelation with ljung-Box Q
        % test for residuals, check if normal with Kolgomorov-Smirnov test.
        % mayby I could use jbtest for normality or even Bia-Ng tests for
        % normlaity/skewness (2005, journal of business and economical statistics)
        % Of course only if t-tests coefficiente of the model are correct
        % those who survive, compare with BIC or AIC. (bayesian inference
        % methods), the lowest BIC/AIC the better.
        
        % By default, the time series errors (also called unconditional
        % disturbances) are independent, identically distributed, mean 0
        % Gaussian random variables. If not combine Regression with Arima
        % for residuals/errors
        
        
        %arimax models, estimation of betas and arima at the same time, so
        %arima is not about residuals but about the signal
        
        y0=[0;0;0];  % for arimax  models, you need to have len(dependent)+ order of arima: this implies to add a initial value for each regressor
        % thus, if I(0), one y0 val, if I(2), 3 values. I
        % assume that the value before experiment should be
        % zero, there is always some time between experimetns
        % that allows for this assumption.
        
        % %          X(:,2)= (val_rat_rat(:,floor(jj/2))./3); % which covariate should be included
        
%         [est1x, covar1x, log1x, info1x]= estimate(mdl1, Y, 'X', X(:,2),'Display','params', 'Y0', y0);
%         [est2x, covar2x, log2x, info2x]= estimate(mdl2, Y, 'X', X(:,2),'Display','params', 'Y0', y0);
%         [est3x, covar3x, log3x, info3x]= estimate(mdl3, Y, 'X', X(:,2),'Display','params', 'Y0', y0);
        [est4x, covar4x, log4x, info4x]= estimate(mdl4, Y, 'X', X(:,2),'Display','params', 'Y0', y0);
        [est5x, covar5x, log5x, info5x]= estimate(mdl5, Y, 'X', X(:,2),'Display','params', 'Y0', y0);
        [est6x, covar6x, log6x, info6x]= estimate(mdl6, Y, 'X', X(:,2),'Display','params', 'Y0', y0);
        [est7x, covar7x, log7x, info7x]= estimate(mdl7, Y, 'X', X(:,2),'Display','params', 'Y0', y0);
        [est8x, covar8x, log8x, info8x]= estimate(mdl8, Y, 'X', X(:,2),'Display','params', 'Y0', y0);
%         [est9x, covar9x, log9x, info9x]= estimate(mdl9, Y, 'X', X(:,2),'Display','params', 'Y0', y0);
%         [est10x, covar10x, log10x, info10x]= estimate(mdl10, Y, 'X', X(:,2),'Display','params', 'Y0', y0);
%         
        %regArima models, assumes gaussian residuals, performs lin reg to
        %estimate betas, and then uses residals to estimate arimas, so only
        %residuals are modelled as time-dependent.
        
        mdl11 = regARIMA(1,0,1);
        mdl12 = regARIMA(2,0,2);
        mdl13 = regARIMA(3,0,3);
        mdl14= regARIMA(1,0,0);
        mdl15= regARIMA(2,0,0);
        mdl16= regARIMA(0,0,1);
        mdl17= regARIMA(0,0,2);
        mdl18= regARIMA(0,0,3);
        mdl19= regARIMA(0,1,1);
        mdl20= regARIMA(0,2,1);
%         [est11, covar11, log11, info11]= estimate(mdl11, Y, 'X', X(:,2),'Display','params');   % estimate parameters, not included from linear regression but estimated here
%         [est12, covar12, log12, info12]= estimate(mdl12, Y, 'X', X(:,2),'Display','params');
%         [est13, covar13, log13, info13]= estimate(mdl13, Y, 'X', X(:,2), 'Display','params');   %  intercept not included as varable... be careful, maybe necessary
        [est14, covar14, log14, info14]= estimate(mdl14, Y, 'X', X(:,2), 'Display','params');
        [est15, covar15, log15, info15]= estimate(mdl15, Y, 'X', X(:,2), 'Display','params');
        [est16, covar16, log16, info16]= estimate(mdl16, Y, 'X', X(:,2), 'Display','params');
        [est17, covar17, log17, info17]= estimate(mdl17, Y, 'X', X(:,2), 'Display','params');   %  intercept not included as varable... be careful, maybe necessary
        [est18, covar18, log18, info18]= estimate(mdl18, Y, 'X', X(:,2), 'Display','params');
%         [est19, covar19, log19, info19]= estimate(mdl19, Y, 'X', X(:,2), 'Display','params');
%         [est20, covar20, log20, info20]= estimate(mdl20, Y, 'X', X(:,2), 'Display','params');
    catch
        warning('error somewhere in your models. Bro!')
    end
    
    % be careful while modelling, the innovation process(error) must
    % be normal and have mean 0.
    
    
    %Properties of a MA process or a AR process : is stationary  AND
    %weakly coupled. If not, differenciate it until is stationary.
    %Remember that weakly dependent means that the correlation of x(t)
    %wit x(t-h) tends to zero as h tends to infinity
    



end



end




%% analysis using rating


%
%         %paramters
%         for j=1:5
%
%             window_size=j;
%             starting_point=window_size + 1;
%
%             if starting_point <= window_size
%                 error('window size has to be smaller than the starting point of the loop!!')
%             end
%
%
%             shift_starting=starting_point - 1;
%             window_bef=zeros(size(satcue{1,2}(:,14),1)-starting_point,1);
%             rating_bef=zeros(size(satcue{1,2}(:,14),1)-starting_point,1);
%             window_aft=zeros(size(satcue{1,5}(:,14),1)-starting_point,1);
%             rating_aft=zeros(size(satcue{1,5}(:,14),1)-starting_point,1);
%
%             for i=starting_point:size(satcue{1,2}(:,14),1) %skip first 3 because no previous window to compare
%
%                 window_bef(i-shift_starting)= mean(cell2mat(satcue{2}(i-window_size:i-1,14)));
%                 rating_bef(i-shift_starting)= cell2mat(satcue{2}(i,14));
%
%             end
%             [rho_bef(j),pvalrho_bef(j)]=corr(window_bef, rating_bef, 'type', 'Spearman');
%             [tau_bef(j),pvaltau_bef(j)]=corr(window_bef, rating_bef, 'type', 'Kendall');
%
%
%
%
%             for i=starting_point:size(satcue{1,5}(:,14),1) %skip first 3 because no previous window to compare
%
%                 window_aft(i-shift_starting)= mean(cell2mat(satcue{5}(i-window_size:i-1,14)));
% %                 rating_aft(i-shift_starting)= cell2mat(satcue{5}(i,14));
% %
% %             end
% %             [rho_aft(j),pvalrho_aft(j)]=corr(window_aft, rating_aft, 'type', 'Spearman');
% %             [tau_aft(j),pvaltau_aft(j)]=corr(window_aft, rating_aft, 'type', 'Kendall');
% %         end
% %
% %         figure;
% %         plot((1:length(pvalrho_bef)),pvalrho_bef, 'k-', (1:length(pvalrho_aft)), pvalrho_aft , 'b-', (1:length(pvalrho_aft)), ones(length(pvalrho_aft),1)*0.05 , 'r-',(1:length(pvalrho_aft)), ones(length(pvalrho_aft),1)*0.1 , 'r-');
% %         % figure;
% %         % plot((1:length(pvaltau_bef)),pvaltau_bef, 'k-', (1:length(pvaltau_aft)), pvaltau_aft, 'b-');
%
%
%
%
%         %% trying stuff
%         for j=1:5
%             window_size=j;
%             window_bef2=zeros(size(satcue{1,2}(:,14),1)-1,1);
%             rating_bef2=zeros(size(satcue{1,2}(:,14),1)-1,1);
%             window_aft2=zeros(size(satcue{1,5}(:,14),1)-1,1);
%             rating_aft2=zeros(size(satcue{1,5}(:,14),1)-1,1);
%
%             for i=1:size(satcue{1,2}(:,14),1) %skip first 3 because no previous window to compare
%
%                 if  i == 1
%
%                 elseif i <= window_size
%                     window_bef2(i-1)= mean(cell2mat(satcue{2}(1:i-1,14)));
%                     rating_bef2(i-1)= cell2mat(satcue{2}(i,14));
%
%                 elseif i > window_size
%                     window_bef2(i-1)= mean(cell2mat(satcue{2}(i-window_size:i-1,14)));
%                     rating_bef2(i-1)= cell2mat(satcue{2}(i,14));
%                 end
%
%
%
%
%                 if  i == 1
%
%                 elseif i <= window_size
%                     window_aft2(i-1)= mean(cell2mat(satcue{5}(1:i-1,14)));
%                     rating_aft2(i-1)= cell2mat(satcue{5}(i,14));
%
%                 elseif i > window_size
%                     window_aft2(i-1)= mean(cell2mat(satcue{5}(i-window_size:i-1,14)));
%                     rating_aft2(i-1)= cell2mat(satcue{5}(i,14));
%                 end
%
%             end
%
%             [rho_bef2(j),pvalrho_bef2(j)]=corr(window_bef2, rating_bef2, 'type', 'Spearman');
%             [tau_bef2(j),pvaltau_bef2(j)]=corr(window_bef2, rating_bef2, 'type', 'Kendall');
%
%
%             [rho_aft2(j),pvalrho_aft2(j)]=corr(window_aft2, rating_aft2, 'type', 'Spearman');
%             [tau_aft2(j),pvaltau_aft2(j)]=corr(window_aft2, rating_aft2, 'type', 'Kendall');
%
%
%         end
%
%         figure;
%         plot((1:length(pvalrho_bef2)),pvalrho_bef2, 'k-', (1:length(pvalrho_aft2)), pvalrho_aft2 , 'b-',(1:length(pvalrho_aft)), ones(length(pvalrho_aft),1)*0.05 , 'r-',(1:length(pvalrho_aft)), ones(length(pvalrho_aft),1)*0.1 , 'r-');
%
%

