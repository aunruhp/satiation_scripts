function attention_bias_multilevel



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
    
    [rankRat_bias{i},rankAFC_bias{i},rankTot_bias{i},meanRat_bias{i},meanAFC_bias{i},meanTot_bias{i}]= attention_bias(f{i});
    
end



for k=1:length(f)
    
    for l=1:2
        
        rankRatnum(k,l)= rankRat_bias{k}{l};
        rankAFCnum(k,l)= rankAFC_bias{k}{l};
        rankTotnum(k,l)= rankTot_bias{k}{l};
        meanRatnum(k,l)= meanRat_bias{k}{l};
        meanAFCnum(k,l)= meanAFC_bias{k}{l};
        meanTotnum(k,l)= meanTot_bias{k}{l};
    end
    
end

%ranksum analysis (or ttest2)

for k=1:2
    
    % One-tail Binomial test.  Attention bias should be bigger than chance
    
    bias_num= sum(rankRatnum(:,1));
    p_binoRat = sum(binopdf((bias_num:220), 220, 0.05));
    

    
    [pBinoRat(k)]=myBinomTest(sum(rankRatnum(:,k)), 220, 0.1, 'one');
    [pBinoAFC(k)]=myBinomTest(sum(rankAFCnum(:,k)), 220, 0.1, 'one');
    [pBinoTot(k)]=myBinomTest(sum(rankTotnum(:,k)), 220, 0.1, 'one');
    

        uMLE= mle(sum(rankRatnum(:,k)), 'distribution', 'binomial', 'ntrials', 220);
        
        uLogL=  sum(log(binopdf(sum(rankRatnum(:,k)), 220, uMLE)));
        rLogL=  sum(log(binopdf(sum(rankRatnum(:,k)), 220, 0.1)));
        
        [h_lratioLR, p_lratioRat(k), RatioLR, CriticalValueLR] = lratiotest(uLogL, rLogL,1);  % make the test likelihood ratio with the negative loglikelihood and 1 degree of freedom.
        
        
        if h_lratioLR>0
            disp('This patient has a significat bias in likelihood ratio test')
            
        end

        uMLE= mle(sum(rankAFCnum(:,k)), 'distribution', 'binomial', 'ntrials', 220);
        
        uLogL=  sum(log(binopdf(sum(rankAFCnum(:,k)), 220, uMLE)));
        rLogL=  sum(log(binopdf(sum(rankAFCnum(:,k)), 220, 0.1)));
        
        [h_lratioLR, p_lratioAFC(k), RatioLR, CriticalValueLR] = lratiotest(uLogL, rLogL,1);  % make the test likelihood ratio with the negative loglikelihood and 1 degree of freedom.
        
        
        if h_lratioLR>0
            disp('This patient has a significat bias in likelihood ratio test')
            
        end
  
    
        uMLE= mle(sum(rankTotnum(:,k)), 'distribution', 'binomial', 'ntrials', 220);
        
        uLogL=  sum(log(binopdf(sum(rankTotnum(:,k)), 220, uMLE)));
        rLogL=  sum(log(binopdf(sum(rankTotnum(:,k)), 220, 0.1)));
        
        [h_lratioLR, p_lratioTot(k), RatioLR, CriticalValueLR] = lratiotest(uLogL, rLogL,1);  % make the test likelihood ratio with the negative loglikelihood and 1 degree of freedom.
        
        
        if h_lratioLR>0
            disp('This patient has a significat bias in likelihood ratio test')
            
        end
    
end


% mean analysis

for k=1:2
    
    % One-tail Binomial test.  Attention bias should be bigger than chance
    
    bias_num= sum(meanRatnum(:,k));
    pbinoRatLefttail(k) = sum(binopdf((bias_num:220), 220, 0.5));
    bias_num= sum(meanAFCnum(:,k));
    pbinoAFCLefttail(k) = sum(binopdf((bias_num:220), 220, 0.5));
    bias_num= sum(meanTotnum(:,k));
    pbinoTotLefttail(k) = sum(binopdf((bias_num:220), 220, 0.5));
    

    
    [pBinoRat(k)]=myBinomTest(sum(meanRatnum(:,k)), 220, 0.5, 'one');
    [pBinoAFC(k)]=myBinomTest(sum(meanAFCnum(:,k)), 220, 0.5, 'one');
    [pBinoTot(k)]=myBinomTest(sum(meanTotnum(:,k)), 220, 0.5, 'one');
    

        uMLE= mle(sum(meanRatnum(:,k)), 'distribution', 'binomial', 'ntrials', 220);
        
        uLogL=  sum(log(binopdf(sum(meanRatnum(:,k)), 220, uMLE)));
        rLogL=  sum(log(binopdf(sum(meanRatnum(:,k)), 220, 0.5)));
        
        [h_lratioLR, p_lratioRat(k), RatioLR, CriticalValueLR] = lratiotest(uLogL, rLogL,1);  % make the test likelihood ratio with the negative loglikelihood and 1 degree of freedom.
        
        
        if h_lratioLR>0
            disp('This patient has a significat bias in likelihood ratio test')
            
        end

        uMLE= mle(sum(meanAFCnum(:,k)), 'distribution', 'binomial', 'ntrials', 220);
        
        uLogL=  sum(log(binopdf(sum(meanAFCnum(:,k)), 220, uMLE)));
        rLogL=  sum(log(binopdf(sum(meanAFCnum(:,k)), 220, 0.5)));
        
        [h_lratioLR, p_lratioAFC(k), RatioLR, CriticalValueLR] = lratiotest(uLogL, rLogL,1);  % make the test likelihood ratio with the negative loglikelihood and 1 degree of freedom.
        
        
        if h_lratioLR>0
            disp('This patient has a significat bias in likelihood ratio test')
            
        end
  
    
        uMLE= mle(sum(meanTotnum(:,k)), 'distribution', 'binomial', 'ntrials', 220);
        
        uLogL=  sum(log(binopdf(sum(meanTotnum(:,k)), 220, uMLE)));
        rLogL=  sum(log(binopdf(sum(meanTotnum(:,k)), 220, 0.5)));
        
        [h_lratioLR, p_lratioTot(k), RatioLR, CriticalValueLR] = lratiotest(uLogL, rLogL,1);  % make the test likelihood ratio with the negative loglikelihood and 1 degree of freedom.
        
        
        if h_lratioLR>0
            disp('This patient has a significat bias in likelihood ratio test')
            
        end
    
end






end
