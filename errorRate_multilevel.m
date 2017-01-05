function errorRate_multilevel 



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
[final_coeff{i}, final_corrfreq{i}, final_noncorr_freq{i}, final_corrCoeff{i}, corrDiffError{i}, pvalDiffError{i}, corrInOut{i}, pvalcorrInOut{i}, corrDiffErrorAbs{i}, pvalDiffErrorAbs{i}, corrInOutAbs{i}, pvalcorrInOutAbs{i},corrRTError{i}, pvalcorrRTError{i},corrRTErrorFreq{i}, pvalcorrRTErrorFreq{i}]= error_rate2(f{i});
% error_rate2(f{i});
end


jjj=1;
for jj=[1:7,9:11]
   A1(jjj,:)=final_coeff{jj}{1}; B1(jjj,:)=final_corrfreq{jj}{1}; C1(jjj,:)=final_noncorr_freq{jj}{1}; D1(jjj,:)=final_corrCoeff{jj}{1}; 
   A2(jjj,:)=final_coeff{jj}{2}; B2(jjj,:)=final_corrfreq{jj}{2}; C2(jjj,:)=final_noncorr_freq{jj}{2}; D2(jjj,:)=final_corrCoeff{jj}{2}; 
   jjj=jjj+1;
end

% % beforeA=nansum(A1,1)/10; beforeB=nansum(B1,1)/10; beforeC=nansum(C1,1)/10; beforeD=nansum(D1,1)/10;
% % afterA=nansum(A2,1)/10; afterB=nansum(B2,1)/10; afterC=nansum(C2,1)/10; afterD=nansum(D2,1)/10;
beforeA=nanmean(A1,1); beforeB=nanmean(B1,1); beforeC=nanmean(C1,1); beforeD=nanmean(D1,1);
afterA=nanmean(A2,1); afterB=nanmean(B2,1); afterC=nanmean(C2,1); afterD=nanmean(D2,1);


figure;
bar([beforeA;afterA]);
title('Partial Error Coefficient'); ylabel('Mean Partial Error Coefficient')

figure;
bar([beforeC;afterC]);
title('Partial Error Frequency');ylabel('Mean Partial  Error Frequency' )

figure;
bar([beforeB;afterB]);
title('Partial Error Frequency Corrected');ylabel('Mean Partial Error Frequency' )

figure;
bar([beforeD;afterD]);
title('Partial Error Coefficient Corrected');ylabel('Mean Partial Error Coefficient' )


%% intransitivities analysis 

for i=1:length(f)
[pval_intrans{i}, intransitivity_ratio{i}, diff_intran{i}]=intransitivities_analysis(f{i});
end



% step= 1200/9;
% for jj=[1:7, 9:11]
%     for k=1:2
%         E1(jj,k)=(sum(diff_intran{jj}{k}< (step*1 -600)) + sum(diff_intran{jj}{k}> step*8 -600))/190;
%         E2(jj,k)=(sum(diff_intran{jj}{k}< (step*2 -600) & diff_intran{jj}{k}> step*1 -600) + sum(diff_intran{jj}{k}< (step*8 -600) & diff_intran{jj}{k}> step*7 -600))/190;
%         E3(jj,k)=(sum(diff_intran{jj}{k}< (step*3 -600) & diff_intran{jj}{k}> step*2 -600) + sum(diff_intran{jj}{k}< (step*7 -600) & diff_intran{jj}{k}> step*6 -600))/190;
%         E4(jj,k)=(sum(diff_intran{jj}{k}< (step*4 -600) & diff_intran{jj}{k}> step*3 -600) + sum(diff_intran{jj}{k}< (step*6 -600) & diff_intran{jj}{k}> step*5 -600))/190;
%         E5(jj,k)=(sum(diff_intran{jj}{k}< (step*5 -600) & diff_intran{jj}{k}> step*4 -600))/190;
%     end
% 
% end
% 
% 
% BeforeF=[];AfterF=[];
% for kk=5:-1:1 
%     eval(['BeforeF=[BeforeF, sum(E' num2str(kk) '(:,1))/10];' ])
%     eval(['AfterF=[AfterF, sum(E' num2str(kk) '(:,2))/10];' ])
% end



step= 1200/9;
jjj=1;
for jj=[1:7, 9:11]
    
    for ll= 1:9
        
        E1(jjj,ll)=sum(diff_intran{jj}{1}< (step*ll -600) & diff_intran{jj}{1}> (step*(ll-1) -600)) /190;
        E2(jjj,ll)=sum(diff_intran{jj}{2}< (step*ll -600) & diff_intran{jj}{2}> (step*(ll-1) -600)) /190;
    end
    jjj=jjj+1;
end


beforeE=nanmean(E1,1);
afterE=nanmean(E2,1);

figure;
bar([BeforeE;AfterE]);
title('Absolute Error Frequency');ylabel('Mean Absolute Error Frequency' )

