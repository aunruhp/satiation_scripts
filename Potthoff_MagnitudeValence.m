function [pottlm, compFval, interceptPval]=Potthoff_MagnitudeValence(RTin, valin)


dbstop if error

g=0;
for hh=1:2

    RT2=RTin(:,hh);
    val2=valin(:,hh);  
  

    ind=find(val2>=0);
    valence=zeros(size(val2));
    valence(ind)=1;
    magnitude= abs(val2);




    tbl= table(zscore(RT2), zscore(magnitude), valence, 'VariableNames', {'RT', 'magnitude', 'valence' });

    pottlm{1+g} = fitlm(tbl, 'RT~magnitude');
    pottlm{2+g} = fitlm(tbl, 'RT~magnitude*valence');
    dof1=2;
    dof2=4;


    num= pottlm{2}.SSR-pottlm{1}.SSR;
    MSE= pottlm{2}.SSE/pottlm{2}.DFE;  % sum squares of error / dof of error
    compModel= num/((dof2-dof1) * MSE);
    pvalF = 1-fcdf(compModel,((dof2-dof1)), pottlm{2}.DFE );% F will have an F distribution, with (p2−p1, n−p2) degrees of freedom.
    compFval{hh}= [compModel, pvalF];


    shift=unique(magnitude);
    
    for j=1:length(shift)
        tbl2= table(RT2, magnitude-shift(j), valence, 'VariableNames', {'RT', 'magnitude', 'valence' }); %shift value to have intercept at-100
        potlm2 = fitlm(tbl2, 'RT~magnitude*valence');
        interceptPval{hh}(j)=potlm2.Coefficients.pValue(3); 
    end
    g=g+2;
end



% multivariate analysis
    
    
% tbl= table(RT2, magnitude, valence, group, 'VariableNames', {'RT', 'magnitude', 'valence', 'group'});
% 
% pottlm{1} = fitlm(tbl, 'RT~magnitude');
% pottlm{2} = fitlm(tbl, 'RT~magnitude*valence');
% pottlm{3} = fitlm(tbl, 'RT~magnitude*valence*group');
% pottlm{4} = fitlm(tbl, 'RT~magnitude*valence*group-magnitude:valence:group ');
% dof1=2;
% dof2=4;
% dof3=6;
% 
% 
% num= pottlm{3}.SSR-pottlm{2}.SSR;
% MSE= pottlm{3}.SSE/pottlm{3}.DFE;  % sum squares of error / dof of error
% compModel= num/((dof3-dof2) * MSE);
% pvalF = 1-fcdf(compModel,((dof3-dof2)), pottlm{3}.DFE );% F will have an F distribution, with (p2−p1, n−p2) degrees of freedom.
% compFval= [compModel, pvalF];
% 
% 
% shift=unique(magnitude);
% for j=1:length(shift)
%     tbl2= table(RT2, magnitude-shift(j), valence, 'VariableNames', {'RT', 'magnitude', 'valence' }); %shift value to have intercept at-100
%     potlm2 = fitlm(tbl2, 'RT~magnitude*valence');
%     interceptPval(j)=potlm2.Coefficients.pValue(3); 
% end
% 
% 



