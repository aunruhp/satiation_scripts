function [potlm, compFval, interceptPval]=Potthoff_MagitudeValence(RTin, valin)


dbstop if error



RTtest= ttest2(RTin(:,1), RTin(:,2));
Valtest= ttest2(valin(:,1), valin(:,2));

RT2=[RTin(:,2); RTin(:,2)];
val2=[valin(:,2); valin(:,2)];  
group=[zeros(size(valin(:,1))); ones(size(valin(:,2)))]; % include dummy variable

ind=find(val2>=0);
valence=zeros(size(val2));
valence(ind)=1;
magnitude= abs(val2);


    
    
tbl= table(RT2, magnitude, valence, 'VariableNames', {'RT', 'magnitude', 'valence' });

potlm{1} = fitlm(tbl, 'RT~magnitude');
potlm{2} = fitlm(tbl, 'RT~magnitude*valence');
dof1=2;
dof2=4;


num= potlm{2}.SSR-potlm{1}.SSR;
MSE= potlm{2}.SSE/potlm{2}.DFE;  % sum squares of error / dof of error
compModel= num/((dof2-dof1) * MSE);
pvalF = 1-fcdf(compModel,((dof2-dof1)), potlm{2}.DFE );% F will have an F distribution, with (p2−p1, n−p2) degrees of freedom.
compFval= [compModel, pvalF];


shift=unique(magnitude);
for j=1:length(shift)
    tbl2= table(RT2, magnitude-shift(j), valence, 'VariableNames', {'RT', 'magnitude', 'valence' }); %shift value to have intercept at-100
    potlm2 = fitlm(tbl2, 'RT~magnitude*valence');
    interceptPval(j)=potlm2.Coefficients.pValue(3); 
end


% multivariate analysis
    
    
tbl= table(RT2, magnitude, valence, group, 'VariableNames', {'RT', 'magnitude', 'valence', 'group'});

potlm{1} = fitlm(tbl, 'RT~magnitude');
potlm{2} = fitlm(tbl, 'RT~magnitude*valence');
potlm{3} = fitlm(tbl, 'RT~magnitude*valence*group');
potlm{4} = fitlm(tbl, 'RT~magnitude*valence*group-magnitude:valence:group ');
dof1=2;
dof2=4;
dof3=6;


num= potlm{3}.SSR-potlm{2}.SSR;
MSE= potlm{3}.SSE/potlm{3}.DFE;  % sum squares of error / dof of error
compModel= num/((dof3-dof2) * MSE);
pvalF = 1-fcdf(compModel,((dof3-dof2)), potlm{3}.DFE );% F will have an F distribution, with (p2−p1, n−p2) degrees of freedom.
compFval= [compModel, pvalF];


shift=unique(magnitude);
for j=1:length(shift)
    tbl2= table(RT2, magnitude-shift(j), valence, 'VariableNames', {'RT', 'magnitude', 'valence' }); %shift value to have intercept at-100
    potlm2 = fitlm(tbl2, 'RT~magnitude*valence');
    interceptPval(j)=potlm2.Coefficients.pValue(3); 
end





