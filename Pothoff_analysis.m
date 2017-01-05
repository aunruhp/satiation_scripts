function [potlm, compFval, interceptPval]=Potthoff_analysis(RTin, valin)

% This analysis is thought to compare regrssion slopes in independent
% groups, I'll use it for within slopes before and after lunch, though.
% This is better than comparing betas when the original data are present as it 
% does not gave rounding effects as it uses all data for comparisson
%
% It consists in mixing the data form before and after lunch of value_RT
% Then make a Restricted model, RT = a + beta * Difference
% Then in a second step RT= a + beta1*Difference + beta2* Before_After + beta3* Interaction

% First test if difference in RT between groups, and values
% second build 2 step models
% Third, compare the R-squared of nested model with F-test
% Fourth, beta2 = difference in intercept
% Fifth, beta3= difference in slope
% Finally, shift the intercept, so that 

dbstop if error



RTtest= ttest2(RTin(:,1), RTin(:,2));
Valtest= ttest2(valin(:,1), valin(:,2));

RT=[RTin(:,1); RTin(:,2)];
val=[valin(:,1); valin(:,2)];
group=[zeros(size(valin(:,1))); ones(size(valin(:,2)))]; % include dummy variable


    
    
tbl= table(RT, val, group, 'VariableNames', {'RT', 'val', 'group' });

potlm{1} = fitlm(tbl, 'RT~val');
potlm{2} = fitlm(tbl, 'RT~val*group');
dof1=2;
dof2=4;


num= potlm{2}.SSR-potlm{1}.SSR;
MSE= potlm{2}.SSE/potlm{2}.DFE;  % sum squares of error / dof of error
compModel= num/((dof2-dof1) * MSE);
pvalF=fcdf(compModel,((dof2-dof1)), potlm{2}.DFE );% F will have an F distribution, with (p2−p1, n−p2) degrees of freedom.
compFval(1:2)= [compModel, pvalF];

shift=unique(val);
for j=1:length(shift)
    tbl2= table(RT, val-shift(j), group, 'VariableNames', {'RT', 'val', 'group' }); %shift value to have intercept at-100
    potlm2 = fitlm(tbl2, 'RT~group*val');
    interceptPval(j)=potlm2.Coefficients.pValue(3); 
end

end



