function [t,p]=ttestbetas(lm, dof,k)

  beta1= lm{1}.Coefficients.Estimate(k);
  beta2= lm{2}.Coefficients.Estimate(k);
  se1= lm{1}.Coefficients.SE(k);
  se2= lm{2}.Coefficients.SE(k);
  
  num= (beta1- beta2);
  den= sqrt(se1^2+se2^2);
  t= num/den;

  p = 2 * tcdf(-abs(t), dof);

  
end