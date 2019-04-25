% Evaluates the Cauchy Stress Tensor for Material 26 %

function Cauchy = stress26(kinematics,properties,dim)

mu_1 = properties(2);
mu_2 = properties(2); % same as mu_1
K = properties(2); % same as mu_1 and mu_2
F = kinematics.F;
J = kinematics.J;
C = transpose(F)*F;
I1 = trace(C);
I2 = (1/2)*(I1^2 - trace(C^2));
I1_bar = (J^(-2/3))*I1;
I2_bar = (J^(-4/3))*I2;
C_bar = (J^(-2/3))*C;
C_bar_inv = inv(C_bar);
I = eye(length(J));
I_double = eye(dim,dim,dim,dim);
c = zeros(dim,dim,dim,dim);
c_vol = zeros(dim,dim,dim,dim);
c_iso = zeros(dim,dim,dim,dim);

for l = 1:dim
  for k = 1:dim
    for j = 1:dim
      for i = 1:dim
        c_vol(i,j,k,l) = (J^(-1/3))*(2*J - 1)*K*(C_bar_inv(i,j)*C_bar_inv(k,l)) - 2*(J^(-1/3))*(J-1)*K*(C_bar_inv(i,l)*C_bar_inv(j,k));
        c_iso(i,j,k,l) = 2*(J^(-4/3))*mu_2*(I(i,j)*I(k,l) - I_double) - (2/3)*(J^(-4/3))*(mu_1 + 2*mu_2*I1_bar)*(C_bar_inv(i,j)*I(k,l) + I(i,j)*C_bar_inv(k,l)) + (4/3)*(J^(-4/3))*mu_2*(C_bar_inv(i,j)*C_bar(k,l) + C_bar(i,j)*C_bar_inv(k,l)) + (2/9)*(J^(-4/3))*(mu_1*I1_bar = 4*mu_2*I2_bar)*(C_inv_bar(i,j)*C_inv_bar(k,l)) + (2/3)*(J^(-4/3))*(mu_1*I1_bar + 2*mu_2*I2_bar)*(C_inv_bar(i,l)*C_inv_bar(j,k));
        
        c(i,j,k,l) = c_vol(i,j,k,l) + c_iso(i,j,k,l);
      end
    end
  end
end
