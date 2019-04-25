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

S_iso = (J^(-2/3))*((-1/3)*(mu_1*I1_bar + 2*mu_2*I2_bar)*C_bar_inv + (mu_1 + mu_2*I1_bar)*I - mu_2*C);
S_vol = (J^(1/3))*(J-1)*k*C_bar_inv;

S_Kirchoff = S_iso + S_vol;
