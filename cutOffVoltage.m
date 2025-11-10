function [val,term,dir] = cutOffVoltage(t, u, I, data, spm_matrices, N, vlim)

[V,~] = calculate_volt_temp(t, u, I, data, spm_matrices, N);

Vsparse_dis = V - vlim(1); % Voltage left on discharge
Vsparse_chg = vlim(2) - V; % Voltage left on charge

val = [Vsparse_dis; Vsparse_chg];
term = [1;1];
dir = [-1;-1];

end