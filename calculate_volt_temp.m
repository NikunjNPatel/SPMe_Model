function [V,T] = calculate_volt_temp(t, u, I, cellData, spm_matrices, N)

% Extract anode, cathode state and temperature 
N = spm_matrices.N;
M_a = spm_matrices.M_a;
K_s = spm_matrices.K_s;
Q_c = spm_matrices.Q_c;
us3 = u(1:N-1, :); % Cathode U
us1 = u(N:2*N-2,:); % Anode U
ue = u(2*N-1:(2*N-2+M_a+K_s+Q_c-3),:);
T = u((2*N-2+M_a+K_s+Q_c-3)+1,:); % Temperature
% T = 298.15;

% Calculate applied current density [A.m-2]
i_app = I(t)/cellData.As;

% Calculate molar flux on the electrode [mol.m-2.s-1]
j3 = -i_app/(cellData.as3*cellData.F*cellData.thick3);  % Cathode
j1 = i_app/(cellData.as1*cellData.F*cellData.thick1); % Anode

% Calculate temperature dependent diffusion coefficient
% D = Dref*exp(E/R*(1/Tref - 1/T))
Ds3 = cellData.Ds3_ref*exp((cellData.Ea_Ds3/cellData.R)*(1/cellData.T_ref - 1/T)); % Cathode Diffusion Coefficient
Ds1 = cellData.Ds1_ref*exp((cellData.Ea_Ds1/cellData.R)*(1/cellData.T_ref - 1/T)); % Anode Diffusion Coefficient


% Calculate surface concentration
css3 = spm_matrices.Cp(1,:)*us3 + (spm_matrices.Dp(1,:)*j3)/Ds3;    % Cathode surface concentration
css1 = spm_matrices.Cn(1,:)*us1 + (spm_matrices.Dn(1,:)*j1)/Ds1;    % Anode surface concentration

% k = kref*exp(E/R*(1/Tref - 1/T))
k3 = cellData.k3_ref*exp((cellData.Ea_k3/cellData.R)*(1/cellData.T_ref - 1/T)); % Cathode chemical rate constant
k1 = cellData.k1_ref*exp((cellData.Ea_k1/cellData.R)*(1/cellData.T_ref - 1/T)); % Anode chemical rate constant

% Calculate electrolyte concentration at anode and cathode current
% collector interface
ce = spm_matrices.M*ue;
cq0 = ce(1,:); % Electrolyte concentration at cathode current collector interface
cm0 = ce(4,:);  % Electrolyte concentration at anode current collector interface
ckq = ce(2,:);
cmk = ce(3,:);

ce3 = [cq0;ue(1:Q_c-1,:);ckq];
% ce2 = [ckq;ue(Q_c:Q_c+K_s-2,:);cmk];
ce1 = [cmk;ue(Q_c+K_s-1:Q_c+K_s+M_a-3,:);cm0];

ce3_avg = (spm_matrices.wq_c*ce3)/2;
ce1_avg = (spm_matrices.wm_a*ce1)/2;

% i0 = k * F * Ce^alpha_a * Cs_surf^alpha_c * (Cs_max - Cs_surf)^alpha_a
% k - temperature dependent chemical rate constant
% Ce - Electrolyte concentration (assumed to be constant as it is not
% modelled in SPM)
% alpha_a - transfer coefficient of anode (assumed to be 0.5 for symmetric
% reaction kinetics)
% alpha_c = 1 - alpha_a transfer coefficient of cathode
i03_e = k3.*(ce3.^0.5).*(css3.^0.5).*(cellData.cs3_max - css3).^0.5;
i01_e = k1.*(ce1.^0.5).*(css1.^0.5).*(cellData.cs1_max - css1).^0.5;

i03 = (spm_matrices.wq_c*i03_e)/2;
i01 = (spm_matrices.wm_a*i01_e)/2;

% Overpotential at each electrode
% eta = (2*R*T/F)*asinh(j*F/(2*i0))
eta3 = (2*cellData.R*T/cellData.F)*asinh(j3*cellData.F/(2*i03));
eta1 = (2*cellData.R*T/cellData.F)*asinh(j1*cellData.F/(2*i01));

% Overpotential due to electrolyte concentration change
eta_e = (2*cellData.R*T/(cellData.F*cellData.ce_init))*(1-cellData.t_plus_0)*(ce3_avg - ce1_avg);

% Electrolyte ohmic loss
Ke_ce0 = cellData.Ke(cellData.ce_init);
delphi_e = (i_app/Ke_ce0)*(cellData.thick1/(3*cellData.eps1e^cellData.brug) + cellData.thick2/(cellData.eps2e^cellData.brug) + cellData.thick3/(3*cellData.eps3e^cellData.brug));

% Solid ohmic loss
delphi_s = (i_app/3)*(cellData.thick1/cellData.sigma1 + cellData.thick3/cellData.sigma3);

x3 = css3/cellData.cs3_max; % Li fraction in cathode
x1 = css1/cellData.cs1_max; % Li fraction in anode
[OCP1,~,OCP3,~] = get_openCircuitPotentialLGM50(x1,x3,T,cellData);

V = real((OCP3 + eta3) - (OCP1 + eta1) + eta_e - delphi_s - delphi_e);
end