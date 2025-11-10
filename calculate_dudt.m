function du = calculate_dudt(t, u, I, cellData, spm_matrices)
    
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

% Calculate du/dt = D*A*u + B*j
du3 = (ones(N-1,1)*Ds3).*(spm_matrices.Ap*us3) + spm_matrices.Bp*j3;
du1 = (ones(N-1,1)*Ds1).*(spm_matrices.An*us1) + spm_matrices.Bn*j1;

%% Calculate dce/dt electrolyte
A3 = ((4*cellData.eps3e^(cellData.brug - 1)/cellData.thick3^2).*ones(Q_c-1,1));
A2 = ((4*cellData.eps2e^(cellData.brug - 1)/cellData.thick2^2).*ones(K_s-1,1));
A1 = ((4*cellData.eps1e^(cellData.brug - 1)/cellData.thick1^2).*ones(M_a-1,1));

A = diag([A3;A2;A1]'); 
B3 = ones(Q_c-1,1).* -(1 - cellData.t_plus_0)*i_app/(cellData.F*cellData.thick3*cellData.eps3e);
B2 = zeros(K_s-1,1);
B1 = ones(M_a-1,1).* (1 - cellData.t_plus_0)*i_app/(cellData.F*cellData.thick1*cellData.eps1e);
B = [B3; B2; B1];
De = diag(cellData.De(ue));
De_1 = diag(cellData.De_1(ue));

due = A*De*spm_matrices.D2_321*ue + A*De_1*(spm_matrices.D1_321*ue).^2 + B;

%% Calculate temperature using thermal model
% 1. calculate surface concentration
css3 = spm_matrices.Cp(1,:)*us3 + (spm_matrices.Dp(1,:)*j3)/Ds3;    % Cathode surface concentration
css1 = spm_matrices.Cn(1,:)*us1 + (spm_matrices.Dn(1,:)*j1)/Ds1;    % Anode surface concentration

% 2. calculate temperature dependent chemical rate constant
% k = kref*exp(E/R*(1/Tref - 1/T))
k3 = cellData.k3_ref*exp((cellData.Ea_k3/cellData.R)*(1/cellData.T_ref - 1/T)); % Cathode chemical rate constant
k1 = cellData.k1_ref*exp((cellData.Ea_k1/cellData.R)*(1/cellData.T_ref - 1/T)); % Anode chemical rate constant

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
% 3. Exchange current density i0
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

% 4. Overpotential at each electrode
% eta = (2*R*T/F)*asinh(j*F/(2*i0))
eta3 = (2*cellData.R*T/cellData.F)*asinh(j3*cellData.F/(2*i03));
eta1 = (2*cellData.R*T/cellData.F)*asinh(j1*cellData.F/(2*i01));

% 5. Calculate heat generation in model
% qtot = qrev + qrea + Qc
% qtot = I_vol*T*dU/dT + I_vol*(eta_a - eta_c) + i_app^2*Rc*As/Vc
% qtot - Volumentric heat generation [W/m3]
% I_vol - Current per unit of battery volume [i_app/L]
% dU/dT - entropic coefficient of reaction at the current
x3 = css3/cellData.cs3_max; % Li fraction in cathode
x1 = css1/cellData.cs1_max; % Li fraction in anode
[~,dOCP1_dT,~,dOCP3_dT] = get_openCircuitPotentialLGM50(x1,x3,T,cellData);
qtot = (-i_app.*T.*(dOCP3_dT - dOCP1_dT))/cellData.L + (i_app.*(eta1 - eta3))/cellData.L + i_app.^2*cellData.Rc*cellData.As/cellData.Vc;

% qconv = -h*A*dT
qconv = -cellData.h*cellData.SA_V*(T - cellData.T_amb);

% rho*cp*dT/dt = qtot - qconv
dT = (1/(cellData.rho*cellData.Cp))*(qtot + qconv);

du = [du3; du1; due; dT];

end