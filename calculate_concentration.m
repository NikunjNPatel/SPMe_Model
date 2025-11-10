function [cs1,cs3,cs1avg,cs3avg,ce1,ce2,ce3,ce1_avg,ce3_avg] = calculate_concentration(t, u, I, cellData, spm_matrices, N, nodes)

% Extract anode, cathode state and temperature 
N = spm_matrices.N;
M_a = spm_matrices.M_a;
K_s = spm_matrices.K_s;
Q_c = spm_matrices.Q_c;
us3 = u(1:N-1, :); % Cathode U
us1 = u(N:2*N-2,:); % Anode U
ue = u(2*N-1:(2*N-2+M_a+K_s+Q_c-3),:);
T = u((2*N-2+M_a+K_s+Q_c-3)+1,:); % Temperature
% T = 25 + 273.15;

% Calculate applied current density [A.m-2]
i_app = I(t)/cellData.As;

% Calculate molar flux on the electrode [mol.m-2.s-1]
j3 = -i_app/(cellData.as3*cellData.F*cellData.thick3);  % Cathode
j1 = i_app/(cellData.as1*cellData.F*cellData.thick1); % Anode

% Calculate temperature dependent diffusion coefficient
% D = Dref*exp(E/R*(1/Tref - 1/T))
Ds3 = cellData.Ds3_ref*exp((cellData.Ea_Ds3/cellData.R)*(1/cellData.T_ref - 1./T)); % Cathode Diffusion Coefficient
Ds1 = cellData.Ds1_ref*exp((cellData.Ea_Ds1/cellData.R)*(1/cellData.T_ref - 1./T)); % Anode Diffusion Coefficient

% Calculate concentration at inner nodes (except at center r=0)
cs1 = spm_matrices.Cn*us1 + (spm_matrices.Dn*j1)./(ones(N,1)*Ds1);
cs3 = spm_matrices.Cp*us3 + (spm_matrices.Dp*j3)./(ones(N,1)*Ds3);

xr1 = nodes.xr(1:N).*cellData.Rs1;
xr3 = nodes.xr(1:N).*cellData.Rs3;

cs1avg = ((3/cellData.Rs1^3)*cellData.Rs1*spm_matrices.wn(1:N))*(xr1.^2.*cs1)/2;
cs3avg = ((3/cellData.Rs3^3)*cellData.Rs3*spm_matrices.wn(1:N))*(xr3.^2.*cs3)/2;

% Calculating concentration in electrolyte
ce = spm_matrices.M*ue;
cq0 = ce(1,:);
ckq = ce(2,:);
cmk = ce(3,:);
cm0 = ce(4,:);

ce3 = [cq0;ue(1:Q_c-1,:);ckq];
ce2 = [ckq;ue(Q_c:Q_c+K_s-2,:);cmk];
ce1 = [cmk;ue(Q_c+K_s-1:Q_c+K_s+M_a-3,:);cm0];

ce3_avg = (spm_matrices.wq_c*ce3)/2;
ce1_avg = (spm_matrices.wm_a*ce1)/2;
end