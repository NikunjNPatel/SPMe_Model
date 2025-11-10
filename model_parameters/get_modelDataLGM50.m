function data = get_modelDataLGM50
    
%% Cell Parameters
data.C_nom = 5; % [A.h]

%% CONSTANTS
data.R = 8.314;                     % Gas constant [J.K-1.mol-1]
data.F = 96487;                     % Faraday constant [C.mol-1]
data.T_ref = 25 + 273.15;           % Reference temperature [K]

%% CURRENT COLLECTOR
data.Rc     = 0;                % Current collector resistance [ohm.m2]

%% BATTERY GEOMETRY
data.thick1 = 8.52e-5;             % Thickness of anode [m]
data.thick2 = 1.2e-5;             % Thickness of separator [m]
data.thick3 = 7.56e-5;             % Thickness of cathode [m]
data.L	= data.thick1 + data.thick2 + data.thick3;  % Thickness of the cell [m]

% Electrode active surface area [m2] (assumed equal for anode and cathode)
data.As = 0.1027; 

%% ACTIVE POROUS MATERIAL
% Solid particles of the electrodes
data.Rs1 = 5.86e-6;                 % Anode solid particles' radius [m]
data.Rs3 = 5.22e-6;                 % Cathode solid particles' radius [m]

% Volume fraction of active material [-]
data.eps1s = 0.75;                % anode
data.eps3s = 0.665;                % cathode

% Specific interfacial surface area in porous electrodes [m2/m3]
    % due to the porosity, the 'effective' surface is larger than the
    % 'geometric' surface of the electrode.
    % A_effective = as * V = specific surface * volume
    % as = 3*volume fraction / radius (Ning, Popov, 'Cycle life modeling of
    % Lithium-Ion batteries', 2004)
data.as1 = (3*data.eps1s)/data.Rs1; % anode   
data.as3 = (3*data.eps3s)/data.Rs3; % cathode

% Max Solid phase concentration [mol.m-3]
data.cs1_max = 33133;           
data.cs3_max = 63104;

% Stoichiometry limits [-]
data.x1_soc0 = 0.027;              % at 0% soc in anode
data.x1_soc1 = 0.9014;              % at 100% soc in anode
data.y3_soc0 = 0.8536;              % at 0% soc in cathode
data.y3_soc1 = 0.27;          	% at 100% soc in cathode

% Diffusion coefficient of Li in active material [m2.s-1]
% data.Ds1_ref = @(x1) 3.0321*10^(11.17*x1 - 15.11 - 1.553*exp(-((x1-0.2031)^2)/0.0006091) - 6.136*exp(-((x1-0.5375)^2)/0.06438) - 9.725*exp(-((x1-0.9144)^2)/0.0578) + 1.85*exp(-((x1-0.5953)^2)/0.001356));
data.Ea_Ds1 = 2092*data.R;

% data.Ds3_ref = @(x3) 2.7*10^(-13.96 - 0.9231*exp(-((x3-0.3216)^2)/0.002534) - 0.4066*exp(-((x3-0.4532)^2)/0.003926) - 0.993*exp(-((x3-0.8098)^2)/0.09924));
data.Ea_Ds3 = 1449*data.R;

data.Ds1_ref = 3.3e-14;     % Anode diffusion coeff [m2.s-1]
data.Ds3_ref = 4e-15;     % Cathode diffusion coeff [m2.s-1]

% Charge transfer coefficient in active material
data.alpha_1 = 0.5;
data.alpha_3 = 0.5;

% Reaction rate constant [m2.5 mol-0.5 s-1]
data.k1_ref = 6.48e-7;
data.k3_ref = 3.42e-6;

data.Ea_k1 = 35000;      % Activation energy [J/mol]
data.Ea_k3 = 17800;      % Activation energy [J/mol]

%{
    The reaction rate constant at a given temperature T is computed using
    the following Arrhenius equation:
                k(T) = k_ref*exp(Ea_k/R*(1/T_ref - 1/T))
%}

data.ce_avg = 1.0e3;    % Average electrolyte concentration

%% THERMAL MODEL PARAMETERS WITH CONVECTION BOUNDARY CONDITIONS
data.height =   70e-3;            	% 21700 height [m]
data.diam   =   21e-3;           	% 21700 diameter [m]

% Surface area to volume ratio for an 21700 cell [m-1]
data.SA_V   =   4*(1 + data.diam/data.height/2)/data.diam;
% Volume of the cell
data.Vc     =   pi*(data.diam/2)^2*data.height;

% Convective boundary condition
data.h      =   30;         % Convection heat transfer coefficient [W/m2/K]                  
data.T_amb  = 25 + 273.15;  % Ambient temperature [K]

% Thermal properties of cell in detail
L_n = 85.2e-6; % Length of negative electrode [m]
L_s = 12e-6;    % Length of separator [m]
L_p = 75.6e-6;  % Length of positive electrode [m]
L_ncc = 12e-6;  % Length of negative current collector [m]
L_pcc = 16e-6;  % Length of positive current collector [m]
kT_p = 2.331;    % Positive electrode thermal conducitivity [W/m.K]
kT_n = 1.7;    % Negative electrode thermal conductivity [W/m.K]
kT_ncc = 401;   % Negative current collector thermal conductivity [W/m.K]
kT_pcc = 237;   % Positive current collector thermal conductivity [W/m.K]
kT_s = 0.3344;    % Separator thermal conductivity [W/m.K]
rho_n = 2060;   % Negative electrode density [kg/m3]
rho_p = 3699;   % Positive electrode density [kg/m3]
rho_s = 1548;   % Separator density [kg/m3]
rho_pcc = 2702; % Positive current collector density [kg/m3]
rho_ncc = 8933; % Negative current collector density [kg/m3]
Cp_p = 700;   % Positive electrode heat capacity [J/kg.K]
Cp_n = 700;  % Negative electrode heat capacity [J/kg.K]
Cp_pcc = 897;   % Positive current collector heat capacity [J/kg.K]
Cp_ncc = 385;   % Negative current collector heat capacity [J/kg.K]
Cp_s = 700;  % Separator heat capacity [J/kg.K]

totalLength = L_n + L_s + L_p + L_pcc + L_ncc;
data.Cp     = (L_n*Cp_n + L_s*Cp_s + L_p*Cp_p + L_pcc*Cp_pcc + L_ncc*Cp_ncc)/totalLength;                  % Heat capacity [J/kg/K]
data.rho    = (L_n*rho_n + L_s*rho_s + L_p*rho_p + L_pcc*rho_pcc + L_ncc*rho_ncc)/totalLength;                 % Density [kg/m3]
data.kT = (L_n*kT_n + L_s*kT_s + L_p*kT_p + L_pcc*kT_pcc + L_ncc*kT_ncc)/totalLength; % Thermal conductivity [W/m.K]

%% Electrolyte Diffusion Parameters

data.eps1e = 0.25;  % Anode porosity
data.eps3e = 0.335; % Cathode porosity
data.eps2e = 0.47;  % Separator porosity
data.brug = 1.5;    % Bruggeman exponent
data.t_plus_0 = 0.2594; % Cation transference number with respect to solvent
data.ce_init = 1000; % initial concentration [mol.m-3]
data.sigma3 = 0.18; % Cathode conductivity [S.m-1]
data.sigma1 = 215;  % Anode conductivity [S.m-1]

data.De = @(ce) 8.794e-17.*ce.^2 - 3.972e-13.*ce + 4.862e-10;   % Electrolyte diffusivity
data.De_1 = @(ce) 2*8.794e-17.*ce - 3.972e-13;   % Electrolyte diffusivity differentiation wrt to electrolyte concentration
data.Ke = @(ce) 1.297e-10*ce.^3 - 7.94e-5*ce.^1.5 + 3.329e-3*ce;   % Ionic conductivity of electrolyte
end