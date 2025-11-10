function [ V1,dV1dT,V3,dV3dT] = get_openCircuitPotentialLGM50( x1,x3,T,data )

% Anode OCP [Graphite]
V1_ref = 1.9793 * exp(-39.3631 * x1) ...
        + 0.2482 - 0.0909 * tanh(29.8538 * (x1 - 0.1234)) ...
        - 0.04478 * tanh(14.9159 * (x1 - 0.2769)) ...
        - 0.0205 * tanh(30.4444 * (x1 - 0.6103));

% Anode Entropy [dOCP/dT]
dV1dT = (-0.1112*x1 + 0.02914 + 0.3561*exp(-((x1-0.08309)^2)/0.004616) ...
        - tanh(63.9*(x1 - (0.4955 + 0.1122))))/1000;

V1 = V1_ref + (ones(size(x1,1),1)*T-data.T_ref).*dV1dT;

% Cathode OCP [NMC]
V3_ref = -0.809 * x3 + 4.4875 ...
        - 0.0428 * tanh(18.5138 * (x3 - 0.5542)) ...
        - 17.7326 * tanh(15.789 * (x3 - 0.3117)) ...
        + 17.5842 * tanh(15.9308 * (x3 - 0.312));

dV3dT = (0.04006*exp(-((x3-0.2828)^2)/0.0009855) - 0.06656*exp(-((x3-0.8032)^2)/0.02179))/1000;

% Open-circuit potential at temperature T
V3 = V3_ref + (ones(size(x3,1),1)*T-data.T_ref).*dV3dT;
end