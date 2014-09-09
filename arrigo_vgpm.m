% Arrigo 2008 VGPM

% Input fields: SST, PAR(z), Chl, daylength

% Calculating Depth Dependent Chlorophyll

% Chl0 is the surface Chl field read in from the model

% z_base = 100 % Set base of integration, e.g. 100 m euphotic zone

% for z = 1:z_base % z_base, to bottom of integration depth
%     if z = 1:20  % where z is depth in metres and we are assuming MLD = 20
% (following Pabi 2008)
%         Chl(z) = Chl0
%     else
%         Chlz(z) = Chl0 * exp (0.033 * (-z - 20))
% end

% 
 

% N.B. you don't need the algorithm to compare model-forced and
% model-output pan-Arctic Chl(z). All you need is the above exponential
% equation.

% Calculating PP from Chlorophyll

% PP = Chl(z) * (C_Chl) * G dz dt

% where

% C_Chl = 90    & C_Chl is carbon to chlorophyll ratio (Pabi 2008)

% and

% G(z,t) = 0.59 * exp(0.0633 * T) * L    % r = 0.0633 /C; G0 = 0.59 /d, G
% at T=0°C (Pabi 2008)

% where

% L = 1 - exp(-PUR / Ek)

% where

% PUR = PAR * (achlspec(lambda) / achlspecmax) % where the expression is
% integrated is over 400-700 nm at each t and depth interval and
% achlspecmax is from the dataset of interest: achlspecmax = max(achl/[chl](data))

% Assuming PUR is constant over the mixed layer depth so take av (Arrigo 2008)
% MLD = 20 m % assuming constant MLD

% PUR0 = 0 <==== junk at the moment

% For z 1:z_base
%     if z <= MLD
%        PUR*(z) = PUR
%        PUR*_MLD = PUR0 + PUR(z)

%        PUR0 = PUR_MLD
%     else if z > MLD
%        PUR*(z) = PUR / daylength
%
% 
%
% and

% Ekmax = 80 % Arrigo 2008, n.b. a southern ocean value

%
% Ek(z) = Ekmax / (1 + 2 * exp (- b * PUR*) a = 2 and b = -0.12, fitted
% constants (Arrigo 2008)
%
% where b = exp(1.089 - 2.12 * log(Ekmax))


% achl




%