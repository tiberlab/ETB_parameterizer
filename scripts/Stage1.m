% STAGE 1: ESTIMATING JANCU-LIKE PARAMETERS

%=================
% Input for ZB GaN
%=================
% From Table 3 of [Vogl1983Semiempirical]
wsc = -11.5490; %Ga
wpc = -5.6712; %Ga
wec = -2.60; %Ga
wsa = -25.7130; %N
wpa = -15.4388; %N
wea = -4.99; %N

% Target energies at critical k-points (see Figure 1 of [Vogl1983Semiempirical])
G1v = -16.214451;
G1c = 3.310952;
G15v = 0.000;
G15c = 11.984989;
X1v = -13.675603;
X3v = -5.533005;
X5v = -2.368161;

% Other info
a = 4.5030; %(*in Angstrom*)

%=======================================
% Step 1: follow [Vogl1983Semiempirical]
%=======================================
betas = 0.8;
betap = 0.6;
betae = 0.8;

betaws = betas*(wsc - wsa);
betawp = betap*(wpc - wpa);
betawe = betae*(wec - wea);

Esc = 0.5*(G1c + G1v + betaws);
Esa = 0.5*(G1c + G1v - betaws);
Epc = 0.5*(G15c + G15v + betawp);
Epa = 0.5*(G15c + G15v - betawp);

Vsss = -0.125*sqrt((G1c - G1v)^2 - betaws^2); %(*Vsss < 0*)
Vscpas = sqrt(3)/8*sqrt((Esc + Epa - 2*X3v)^2 - (Esc - Epa)^2); %(*Vsps > 0*)
Vsapcs = sqrt(3)/8*sqrt((Esa + Epc - 2*X1v)^2 - (Esa - Epc)^2); %(*Vpss > 0*)
Vpps = 0.125*sqrt((G15c - G15v)^2 - betawp^2) + 0.250*sqrt((G15c + G15v - 2*X5v)^2 - betawp^2); %(*Vpps > 0*)
Vppp = 0.125*sqrt((G15c - G15v)^2 - betawp^2) - 0.125*sqrt((G15c + G15v - 2*X5v)^2 - betawp^2); %(*Vppp < 0*)

%====================================
% Step 2: follow [Jancu1998Empirical]
%====================================
hbar = 1.054571817*10^(-34);
m0 = 9.109*10^(-31);
Eunit = hbar^2*(2*pi/(a*10^(-10)))^2/(2*m0)*6.2415*10^18; %(*in eV*)

Edc = (Esc + Esa)/2 + (4.9823-1.9378)*Eunit;
Eda = (Esc + Esa)/2 + (4.9823-1.9378)*Eunit;
Eec = (Esc + Esa)/2 + (6.0733-1.9378)*Eunit + betawe/2;
Eea = (Esc + Esa)/2 + (6.0733-1.9378)*Eunit - betawe/2;
% Other couplings are initialized to be zero
Eunit = 0.0;
Vees = -0.9317*Eunit;
Vsceas = -0.3093*Eunit;
Vsaecs = -0.3093*Eunit;
Vecpas = 0.2199*Eunit; %(*0 in free-electron model*)
Veapcs = 0.2199*Eunit; %(*0 in free-electron model*)
Vscdas = -0.3837*Eunit;
Vsadcs = -0.3837*Eunit;
Vecdas = -0.0198*Eunit; %(*0 in free-electron model*)
Veadcs = -0.0198*Eunit; %(*0 in free-electron model*)
Vpcdas = -0.1849*Eunit;
Vpadcs = -0.1849*Eunit;
Vpcdap = 0.6230*Eunit;
Vpadcp = 0.6230*Eunit;
Vdds = -0.2311*Eunit;
Vddp = 0.7309*Eunit;
Vddd = -0.7248*Eunit;

%====================================================
% Step 3: shift all onsite parameters so that Esc = 0
%====================================================
shift = -Esc;
Esc = Esc + shift;
Epc = Epc + shift;
Eec = Eec + shift;
Edc = Edc + shift;
Esa = Esa + shift;
Epa = Epa + shift;
Eea = Eea + shift;
Eda = Eda + shift;

%======================
% Step 4: write to file
%======================
%estimation = [Epc, Eec, Edc, Esa, Epa, Eea, Eda, ...
%              Vsss, Vsapcs, Vsaecs, Vsadcs, Vscpas, ...
%              Vpps, Vecpas, Vpadcs, Vsceas, Veapcs, ...
%              Vees, Veadcs, Vscdas, Vpcdas, Vecdas, ...
%              Vdds, Vppp, Vpadcp, Vpcdap, Vddp, Vddd];
              
%estimation = [Epc, Eec, Edc, Esa, Epa, Eea, ...
%              Vsss, Vsapcs, Vsaecs, Vsadcs, Vscpas, ...
%              Vpps, Vecpas, Vpadcs, Vsceas, Veapcs, ...
%              Vees, Veadcs, ...
%              Vppp, Vpadcp];

estimation = [Epc, Eec, Esa, Epa, Eea, ...
              Vsss, Vsapcs, Vsaecs, Vscpas, ...
              Vpps, Vecpas, Vsceas, Veapcs, ...
              Vees, Vppp];

%estimation = [Eec, Eea, ...
%              Vsaecs, Vecpas, Vsceas, Veapcs, Vees];

%timestamp = datestr(now, 'HH:MM-ddmmmyyyy');
writematrix(estimation, '../temp/estimated-point_stage1_spe.csv');
