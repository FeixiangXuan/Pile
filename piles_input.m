%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  MATRICE DEI DATI DI INPUT DEL PROGRAMMA  %%
%%        Pali caricati orizzontalmente      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc
format short

% **********************************
% PILE Properties
% **********************************
Lp = 14;     % lenght of the pile [m]
Dp = 1;   % diameter of the pile [m]
tp = 0.5; % thicnkess of pile section [m]
Ep = 30e3;  % Young modulus of pile section [MPa]
Ktet = 1e5;  % rotational stiffness at pile head [MN*m]

% **************************************
% imposed displacement path at pile head
% **************************************
% u0 = final imposed horizontal displacement at pile head [m]
% w0 = final imposed vertical displacement at pile head [m]

% Nstep [-]   w0 [m]    u0 [m]
displ = [
          0    0           0;
         10    0.0045      0;
         20    0.0045      0.18;
        ];

% **********************************
% NUMERICAL PARAMETERS
% **********************************
N = 20;            % number of finite element discretizing the pile [-]

antype = 'NL';   % type of analysis: 'LE'= linear elastic; 'NL' = non linear

tol = 1e-3;       % numeric tolerance for iteration
max_it=100;       % maximum number of iterations

% **************************************************************
% reading of P-y curves from input file
% **************************************************************
filePy='Py_curves.txt';

% **************************************************************
% M-teta curves for pile section
% **************************************************************
% in case of non linear behaviour of pile section, plese specify
% these values along the pile axis
% depth [m]    EJ [MN*m2]    My [MN*m]
Data_pile=[
    0             1.473           0.608;
    Lp            1.473           0.608
    ];


% **************************************************************
% SOIL MECHANICAL PARAMETERS for vertical shaft-soil interaction
% **************************************************************
% depth = depth below ground surface
% gter = unit weight of the soil
% E = representative Young of the soil
% nu = representative Poisson coefficient of the soi
% fs = unit shaft resistance

% depth [m]    gter [kN/m3]    E [MPa]     nu [-]   fs [kPa]
Data_soil=[
    0              20           0           0.25       18;
    1.4            20           5.92         0.25       18
    2.8            20           8.37         0.25       18
    4.2            20           10.25         0.25       18
    5.6            20           11.83         0.25       18
    7.0            20           13.23        0.25       18
    8.4            20           14.49        0.25       18
    9.8            20           15.65        0.25       18
    11.2           20           16.73        0.25       18
    12.6           20           17.75        0.25       18 
    Lp             20           18.71          0.25       18;
];

zw = 0;      % depht of the water table [m]
gw = 0;   % unit weight of the water [kN/m3]

% **************************************************************
% SOIL MECHANICAL PARAMETERS for vertical base-soil interaction
% **************************************************************
% vertical bearing capacity factors
Nq = 0;
Nc = 9;
% base soil cohesion
cb = 30; %[kPa]
% equivalent elastic vertical stiffness at the base
Eb = 21;  % [MPa]
% Poisson coefficient of the soil at the base
nub = 0.25;  % [MPa]
% **************************************************************


disp(' ')
disp(' ')
disp(['START:  ',datestr(now)])
disp(' ')
piles_rigid
piles_solve_incr
disp(' ')
disp(['FINISH: ',datestr(now)])
