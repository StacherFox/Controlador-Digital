clear;
clc;


R1 = 10e3;
R2 = 10e3;
R3 = 10e3;
R4 = 10e3;
R5 = 68e3;
R6 = 24e3;
C1 = 100e-9;
C2 = 680e-9;
C3 = 15e-9;

A = [(-1/(C1*R6))                               1/(C1*R6)               0;
    (C1 - C2)/(C1*C2*R6), (-(C1*R5 + C1*R6 - C2*R5)/(C1*C2*R5*R6))      (R4)/(R3*C2*R5);
    0,                               0,                                -1/(C3*R2)];
    
B = [    0;
         0;
         1/(C3*R1)];
            
C = [1, 0, 0];

D = 0;

sys = ss(A,B,C,D);

sys_obsv = rank(obsv(A,C));
sys_ctrb = rank(ctrb(A,B));

Mp = 0.218;
ts = 55.3e-3;
zeta = -log(Mp)/sqrt(pi^2 + log(Mp)^2);
wn = -log(0.05*sqrt(1-zeta^2))/(ts*zeta);

Mp_desejado = 0.05;
ts_desejado = 20e-3;
zeta_desejado = -log(Mp_desejado)/sqrt(pi^2 + log(Mp_desejado)^2);
wn_desejado = -log(0.05*sqrt(1-zeta_desejado^2))/(ts_desejado*zeta_desejado);
sigma = zeta_desejado*wn_desejado;


s1 = -sigma + 1i*wn_desejado*sqrt(1-zeta_desejado^2);
s2 = s1';
s3 = -2*sigma;
s4 = s3;

%% Controlador de estados

Z = zeros(rank(A),1);
Ac = [A Z; -C 0];
Bc = [B; 0];
Cc = [C 0];

P = [s1 s2 s3 s4];

K = acker(Ac,Bc,P);

%% Observador

sl = -40*sigma;
Pl = [sl sl sl];
L = acker(A',C',Pl)';

%% Discretização

T_desejado = 1/((2*pi*abs(sl)));
T = 2e-4;

