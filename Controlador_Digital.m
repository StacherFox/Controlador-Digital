clear;
clc;

%% Especifications



s = tf('s');


% Planta Te�rica
G1t = 1/(150e-6*s+1);
G2t = 9010.83/(s^2 + 82.8763*s + 9010.83);
Gt = minreal(G1t*G2t);
[numt, dent] = tfdata(Gt,'v');

% Especifica��es da planta pr�tica
ts1 = 688e-6;
Mp = 0.22;
tp = 37e-3;

% Bloco primeira ordem
tau = ts1/5;
G1 = 1/(tau*s+1);

% Bloco segunda ordem
zeta = fzero(@(x) ((log(Mp)/pi) + (x/sqrt(1-x^2))), 0.5);
wn = pi/(tp*sqrt(1-zeta^2));
G2 = (wn^2)/(s^2 + 2*zeta*wn*s + wn^2);
ts2 = 3/(zeta*wn);

% Planta Pr�tica Cont�nua
G = G1*G2;
G = minreal(G);

% % Mostrar localiza��o dos polos dominantes
% pzmap(G,Gt);
% xlim([-50 0]);


%% Controlador

% Requisitos
Mp = Mp/2;
ts = ts2/2;
zeta = fzero(@(x) ((log(Mp)/pi) + (x/sqrt(1-x^2))), 0.5);
wn = 3/(ts*zeta);

Fs_required = (10*wn*sqrt(1-zeta^2))/(2*pi);
Fs = 185;
Ts = 1/Fs;
z = tf('z',Ts);

z1 = 0.228 + 0.08i; % ponto z que deseja ser parte do lugar das ra�zes

% Planta Pr�tica Discreta
Gz = c2d(G,Ts);

% Mostrar Requisitos
figure;
pzmap(Gz,(1/((z-z1)*(z-z1'))));
zgrid(zeta,wn*Ts);
legend('Planta','Desejado')


angle_required = -(angle(evalfr(Gz,z1)) - pi);

% O controlador anula os polos complexos da planta e cont�m um integrador,
% portanto s� precisamos encontrar o valor de um p�lo para satisfazer a
% condi��o de �ngulo

zeros_planta = pole(Gz);
angulo_polo = angle(z1 - zeros_planta(1)) + angle(z1 - zeros_planta(2)) - angle(z1 - 1) - angle_required;
polo_desejado = real(z1) - imag(z1)/tan(angulo_polo);

Cz_semK = ((z-zeros_planta(1))*(z-zeros_planta(2)))/((z-polo_desejado)*(z-1));
K = abs(1/(evalfr(minreal(Cz_semK*Gz), z1)));
Cz = K*Cz_semK;

FTMA = minreal(Cz*Gz);
FTMF = minreal(feedback(FTMA,1));

[numg, deng]    = tfdata(G, 'v');
[numgz, dengz]  = tfdata(Gz,'v');
[numc, denc]    = tfdata(Cz,'v');

