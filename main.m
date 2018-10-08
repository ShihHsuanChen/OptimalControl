% main function

clear all;
figure(1);
% Set fix parameters
global B0 gamma_e A1 A2 J T0 W0;
% Energy  : in unit [meV]
% h-bar,c : set to be 1 
% B0      : strength of stationary magnetic field (alone Z-axis)
%           in unit [(eV)^2/e/h-bar/c^2]  (1 (eV)^2/e/h-bar/c^2 = 0.0169 T = 169 G)
% gamma_e : gyromagnetic ratio
%           gamma_e = g*q/2/m_e = g*(-e)/2/0.511MeV, where g = -2.0023
%                   = 1.959 [1/MV]
% J       : electron-electron interaction term: Hee = J * Se1。Se2
%           in unit meV/h-bar^2
% A1, A2  : electron-nucleus interaction term : Hen = A * Se1。Sn1
%           in unit meV/h-bar^2
% T0      : characteristic time: use 4pi/gamma_e/B0

B0 = 118.34; % 2T
gamma_e = 1.959;
A1 = B0*gamma_e/100;
A2 = B0*gamma_e/100;
J  = B0*gamma_e/10000;
T0 = 4*pi/gamma_e/B0;
W0 = gamma_e*B0/2;

% Parameter set (to be optimized)
% B1x = B0( a1x*cos(p*W0*t) + b1x*cos(2*p*W0*t) )
% B1y = B0( a1y*sin(p*W0*t) + b1y*sin(2*p*W0*t) )
% B2x = B0( a2x*cos(p*W0*t) + b2x*cos(2*p*W0*t) )
% B2y = B0( a2x*sin(p*W0*t) + b2x*sin(2*p*W0*t) )
a1x = 0.1;
a1y = 0.1;
a2x = 0.1;
a2y = 0.1;
b1x = 0.05;
b1y = 0.05;
b2x = 0.05;
b2y = 0.05;
p   = 0.8;

Parameters = [ a1x;a1y;a2x;a2y;b1x;b1y;b2x;b2y;p];

%[BestParameters, BestInfidelity] = Optimizer(@calcInfidelity,Parameters);
options = optimset('Display','iter','PlotFcns',@optimplotfval,'TolFun',0.01);
[BestParameters, BestInfidelity] = fminsearch(@calcInfidelity,Parameters,options);

disp(BestInfidelity)
disp(BestParameters)

% re-calculate result
ti = 0;
tf = 2*T0;
U_init = eye(4);
U_target = [ 1 0 0 0 ;
             0 1 0 0 ;
             0 0 0 1 ; 
             0 0 1 0 ];

Ui = reshape(U_init,[],1);
[Tsol,Usol] = ode45(@Schrodinger_H_p_rf,[ti tf],Ui);

T_arr = Tsol';
infidelity_arr = zeros(1,length(Tsol));
B1x_arr = zeros(1,length(Tsol));
B1y_arr = zeros(1,length(Tsol));
B2x_arr = zeros(1,length(Tsol));
B2y_arr = zeros(1,length(Tsol));
P1_arr = zeros(4,length(Tsol));
P2_arr = zeros(4,length(Tsol));

for j = 1:length(Tsol)
    Uf = reshape(Usol(j,:),[],4);
    
    % infidelity
    infidelity = 1-(norm(abs(trace( Uf * U_target') )))^2/16;
    infidelity_arr(j) = infidelity;
    
    % magnetic field
    B1x = MagneticField_X1(Tsol(j));
    B1y = MagneticField_Y1(Tsol(j));
    B2x = MagneticField_X2(Tsol(j));
    B2y = MagneticField_Y2(Tsol(j));
    B1x_arr(j) = B1x;
    B1y_arr(j) = B1y;
    B2x_arr(j) = B2x;
    B2y_arr(j) = B2y;
    
    % 4 cases : (++)>(++), (+-)>(+-), (-+)>(--), (--)>(-+)
    P1 = zeros(4,1);
    P2 = zeros(4,1);
    % P1: probability of the first electron to be +
    % P2: probability of the second electron to be +
    for k = 1:4
        psi = Uf * U_init(:,1);
        P1(k) = psi(1)*psi(1)'+psi(2)*psi(2)';
        P2(k) = psi(2)*psi(2)'+psi(3)*psi(3)';
    end
    P1_arr(:,j) = P1;
    P2_arr(:,j) = P2;
end

% Draw plots
figure(2);

% magnetic field
subplot(5,2,1);
plot(T_arr,B1x_arr,'-k');

subplot(5,2,3);
plot(T_arr,B1y_arr,'-k');

subplot(5,2,5);
plot(T_arr,B2x_arr,'-k');

subplot(5,2,7);
plot(T_arr,B2y_arr,'-k');

% Probibility
subplot(5,2,2);
plot(T_arr,P1_arr(1,j),'-b',T_arr,P2_arr(1,j),'-r');

subplot(5,2,4);
plot(T_arr,P1_arr(2,j),'-b',T_arr,P2_arr(2,j),'-r');

subplot(5,2,6);
plot(T_arr,P1_arr(3,j),'-b',T_arr,P2_arr(3,j),'-r');

subplot(5,2,8);
plot(T_arr,P1_arr(4,j),'-b',T_arr,P2_arr(4,j),'-r');

% infidelity
subplot(5,2,10);
plot(T_arr,infidelity_arr,'-b');


function B = MagneticField_X1(t)
    global B0 W0;
    global ParameterArray;
    a = ParameterArray(1);
    b = ParameterArray(2);
    p = ParameterArray(9);
    B = a*B0*cos(p*W0*t) + b*B0*cos(2*p*W0*t);
end

function B = MagneticField_Y1(t)
    global B0 W0;
    global ParameterArray;
    a = ParameterArray(3);
    b = ParameterArray(4);
    p = ParameterArray(9);
    B = a*B0*sin(p*W0*t) + b*B0*cos(2*p*W0*t) ;
end

function B = MagneticField_X2(t)
    global B0 W0;
    global ParameterArray;
    a = ParameterArray(5);
    b = ParameterArray(6);
    p = ParameterArray(9);
    B = a*B0*cos(p*W0*t) + b*B0*cos(2*p*W0*t) ;
end

function B = MagneticField_Y2(t)
    global B0 W0;
    global ParameterArray;
    a = ParameterArray(7);
    b = ParameterArray(8);
    p = ParameterArray(9);
    B = a*B0*sin(p*W0*t) + b*B0*cos(2*p*W0*t) ;
end
