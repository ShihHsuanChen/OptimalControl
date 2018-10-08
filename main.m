% main function

clear all;
figure(1);
% Set fix parameters

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


B0 = 59.17; % 2T
gamma_e = 1.959;
A1 = B0*gamma_e/1000;
A2 = B0*gamma_e/1000;
J  = B0*gamma_e/10000;

% Parameter set (to be optimized)
% B1x = B0( a1x*cos(p*W0*T0*t) + b1x*cos(2*p*W0*T0*t) )
% B1y = B0( a1y*sin(p*W0*T0*t) + b1y*sin(2*p*W0*T0*t) )
% B2x = B0( a2x*cos(p*W0*T0*t) + b2x*cos(2*p*W0*T0*t) )
% B2y = B0( a2x*sin(p*W0*T0*t) + b2x*sin(2*p*W0*T0*t) )
a1x = 0.001;
a1y = 0.001;
a2x = 0.001;
a2y = 0.001;
b1x = 0.0005;
b1y = 0.0005;
b2x = 0.0005;
b2y = 0.0005;
p   = 2;

Parameters = [a1x;b1x;a1y;b1y;a2x;b2x;a2y;b2y;p];

ti = 0;
tf = 0.2;
Ui = eye(4);
U_target = [ 1 0 0 0 ;
             0 1 0 0 ;
             0 0 0 1 ; 
             0 0 1 0 ];

Model = model(gamma_e,A1,A2,J,B0,Parameters,ti,tf,Ui,U_target);
calcInfidelity = @Model.calcInfidelity;
         
%[BestParameters, BestInfidelity] = Optimizer(@calcInfidelity,Parameters);
options = optimset('Display','iter','PlotFcns',@optimplotfval,'TolFun',1E-6,'TolX',1E-6);
[BestParameters, BestInfidelity] = fminsearch(calcInfidelity,Parameters,options);

disp(BestInfidelity)
disp(BestParameters)

% re-calculate result

%Model = model(gamma_e,A1,A2,J,B0,BestParameters,ti,tf,Ui,U_target);
Schrodinger = @Model.Schrodinger_H_p_rf;
[Tsol,Usol] = ode45(Schrodinger,[ti tf],Ui);

T_arr = Tsol';
dim = length(Tsol);
infidelity_arr  = zeros(1,dim);
B1x_arr         = zeros(1,dim);
B1y_arr         = zeros(1,dim);
B2x_arr         = zeros(1,dim);
B2y_arr         = zeros(1,dim);
P1_arr          = zeros(4,dim);
P2_arr          = zeros(4,dim);

for j = 1:dim
    Uf = reshape(Usol(j,:),[],4);
    
    % infidelity
    infidelity = 1-(norm(abs(trace( Uf * U_target') )))^2/16;
    infidelity_arr(j) = infidelity;
    
    % magnetic field
    B1x = Model.ControlField.MagneticField_X1(Tsol(j));
    B1y = Model.ControlField.MagneticField_Y1(Tsol(j));
    B2x = Model.ControlField.MagneticField_X2(Tsol(j));
    B2y = Model.ControlField.MagneticField_Y2(Tsol(j));
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
        psi = Uf * Ui(:,k);
        P1(k) = psi(1)*psi(1)'+psi(2)*psi(2)';
        P2(k) = psi(1)*psi(1)'+psi(3)*psi(3)';
    end
    P1_arr(:,j) = P1;
    P2_arr(:,j) = P2;
end

% Draw plots
figure(1);

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
plot(T_arr,P1_arr(1,:),'-b',T_arr,P2_arr(1,:),'-r');

subplot(5,2,4);
plot(T_arr,P1_arr(2,:),'-b',T_arr,P2_arr(2,:),'-r');

subplot(5,2,6);
plot(T_arr,P1_arr(3,:),'-b',T_arr,P2_arr(3,:),'-r');

subplot(5,2,8);
plot(T_arr,P1_arr(4,:),'-b',T_arr,P2_arr(4,:),'-r');

% infidelity
subplot(5,2,10);
plot(T_arr,infidelity_arr,'-b');

saveas(gcf,'result.fig')
