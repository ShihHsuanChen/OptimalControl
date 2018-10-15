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


B0 = 59.17; % 1T
gamma_e = 1.959;
A1 = 0.116;%B0*gamma_e/1000;
A2 = 0.116;%B0*gamma_e/1000;
J  = 0.0116;%B0*gamma_e/10000;

% Parameter set (to be optimized)
% Gaussian Field
a1x   = 0.5;
a1y   = -1.;
Wext1 = 0;
phi1  = 0;
a2x   = 2.;
a2y   = -0.8;
Wext2 = 0;
phi2  = 0;
Tw    = 0.5;

Parameters = [a1x;a1y;Wext1;phi1;a2x;a2y;Wext2;phi2;Tw];
ParaLimU = [ 10; 10;Wext1;phi1; 10; 10;Wext2;phi2;Tw];
ParaLimL = [-10;-10;Wext1;phi1;-10;-10;Wext2;phi2;Tw];

ti = 0;
Ui = eye(4);
U_target = [ 1 0 0 0 ;
             0 1 0 0 ;
             0 0 0 1 ; 
             0 0 1 0 ];

Model = model(gamma_e,A1,A2,J,B0,Parameters,ti,Ui,U_target);
calcInfidelity = @Model.calcInfidelity;
         
%[BestParameters, BestInfidelity] = Optimizer(@calcInfidelity,Parameters);
options = optimset('Display','iter','PlotFcns',@optimplotfval,'TolFun',1E-6,'TolX',1E-6);
[BestParameters, BestInfidelity] = fminsearchbnd(calcInfidelity,Parameters,ParaLimL,ParaLimU,options);

disp(BestInfidelity)
disp(BestParameters)

% re-calculate result

Model = model(gamma_e,A1,A2,J,B0,BestParameters,ti,Ui,U_target);
Schrodinger = @Model.Schrodinger_H_p_rf;
[Tsol,Usol] = ode45(Schrodinger,[Model.ti Model.tf],Model.Ui);

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
    if j==dim
        save('Uf.mat','Uf','-mat');
    end
    
    % infidelity
    %infidelity = 1-(norm(abs(trace( Uf * U_target') )))^2/16;
    infidelity = 1-trace( abs( Uf * U_target') )^2/16;
    infidelity_arr(j) = infidelity;
    
    % magnetic field
    %B1 = Model.ControlField.MagneticField_atom1(Tsol(j));
    %B2 = Model.ControlField.MagneticField_atom2(Tsol(j));
    B1 = Model.ControlField.GaussianField_atom1(Tsol(j));
    B2 = Model.ControlField.GaussianField_atom2(Tsol(j));
    B1x_arr(j) = B1.x;
    B1y_arr(j) = B1.y;
    B2x_arr(j) = B2.x;
    B2y_arr(j) = B2.y;
    
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
