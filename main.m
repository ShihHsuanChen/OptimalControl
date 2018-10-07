% main function

% Set fix parameters
global B0 gamma_e A1 A2 J T0;
B0 = 1;
gamma_e = 28024.95164;
A1 = 2200;
A2 = 2300;
J  = 22;
T0 = 400*10^(-6);

a111in = 2300;
a112in = 2300;
a121in = 2300;
a122in = 2300;
a211in = 2300;
a212in = 2300;
a221in = 2300;
a222in = 2300;

Parameters = [ a111in; a112in; a121in; a122in; a211in; a212in; a221in; a222in];

BestParameters = Optimizer(@calcInfidelity,Parameters);

disp(BestParameters)