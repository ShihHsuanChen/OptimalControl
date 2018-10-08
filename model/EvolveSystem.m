% @model/EvolveSystem.m

function  finalState = EvolveSystem(obj,initalState,ti,tf)
% initalState : inital state [4x4 matrix]
% finalState  : final state  [4x4 matrix]
% ti          : start timing [arbitrary unit: multiple of unit time]
% tf          : end timing [arbitrary unit: multiple of unit time]
    
    Ui = reshape(initalState,[],1);
    Schrodinger = @obj.Schrodinger_H_p_rf;
    opts = odeset('RelTol',1e-6,'AbsTol',1e-10);
    [Tsol,Usol] = ode45(Schrodinger,[ti tf],Ui,opts);
    Uf = Usol(end,:);
    Uf = reshape(Uf,[],round(sqrt(length(Uf))));
%    finalState = U_0'*Uf;
    finalState = Uf;
end