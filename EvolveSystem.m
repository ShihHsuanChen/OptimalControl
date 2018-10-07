function  finalState = EvolveSystem(initalState,ti,tf)
% initalState : inital state [4x4 matrix]
% finalState  : final state  [4x4 matrix]
% ti          : start timing [arbitrary unit: multiple of unit time]
% tf          : end timing [arbitrary unit: multiple of unit time]

    Ui = reshape(initalState,[],1);
    [Tsol,Usol] = ode45(@Schrodinger_H_p_rf,[ti tf],Ui);
    Uf =  Usol(end,:);
    finalState = reshape(Uf,[],4);

    
    Ufff_sol = U_0'*U_fin;
end