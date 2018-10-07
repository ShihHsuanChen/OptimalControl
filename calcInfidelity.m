% calcInfidelity.m
function infidelity = calcInfidelity(pars)
% pars : parameters to be optimized

    global ParameterArray;
    ParameterArray = pars;
    
    U_target = [ 1 0 0 0 ;
                 0 1 0 0 ;
                 0 0 0 1 ; 
                 0 0 1 0 ];
    t_init = 0;
    t_fin  = 0.001;
    U_init = eye(4);
    U_fin  = EvolveSystem(U_init,t_init,t_fin);

    infidelity = 1-(norm(abs(trace( U_fin * U_target') )))^2/16;
    %fprintf('%d\t%f\t%f\t%f\n',iterNo,Infidelity,trace(Uff_sol'*Uff_sol-eye(4))/4,trace(Ufff_sol'*Ufff_sol-eye(4))/4);

end


function  finalState = EvolveSystem(initalState,ti,tf)
% initalState : inital state [4x4 matrix]
% finalState  : final state  [4x4 matrix]
% ti          : start timing [arbitrary unit: multiple of unit time]
% tf          : end timing [arbitrary unit: multiple of unit time]

    global U_0;
    Ui = reshape(initalState,[],1);
    [Tsol,Usol] = ode45(@Schrodinger_H_p_rf,[ti tf],Ui);
    Uf = Usol(end,:);
    Uf = reshape(Uf,[],4);
%    finalState = U_0'*Uf;
    finalState = Uf;

end