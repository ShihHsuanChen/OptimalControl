% @model/calcInfidelity.m

function infidelity = calcInfidelity(obj,pars)
% pars : parameters to be optimized
    format long;
    obj.ControlField.ParameterArray = pars;
    obj.tf = pars(9);
    
    U_target = obj.Utar;
    t_init = obj.ti;
    t_fin  = obj.tf;
    U_init = obj.Ui;
    U_fin  = obj.EvolveSystem(U_init,t_init,t_fin);

    %infidelity = 1-(norm(abs(trace( U_fin * U_target') )))^2/16;
    infidelity = 1-trace(abs(U_fin*U_target'))^2/16;

end

