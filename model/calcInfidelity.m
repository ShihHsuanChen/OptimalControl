% @model/calcInfidelity.m

function infidelity = calcInfidelity(obj,pars)
% pars : parameters to be optimized
    format long;
    obj.ControlField.ParameterArray = pars;
    
    U_target = obj.Utar;
    t_init = obj.ti;
    t_fin  = obj.tf;
    U_init = obj.Ui;
    U_fin  = obj.EvolveSystem(U_init,t_init,t_fin);

    %fprintf('%.16f\n',trace(U_fin*U_fin')/4);
    infidelity = 1-(norm(abs(trace( U_fin * U_target') )))^2/16;

end

