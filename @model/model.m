% @model/model.m

classdef model
    properties
        ControlField;
        gamma_e;
        A1;
        A2;
        J;
        B0;
        W0;
        T0;
        ti;
        tf;
        Ui;
        Utar;
    end
    
    methods
        function obj = model(gamma_e,A1,A2,J,B0,ParameterArray,ti,Ui,Utar)
            % W0  : characteristic angular frequency: use gamma_e*B0/2
            % T0  : characteristic time: T0 = 4pi/gamma_e/B0 = 2pi/W0
            obj.gamma_e = gamma_e;
            obj.A1 = A1;
            obj.A2 = A2;
            obj.J = J;
            obj.B0 = B0;
            obj.W0 = gamma_e*B0/2;
            obj.T0 = 2*pi/obj.W0;
            obj.ControlField = ExternalField(B0,obj.W0,obj.T0,ParameterArray);
            obj.ti = ti;
            obj.tf = ParameterArray(9);
            obj.Ui = Ui;
            obj.Utar = Utar;
        end
        
        finalState = EvolveSystem(obj,initalState,ti,tf);
        infidelity = calcInfidelity(obj,pars);
        dUdt       = Schrodinger_H_p_rf(obj,t,U);
        
    end
end