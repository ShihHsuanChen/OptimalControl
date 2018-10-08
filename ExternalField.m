% @model/ExternalField.m
% a class to define external (to be optimized) field

classdef ExternalField
    properties
        B0;
        W0;
        T0;
        ParameterArray;
    end
    
    methods
        function obj = ExternalField(B0,W0,T0,ParameterArray)
            obj.B0 = B0;
            obj.W0 = W0;
            obj.T0 = T0;
            obj.ParameterArray = ParameterArray;
        end
        
        function B = MagneticField_X1(obj,t)
            a = obj.ParameterArray(1);
            b = obj.ParameterArray(2);
            p = obj.ParameterArray(9);
            B = a*obj.B0*cos(p*obj.W0*obj.T0*t) + b*obj.B0*cos(2*p*obj.W0*obj.T0*t);
        end

        function B = MagneticField_Y1(obj,t)
            a = obj.ParameterArray(3);
            b = obj.ParameterArray(4);
            p = obj.ParameterArray(9);
            B = a*obj.B0*sin(p*obj.W0*obj.T0*t) + b*obj.B0*cos(2*p*obj.W0*obj.T0*t) ;
        end

        function B = MagneticField_X2(obj,t)
            a = obj.ParameterArray(5);
            b = obj.ParameterArray(6);
            p = obj.ParameterArray(9);
            B = a*obj.B0*cos(p*obj.W0*obj.T0*t) + b*obj.B0*cos(2*p*obj.W0*obj.T0*t) ;
        end

        function B = MagneticField_Y2(obj,t)
            a = obj.ParameterArray(7);
            b = obj.ParameterArray(8);
            p = obj.ParameterArray(9);
            B = a*obj.B0*sin(p*obj.W0*obj.T0*t) + b*obj.B0*cos(2*p*obj.W0*obj.T0*t) ;
        end
    end
end