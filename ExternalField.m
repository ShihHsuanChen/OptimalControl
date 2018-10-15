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
        
        
        
        function B = GaussianField_atom1(obj,t)
            ax = obj.ParameterArray(1);
            ay = obj.ParameterArray(2);
            wext = obj.ParameterArray(3);
            phi = obj.ParameterArray(4);
            Tw = obj.ParameterArray(9);
            t0 = Tw/2;
            St = Tw/20;
            A = exp(-(t-t0)^2/2/St/St);
            B.x = ax * obj.B0 * A * cos(wext*obj.W0*obj.T0*t);
            B.y = ay * obj.B0 * A * cos(wext*obj.W0*obj.T0*t+phi);
            B.z = 0;
        end
    
        function B = GaussianField_atom2(obj,t)
            ax = obj.ParameterArray(5);
            ay = obj.ParameterArray(6);
            wext = obj.ParameterArray(7);
            phi = obj.ParameterArray(8);
            Tw = obj.ParameterArray(9);
            t0 = Tw/2;
            St = Tw/20;
            A = exp(-(t-t0)^2/2/St/St);
            B.x = ax * obj.B0 * A * cos(wext*obj.W0*obj.T0*t);
            B.y = ay * obj.B0 * A * cos(wext*obj.W0*obj.T0*t+phi);
            B.z = 0;
        end
        
        % Parameter set (to be optimized)
        % B1x = B0( a1x*cos(p*W0*T0*t) + b1x*cos(2*p*W0*T0*t) )
        % B1y = B0( a1y*sin(p*W0*T0*t) + b1y*cos(2*p*W0*T0*t) )
        % B2x = B0( a2x*cos(p*W0*T0*t) + b2x*cos(2*p*W0*T0*t) )
        % B2y = B0( a2x*sin(p*W0*T0*t) + b2x*cos(2*p*W0*T0*t) )
        % a1x = 0.001;
        % a1y = 0.001;
        % a2x = 0.001;
        % a2y = 0.001;
        % b1x = 0.0005;
        % b1y = 0.0005;
        % b2x = 0.0005;
        % b2y = 0.0005;
        % p   = 2;

        % Parameters = [a1x;b1x;a1y;b1y;a2x;b2x;a2y;b2y;p];
    
        function B = SecondOrderField_atom1(obj,t)
            a = obj.ParameterArray(1);
            b = obj.ParameterArray(2);
            c = obj.ParameterArray(3);
            d = obj.ParameterArray(4);
            p = obj.ParameterArray(9);
            B.x = a*obj.B0*cos(p*obj.W0*obj.T0*t) + b*obj.B0*cos(2*p*obj.W0*obj.T0*t);
            B.y = c*obj.B0*sin(p*obj.W0*obj.T0*t) + d*obj.B0*cos(2*p*obj.W0*obj.T0*t);
            B.z = 0;
        end
    
        function B = SecondOrderField_atom2(obj,t)
            a = obj.ParameterArray(5);
            b = obj.ParameterArray(6);
            c = obj.ParameterArray(7);
            d = obj.ParameterArray(8);
            p = obj.ParameterArray(9);
            B.x = a*obj.B0*cos(p*obj.W0*obj.T0*t) + b*obj.B0*cos(2*p*obj.W0*obj.T0*t);
            B.y = c*obj.B0*sin(p*obj.W0*obj.T0*t) + d*obj.B0*cos(2*p*obj.W0*obj.T0*t);
            B.z = 0;
        end
    
    end
end