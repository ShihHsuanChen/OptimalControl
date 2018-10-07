% Schrodinger_H_p_rf.m
% time-dependent Schrodinger equation
function dUdt = Schrodinger_H_p_rf(t,U)
% t    : timing [arbitrary unit: multiple of unit time]
% U    : input state [16x1 matrix]
% dUdt : output dU/dt [16x1 matrix]

    global B0 gamma_e A1 A2 J T0;
    global U_0;
    i = sqrt(-1);
    U_m = reshape(U,[],4);

    B_x1 = MagneticField_X1(t);
    B_y1 = MagneticField_Y1(t);
    B_x2 = MagneticField_X2(t);
    B_y2 = MagneticField_Y2(t);

    delta_A = (A2 - A1)/2;
    ave_A = (A1 + A2)/2;
    E_uu = gamma_e*B0 + J/4 +delta_A/2;
    E_ud = -J/4 - ave_A/2;
    E_du = -J/4 + ave_A/2;
    E_dd = -gamma_e*B0 + J/4 -delta_A/2;
    H = [ E_uu 0    0    0    ;
          0    E_ud J/2  0    ;
          0    J/2  E_du 0    ;
          0    0    0    E_dd ];

    Sigma_x = [0  1; 1  0];
    Sigma_y = [0 -i; i  0];
    Sigma_z = [1  0; 0 -1];
    
    P_1 = B_x1*kron(Sigma_x,eye(2) ) + B_y1*kron(Sigma_y,eye(2) );
    P_2 = B_x2*kron(eye(2) ,Sigma_x) + B_y2*kron(eye(2) ,Sigma_y);
%    H_p = H  + (P_1 + P_2)/2;
    H_p = H  + (P_1 + P_2)/2000;

    E_w = i*(gamma_e*B0+delta_A/2)*2*pi*t/T0;
    U_0 = expm(E_w*(kron(Sigma_z,eye(2))+kron(eye(2),Sigma_z)));
%    dU_0dt = E_w*(kron(Sigma_z,eye(2))+kron(eye(2),Sigma_z))*U_0;
    dU_0dt = i*(gamma_e*B0+delta_A/2)*2*pi/T0*(kron(Sigma_z,eye(2))+kron(eye(2),Sigma_z))*U_0;

%    H_p_rf = U_0'*H_p*U_0 - i*U_0'*dU_0dt;
    H_p_rf = H_p;
    dUdt_m = -i * H_p_rf * U_m ;
    dUdt   = reshape(dUdt_m,[],1);

end


function B = MagneticField_X1(t)
    global T0;
    global ParameterArray;
    a111 = ParameterArray(1);
    a112 = ParameterArray(2);
    B = a111*cos(2*pi*t/T0) + a112*cos(2*pi*2*t/T0);
end

function B = MagneticField_Y1(t)
    global T0;
    global ParameterArray;
    a121 = ParameterArray(3);
    a122 = ParameterArray(4);
    B = a121*sin(2*pi*t/T0) + a122*cos(2*pi*2*t/T0) ;
end

function B = MagneticField_X2(t)
    global T0;
    global ParameterArray;
    a211 = ParameterArray(5);
    a212 = ParameterArray(6);
    B = a211*cos(2*pi*t/T0) + a212*cos(2*pi*2*t/T0) ;
end

function B = MagneticField_Y2(t)
    global T0;
    global ParameterArray;
    a221 = ParameterArray(7);
    a222 = ParameterArray(8);
    B = a221*sin(2*pi*t/T0) + a222*cos(2*pi*2*t/T0) ;
end
