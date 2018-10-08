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
    E_uu = gamma_e*B0/2 + J/4 + delta_A/2;
    E_ud = -J/4 - ave_A/2;
    E_du = -J/4 + ave_A/2;
    E_dd = -gamma_e*B0/2 + J/4 - delta_A/2;
    H = [ E_uu 0    0    0    ;
          0    E_ud J/2  0    ;
          0    J/2  E_du 0    ;
          0    0    0    E_dd ];

    Sigma_x = [0  1; 1  0];
    Sigma_y = [0 -i; i  0];
    Sigma_z = [1  0; 0 -1];
    
    P_1 = B_x1*kron(Sigma_x,eye(2) )/2 + B_y1*kron(Sigma_y,eye(2) )/2;
    P_2 = B_x2*kron(eye(2) ,Sigma_x)/2 + B_y2*kron(eye(2) ,Sigma_y)/2;
    H_p = H  + (P_1 + P_2);

    E_w = i*(gamma_e*B0+delta_A/2)*2*pi*t/T0;
    U_0 = expm(E_w*(kron(Sigma_z,eye(2))+kron(eye(2),Sigma_z)));
    dU_0dt = i*(gamma_e*B0+delta_A/2)*2*pi/T0*(kron(Sigma_z,eye(2))+kron(eye(2),Sigma_z))*U_0;

%    H_p_rf = U_0'*H_p*U_0 - i*U_0'*dU_0dt;
    H_p_rf = H_p;
    dUdt_m = -i * H_p_rf * U_m ;
    dUdt   = reshape(dUdt_m,[],1);

end


function B = MagneticField_X1(t)
    global B0 W0;
    global ParameterArray;
    a = ParameterArray(1);
    b = ParameterArray(2);
    p = ParameterArray(9);
    B = a*B0*cos(p*W0*t) + b*B0*cos(2*p*W0*t);
end

function B = MagneticField_Y1(t)
    global B0 W0;
    global ParameterArray;
    a = ParameterArray(3);
    b = ParameterArray(4);
    p = ParameterArray(9);
    B = a*B0*sin(p*W0*t) + b*B0*cos(2*p*W0*t) ;
end

function B = MagneticField_X2(t)
    global B0 W0;
    global ParameterArray;
    a = ParameterArray(5);
    b = ParameterArray(6);
    p = ParameterArray(9);
    B = a*B0*cos(p*W0*t) + b*B0*cos(2*p*W0*t) ;
end

function B = MagneticField_Y2(t)
    global B0 W0;
    global ParameterArray;
    a = ParameterArray(7);
    b = ParameterArray(8);
    p = ParameterArray(9);
    B = a*B0*sin(p*W0*t) + b*B0*cos(2*p*W0*t) ;
end
