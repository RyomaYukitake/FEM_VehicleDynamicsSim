main();

function main()
timestep=0.001;
Load=62.4*9.8;
Torque=10*13.75;
v_speed0=10;
w_speed0=10.5;
afx_f=[1.736519934348130,-2.903281841026900e-04,2.651682004580759,7.000047010028531e+04,3.270011055825476e+03,-0.001252141095807,1.596836353939939];
mu_tire_F=0.55;
Iwf=0.2;
Rwf=18*(25.4/1000)/2;
M=260;

% w_speed_dot=-100:100;
% Error=w_speed_dot;
% for i=1:length(w_speed_dot)
%     Error(i)=myfun(w_speed_dot(i));
% end
% plot(w_speed_dot, Error)
% grid on

    timestep=0.001;
    x0=[0.1,0.1];
    options=optimoptions('fsolve','Algorithm','Levenberg-Marquardt');
    x = fsolve(@myfun,x0,options);
    w_speed_dot1=x(1)
    v_speed_dot1=x(2)
 

    function Error=myfun(x)
        w_speed_dot=x(1);
        v_speed_dot=x(2);
        w_speed=w_speed0 + w_speed_dot*timestep;
        v_speed=v_speed0 + v_speed_dot*timestep;
        lambda=(w_speed-v_speed)/w_speed;
        C = afx_f(1);
        D = (Load*afx_f(2) + afx_f(3))*Load;
        E = Load*afx_f(6) + afx_f(7);
        B = afx_f(4)*sin(2*atan(Load/afx_f(5)))/(D*C);
        F = D*sin(C*atan(B*lambda-E*(B*lambda-atan(B*lambda))))*mu_tire_F;
        Error1=w_speed_dot-Rwf/Iwf*(Torque - Rwf*F);
        Error2=v_speed_dot-4*F/M;
        Error=[Error1; Error2];
    end

end