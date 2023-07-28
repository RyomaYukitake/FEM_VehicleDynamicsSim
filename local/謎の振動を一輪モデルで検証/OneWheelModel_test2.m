main();

function main()
Load=62.4*9.8;
Torque=3*13.75;
v_speed0=0;
w_speed0=0.001;
afx_f=[1.736519934348130,-2.903281841026900e-04,2.651682004580759,7.000047010028531e+04,3.270011055825476e+03,-0.001252141095807,1.596836353939939];
mu_tire_F=0.55;
Iwf=0.2;
Rwf=18*(25.4/1000)/2;
M=260;

x0=[w_speed0;v_speed0];
[t,x] = ode45(@myfun,[0 0.001],x0);

plot(t,x(:,1),'-o',t,x(:,2),'-o')
w_speed1=x(length(t),1)
v_speed1=x(length(t),2)
title('Equation with ODE45');
xlabel('Time t');
ylabel('Solution y');
legend('w_speed1','v_speed1')

    function dxdt=myfun(t,x)
        w_speed=x(1);
        v_speed=x(2);
        lambda=(w_speed-v_speed)/w_speed;
        C = afx_f(1);
        D = (Load*afx_f(2) + afx_f(3))*Load;
        E = Load*afx_f(6) + afx_f(7);
        B = afx_f(4)*sin(2*atan(Load/afx_f(5)))/(D*C);
        F = D*sin(C*atan(B*lambda-E*(B*lambda-atan(B*lambda))))*mu_tire_F;
        w_speed_dot=Rwf/Iwf*(Torque - Rwf*F);
        v_speed_dot=4*F/M;
        dxdt=[w_speed_dot; v_speed_dot];
    end

end