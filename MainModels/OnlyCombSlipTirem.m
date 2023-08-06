Load=(0.5:0.05:1.5).*(M/4*g); %[N]
SlipAngle=0:0.01:0.7; %[rad]
SlipRate=0:0.01:0.5;

Fx=zeros(length(SlipRate),1);
Fy=Fx;
MaxFx=zeros(length(Load),length(SlipAngle));
MaxFy=MaxFx;
idx=MaxFx;

hold on
for k=1:length(Load)
    for j = 1:length(SlipAngle)
        for i=1:length(SlipRate)
            PureFx = FcnPureFx(Load(k), SlipRate(i), afx_f);
            Fx(i) = FcnFx(Load(k), SlipAngle(j), agx_f, PureFx) * mu_tire_F;
            PureFy = FcnPureFy(Load(k), SlipAngle(j), a_f);
            Fy(i) = FcnFy(Load(k), SlipRate(i), agy_f, PureFy) * mu_tire_F;
        end
%         figure(n)
%         plot(SlipRate,Fx)
        [MaxFx(k,j), idx(k,j)] = max(Fx);
        MaxFy(k,j) = Fy(idx(k,j));
    end

    figure(1)%スリップ角ごとの最適スリップ率
    hold on
    plot(SlipAngle,SlipRate(idx(k,:)))
    figure(2)%摩擦円
    hold on
    plot(MaxFy(k,:),MaxFx(k,:))
end

figure(1)%スリップ角ごとの最適スリップ率
xlabel("SlipAngle")
ylabel("SlipRate FxMax")

figure(2)%摩擦円
xlabel("Fy FxMax")
ylabel("Fx FxMax")

figure(3)%スリップ角ごとの最適スリップ率
hold on
for j = 1:length(SlipAngle)
    plot(Load,SlipRate(idx(:,j)))
end
xlabel("Load")
ylabel("SlipRate FxMax")

hold off

%% Functions
function PureFx = FcnPureFx(Load, SlipRate, afx_f)
    C = afx_f(1);
    D = (Load*afx_f(2) + afx_f(3))*Load;
    E = Load*afx_f(6) + afx_f(7);
    B = afx_f(4)*sin(2*atan(Load/afx_f(5)))/(D*C);
    PureFx = D*sin(C*atan(B*SlipRate-E*(B*SlipRate-atan(B*SlipRate))));
end
function Fx = FcnFx(Load, SlipAngle, agx_f, PureFx)
    C = agx_f(1);
    D = PureFx;
    E = Load*agx_f(4) + agx_f(5);
    B = agx_f(2)*sin(2*atan(Load/agx_f(3)))/(D*C);
    Fx = D*cos(C*atan(B*SlipAngle-E*(B*SlipAngle-atan(B*SlipAngle))));
end
function PureFy = FcnPureFy(Load, SlipAngle, a_f)
    C = a_f(1);
    D = (Load*a_f(2) + a_f(3))*Load;
    E = Load*a_f(6) + a_f(7);
    B = a_f(4)*sin(2*atan(Load/a_f(5)))/(D*C);
    PureFy = D*sin(C*atan(B*SlipAngle-E*(B*SlipAngle-atan(B*SlipAngle))));
end
function Fy = FcnFy(Load, SlipRate, agy_f, PureFy)
    C = agy_f(1);
    D = PureFy;
    E = Load*agy_f(4) + agy_f(5);
    B = agy_f(2)*sin(2*atan(Load/agy_f(3)))/(D*C);
    Fy = D*cos(C*atan(B*SlipRate-E*(B*SlipRate-atan(B*SlipRate))));
end