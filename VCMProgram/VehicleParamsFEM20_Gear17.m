disp("")
disp("現在のワークスペースにはFEM20のパラメータが入っていますが、ギア比だけはFEM17の値です。")

clear
load('efficientmap2.mat');
% Kpower=0;%power limiter parameter
g = 9.80665; %gravity[m/s^2]

%% Vectoring Parameter
% omegan = 3;
% G = 10^(1.1/20);
% tau = 0.01; %0.046;
% 
% hpfDYC = 1;%1/(1000*2*pi);
% KYMO = 1;
% KDYC = 1000-1000;

%% initial param
tspeed = 37.5; % target speed(km/h)
ini_speed=tspeed;%initial kph
v=tspeed;
V = v/3.6;%(m/s)
ini_v=ini_speed/3.6;

%% Vehicle Parameter
M = 200+60; % Vehicle mass (kg)
I = 112; % moment of inertia z-axis (kgm^2)
l = 1.530; % Wheel base (m)
df = 1.225; % Tread Front (m)
dr = 1.225; % Tread Rear (m)
WBf = 0.48; % Weight barance front
hg = 0.28; % Center of gravity(m)
delta = 15.5; % steer angle (deg)
steergearratio=6.5098;%steerangle/tiredelta
% gas = 0.4;%Throttle

Kframe=5.94*10^6;%Frame Stiffness(N/rad)

Frtoe = 0;%Front initial toe(deg) Positive=Toe out
Reartoe = -2;%Rear initial toe(deg) Positive=Toe out

krfr=0.5;%Front Roll Stiff Dist.
rcf=0.000;%Front Roll Center Height(m)
rcr=0.000;%Rear Roll Center Height(m)

castrail_Fr = 0.00;%castrail(mm)
castrail_Rr = 0.00;%castrail(mm)
kk_f = 2800*180/pi;%Fr toe stiffnes (Nm/rad)
kk_r = 4200*180/pi;%Rr toe stiffness (Nm/rad)

%% Suspension
m1f = 21; % 21.3; % 16.34003 - 3.550 - 1.02537; %Frばね下重量[kg]
m1r = 21; % 21.0; % 17.159835; %Rrばね下重量[kg]
m2f = M*WBf-m1f; %Frばね上重量[kg]
m2r = M*(1-WBf)-m1r; %Rrばね上重量[kg]
%m2 = m2f+m2r; %ばね上重量[kg]

kf = 40.0*1000; % 26.56*1000; %Ffホイール端ばね定数[N/m]
kr = 43.2*1000; % 29.42*1000; %Rrホイール端ばね定数[N/m]
kt = 111.2*1000; %タイヤばね定数[N/m]

cf = 0.3*2*sqrt(kf*m2f); %Frダンパー減衰係数[N/(m/s)] % 0.26
cr = 0.3*2*sqrt(kr*m2r); %Rrダンパー減衰係数[N/(m/s)]
%% Battery
HV_voltage=600;

%% aero devices
DFt = 200;%Total Downforce at 40kph [N]
DFb = 0.52;%Down Force Balance Rr
Cd=85;%Drag Force @40kph [N]
Cd=Cd/(40/3.6)^2;

DFf = DFt*(1-DFb);%Downforce Front @ 40kph
DFr = DFt*DFb;%Rear

CLF = DFf/40/40;%Downforce coefficient Front
CLR = DFr/40/40;%Rear

DF_f = CLF * v * v;%Downforce Front[N]
DF_r = CLR * v * v;%Rear[N]

%% gear
gear_f=14.45;%gear ratio
gear_r=14.45;%gear ratio

%% Brake
Pedal_Ratio=3;
mPiston_Area_F=pi*(0.014/2)^2;%[m^2]
mPiston_Area_R=pi*(0.014/2)^2;%[m^2]
cPiston_Area_F=2*pi*(0.0254/2)^2;%[m^2]
cPiston_Area_R=2*pi*(0.0254/2)^2;%[m^2]
eff_Radius_F=(180/2-13.5)/1000;%[m]
eff_Radius_R=(180/2-13.5)/1000;%[m]
Brake_Factor=0.66;

%% Tire & Wheel
Fd=60*2;%Rolling resistance[N]

Iwf=0.2;%wheel inertia [kgm^2]
Rwf=18*(25.4/1000)/2;%wheel rad[m]
Mwf=Iwf/Rwf/Rwf;

Iwr=0.2;%wheel inertia [kgm^2]
Rwr=18*(25.4/1000)/2;%tire rad[m]
Mwr=Iwr/Rwr/Rwr;


%Fr tire Parameter (R10-7.5)
a_f=[1.434106444513048,-3.490394354141787e-04,2.761516722659314,6.000050177239062e+04,2.990012596710455e+03,-1.736376596503330e-04,0.526021936511837];
afx_f=[1.736519934348130,-2.903281841026900e-04,2.651682004580759,7.000047010028531e+04,3.270011055825476e+03,-0.001252141095807,1.596836353939939];
agx_f=[0.230000000000000,2.936116597489246e+04,7.137436879966481e+03,0,0];
agy_f=[0.852564208514954,3.000069064154418e+04,2.070023877794581e+03,-5.089267464508360e-04,0.519364283925805];
b_f=[3.436650779629654,1.566044221042620e-05,0.028598122547093,3.408741969645286e+03,5.995144285148312e+03,-1.736376595720277e-04,0.526021936449385];

%Fr tire Parameter (R13-7.5)
%a_f=[1.62302350628608,-0.000254628770451847,2.60502900323627,80000.4116348567,3290.00847575468,-8.67746278460815e-05,0.530278245761646];
%afx_f=[2.21949035718611,-0.000283502679827773,2.55493371139218,80000.4053645282,3270.00821939796,4.14678545311404e-05,0.968120411249738];
%agx_f=[0.300000000000000,39738.8865217479,5816.38555496722,0,0];
%agy_f=[1.05269490740135,40000.8183967565,3280.03354492566,-3.34308232997439e-05,0.574302080831508];
%b_f=[3.57566496424421,2.05287041273187e-05,0.0210606397565008,5507.15342997233,5993.01226505719,-8.67746278455842e-05,0.530278245761420];

%Rr Tire Parameter
a_r=a_f;
afx_r=afx_f;
agx_r=agx_f;
agy_r=agy_f;
b_r=b_f;

mu_tire_F = 0.587-0.587+0.55;%Tire scaling Fr
mu_tire_R = 0.587-0.587+0.55;%Rr

%TireL_f=0.337222222;%Front Tire relaxation length
TireL_f=0.382916667;%10inch-Front Tire relaxation length
TireL_r=TireL_f;%Rear

%% Parameter calculation

lf = l*(1-WBf);%(m)
lr = l*WBf;%(m)
W0_F = M*WBf/2;%Front initial load
W0_R = M*(1-WBf)/2;%Rear initial load

rch=rcf*WBf+rcr*(1-WBf);%Roll center height at center of gravity

%d =delta*pi/180;%(radians)

Rr_toe = Reartoe*pi/180;%(rad)

zc_f = castrail_Fr/1000;%(mm)
zc_r = castrail_Rr/1000;%(mm)




T_TireR = TireL_r/V;%Tire time constant Rear
T_TireF = TireL_f/V;%Front



% w=8;%Dummy parameter

%TCS
% I_tcs = 10; %繝医Λ繧ｳ繝ｳ縺ｮ遨榊?繧ｲ繧､繝ｳ
% lambda = 0.25;
% Gdis_r = l *WBf; %繝ｪ繧｢繧ｿ繧､繝､縺九ｉ驥榊ｿ?縺ｾ縺ｧ縺ｮ霍晞屬[m]
% launch_fin = 0; %Launch邨ゆｺ?縺ｮ霆贋ｽ馴??
% trq_launch = 1; %Launch譎ゅ?ｮ譛?螟ｧ繝医Ν繧ｯ蛻ｶ髯撰ｼ亥牡蜷茨ｼ?
% mu = 0.7;
% wc_hpf = 100;
% wc = 0.01; %蜉?騾溷ｺｦ縺ｫ縺九￠繧記PF繧ｫ繝?繝医が繝募捉豕｢謨ｰ
% sliprate = 0.25;
% sigma_tcs = 10;%[km/h]
% Ttrq = 0.01;%繝医Ν繧ｯ蜃ｺ蜉帷ｳｻ縺ｮLPF譎ょｮ壽焚

%AMK
% KP = 40;
% TN = 0.020;
% TD = 0;
% Mmax = 21; %[Nm]

%VCM Speed controller
% Kp = 40;
% Tn = 0.020;

%simulation parameters
%sim('three_dof_four_wheel_model_ver2_180530')

%display Max latG
%max(G)


%% Skid Time
M_sim = 220+60:10:230+60;
l_sim = 1.530:0.050:1.530; % 繝ｬ繧ｮ繝･譛?遏ｭ+5mm?ｽ槭ヵ繝ｬ繝ｼ繝?蜈ｨ髟ｷ縺御ｼｸ縺ｳ縺ｪ縺?遞句ｺｦ
                           % 1530+(210919譎らせ縺ｮ繝輔Ξ繝ｼ繝?縺ｮFBH~Fr雜ｳ轤ｹ髢楢ｷ晞屬)
d_sim = 1.250:0.050:1.250;
WBf_sim = 0.50:0.01:0.50;
hg_sim = 0.28:0.01:0.28;
Reartoe_sim = -3:0.5:0;
dyc_sim = 140:-10:70;

skidR = zeros(length(M_sim),length(l_sim),length(d_sim), ...
              length(WBf_sim), length(hg_sim));
skidtime = zeros(length(M_sim),length(l_sim),length(d_sim), ...
                 length(WBf_sim), length(hg_sim));
skidpara = zeros(length(M_sim),length(l_sim),length(d_sim), ...
                 length(WBf_sim), length(hg_sim), 5); % delta,tspeed,Reartoe,dyc,beta

skidR_temp = 0;
skidtime_temp = 100;
delta_temp = 0;
Vnow_temp = 0;
dyc_temp = 0; 
beta_temp = 0;

% rwf = [l*(1-WBf),df/2];
% rwr = [-l*WBf,dr/2];
% rtf = sumsqr([7.5*25.4/1000,Rwf])^(1/2)/2 ...
%       *[cos(atan((7.5*25.4/1000)/(Rwf*2)) + deg2rad(delta)), ...
%         sin(atan((7.5*25.4/1000)/(Rwf*2)) + deg2rad(delta))];
% rtr = sumsqr([7.5*25.4/1000,Rwr])^(1/2)/2 ...
%       *[cos(pi/2 + atan(Rwr*2/(7.5*25.4/1000)) + deg2rad(Reartoe)), ...
%         sin(pi/2 + atan(Rwr*2/(7.5*25.4/1000)) + deg2rad(Reartoe))];

% 蜷?繝代Λ繝｡繧ｿ繝ｼ縺ｮ諢溷ｺｦ隧穂ｾ｡
% for i = 1:length(M_sim)
%     M = M_sim(i);
%     for j = 1:length(l_sim)
%         l = l_sim(j);
%         for k = 1:length(d_sim)
%             df = d_sim(k);
%             dr = d_sim(k);
%             for m = 1:length(WBf_sim)
%                 WBf = WBf_sim(m);
%                 for n = 1:length(hg_sim)
%                     hg = hg_sim(n);
%                     for p = 1:length(Reartoe_sim)
%                         Reartoe = Reartoe_sim(p);
%                         for o = 1:length(dyc_sim)
%                             dyc = dyc_sim(o);
%                             % 繧ｷ繝溘Η繝ｬ繝ｼ繧ｷ繝ｧ繝ｳ
%                             display(num2str([M,l,df,dr,WBf,hg,Reartoe,dyc]))
%                             sim('FourWheelsim_VDCforSkid_210920');
%                             if Div.Data(end) == 0
%                                 % 謖吝虚縺梧険蜍輔?ｻ逋ｺ謨｣縺帙★縺ｫ繧ｷ繝溘Η繝ｬ繝ｼ繧ｷ繝ｧ繝ｳ邨ゆｺ?縺励◆繧会ｼ?
%                                 % 譛?螟ｧdyc繧ｲ繧､繝ｳ縺ｨ縺励※繝ｫ繝ｼ繝礼ｵゆｺ?
%                                 break
%                             end
%                         end
%                         if and(p ~= 1, skidtime_sim.Data(end) >= skidtime_temp)
%                             % 繝医?ｼ螟画峩縺ｫ繧医ｊ繧ｿ繧､繝?縺悟｢怜刈縺励◆繧会ｼ悟燕蝗槭そ繝?繝医ｒ譛?騾溘→縺励※險倬鹸
%                             skidR(i,j,k,m,n) = skidR_temp;
%                             skidtime(i,j,k,m,n) = skidtime_temp;
%                             skidpara(i,j,k,m,n,:) = [delta_temp; Vnow_temp; ...
%                                                      Reartoe_sim(p-1); dyc_temp; beta_temp];
%                         end
%                         skidR_temp =  skidR_sim.Data(end);
%                         skidtime_temp = skidtime_sim.Data(end);
%                         delta_temp = delta_sim.Data(end);
%                         Vnow_temp = Vnow_sim.Data(end);
%                         dyc_temp = dyc;
%                         beta_temp = beta.Data(end);
%                         if and(p ~= 1, skidtime_sim.Data(end) >= skidtime_temp)
%                             % 繝医?ｼ螟画峩縺ｫ繧医ｊ繧ｿ繧､繝?縺悟｢怜刈縺励◆繧会ｼ悟燕蝗槭そ繝?繝医ｒ譛?騾溘→縺励※繝ｫ繝ｼ繝礼ｵゆｺ?&繝ｪ繧ｻ繝?繝?
%                             delta_temp = 0;
%                             Vnow_temp = 0;
%                             dyc_temp = 0; 
%                             beta_temp = 0;
%                             break
%                         elseif p == length(Reartoe_sim)
%                             % 譛?蠕後?ｮ繧ｻ繝?繝医∪縺ｧ繧ｿ繧､繝?縺梧ｸ帛ｰ代＠縺溘ｉ?ｼ後◎縺ｮ繧ｻ繝?繝医ｒ譛?騾溘→縺吶ｋ
%                             skidR(i,j,k,m,n) = skidR_temp;
%                             skidtime(i,j,k,m,n) = skidtime_temp;
%                             skidpara(i,j,k,m,n,:) = [delta_temp; Vnow_temp; ...
%                                                      Reartoe_sim(p-1); dyc_temp; beta_temp];
%                             delta_temp = 0;
%                             Vnow_temp = 0;
%                             dyc_temp = 0; 
%                             beta_temp = 0;
%                         end
% %                         if dycDiv.Data(end) == 2
% %                             break
% %                         end
%                     end
%                 end
%             end
%         end
%     end
% end

%% p.116
% St = ScopeData.signals(1,1).values(:,1);% front tire steer angle (rad)
% B = ScopeData.signals(1,2).values(:,1);% side slip Angle (rad)
% r = ScopeData.signals(1,3).values(:,1);% yaw rate (rad/s)
% % The frequency response is calcilated by using Fourier transform.
% gr = etfe([r,St],[],2^15,0.001);
% % The gain and the phase angle of the frequeny response are calculated.
% [amp,phase,w] = bode(gr);
% % An extra dimension is deleted.
% amp = squeeze(amp);
% pahse = squeeze(phase);
% Drawing 
