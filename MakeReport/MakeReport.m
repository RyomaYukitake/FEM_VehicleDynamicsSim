% How to Make MATLAB Report
% https://jp.mathworks.com/help/rptgen/ug/create-a-presentation-programmatically.html

Param="VehicleParamsFEM20.m";
run(Param);
import mlreportgen.ppt.*;

Model="AccSim_base.slx";

warning('off','all')


%% Default & Test Param
VeloDelay0=80;
VeloDelay=VeloDelay0;
VeloDelayTest=[0, 20, 40, 60, 100, 120, 160, 240, 320];
mu_deff0=1;
mu_deffTest=[0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95];
mu_deff=mu_deff0;
SlipEnergy0=3;
SlipEnergytest=[1, 2, 4, 5, 6, 9, 12];
SlipEnergy=SlipEnergy0;
AllSlide=1+4+length(VeloDelayTest)+length(mu_deffTest)+length(SlipEnergytest);
%% タイトルスライド
Title="Robustness Check of Traction Controll"; %タイトル
SubTille="Tsuyoshi SOGA" + newline + "2023/07/28"; %サブタイトル

ppt = Presentation(Title + ".pptx","MATLAB_Report.potx");
open(ppt);
slide = add(ppt,"Title Slide");
replace(slide,"Title",Title);
replace(slide,"Subtitle",SubTille);
slideNo=1;
slides=repmat(slide,AllSlide);

%% 制御モデルスクショ
VCMmain="VCM_main_ver3_0_0.slx";
TC="TractionControl.slx";
open(VCMmain)
open(TC)

Name=string.empty(4,0);
ModelPath=string.empty(4,0);

Name(1)="AllTractionControl";
ModelPath(1)="-sVCM_main_ver3_0_0/VDC/Traction Control";
Name(2)="TargetSpeed_inAcc";
ModelPath(2)="-sVCM_main_ver3_0_0/VDC/Traction Control/TargetSpeed_inAcc";
Name(3)="FeedForward";
ModelPath(3)="-sVCM_main_ver3_0_0/VDC/Traction Control/FF";
Name(4)="TractionControl";
ModelPath(4)="-sTractionControl";

for i=1:length(Name)
    PathPic="Pictures\Model"+Name(i)+".png";
    print(ModelPath(i), "-dpng", PathPic)
    Snap(i)=Picture(PathPic);
    
    slide = add(ppt,"Title and 1Content");
    replace(slide,"Title",Name(i));
    replace(slide,"Content1", Snap(i));
    slideNo=slideNo+1;
    slide(slideNo)=slide;
end

close_system(VCMmain)
close_system(TC)


%% VeloDelayTest
sim(Model)
PathPlot="Pictures\VeloDelay"+VeloDelay+".png";
PlotScope(PathPlot,"VeloDelay="+VeloDelay,Torque,SlipRate,LongG);
Plot0=Picture(PathPlot);

Plots1=repmat(Plot0, length(VeloDelayTest), 2);

for i=1:length(VeloDelayTest)
    VeloDelay=VeloDelayTest(i);
    sim(Model)
    PathPlot="Pictures\VeloDelay"+VeloDelay+".png";
    PlotScope(PathPlot,"VeloDelay="+VeloDelay,Torque,SlipRate,LongG);
    AddPlot=Picture(PathPlot);
    
    Plots1(i,2)=AddPlot;
    
    slide = add(ppt,"Title and 2Content");
    replace(slide,"Title","VeloDelay="+VeloDelay);
    replace(slide,"Content1",Plots1(i,1));
    replace(slide,"Content2",Plots1(i,2));
    slideNo=slideNo+1;
    slide(slideNo)=slide;
end
VeloDelay=VeloDelay0;%#ok<NASGU> 


%% mu diff test
sim(Model)
PathPlot="Pictures\mu"+mu_deff*mu_tire_F+".png";
PlotScope(PathPlot,"mu="+mu_deff* mu_tire_F,Torque,SlipRate,LongG);
Plot0=Picture(PathPlot);

Plots2=repmat(Plot0, length(mu_deffTest), 2);

for i=1:length(mu_deffTest)
    mu_deff=mu_deffTest(i);
    sim(Model)
    PathPlot="Pictures\mu"+mu_deff*mu_tire_F+".png";
    PlotScope(PathPlot,"mu="+mu_deff* mu_tire_F,Torque,SlipRate,LongG);
    AddPlot=Picture(PathPlot);
    
    Plots2(i,2)=AddPlot;
    
    slide = add(ppt,"Title and 2Content");
    replace(slide,"Title","mu difference x"+mu_deff);
    replace(slide,"Content1",Plots2(i,1));
    replace(slide,"Content2",Plots2(i,2));
    slideNo=slideNo+1;
    slide(slideNo)=slide;
end
mu_deff=mu_deff0; %#ok<NASGU> 

%% Slip Energy Test
sim(Model)
PathPlot="Pictures\SlipEnergy"+SlipEnergy+".png";
PlotScope(PathPlot,"SlipEnergy="+SlipEnergy,Torque,SlipRate,LongG);
Plot0=Picture(PathPlot);

Plots3=repmat(Plot0, length(SlipEnergytest), 2);

for i=1:length(SlipEnergytest)
    SlipEnergy=SlipEnergytest(i);
    sim(Model)
    PathPlot="Pictures\SlipEnergy"+SlipEnergy+".png";
    PlotScope(PathPlot,"SlipEnergy="+SlipEnergy,Torque,SlipRate,LongG);
    AddPlot=Picture(PathPlot);
    
    Plots3(i,2)=AddPlot;
    
    slide = add(ppt,"Title and 2Content");
    replace(slide,"Title","SlipEnergy difference");
    replace(slide,"Content1",Plots3(i,1));
    replace(slide,"Content2",Plots3(i,2));
    slideNo=slideNo+1;
    slide(slideNo)=slide;
end
SlipEnergy=SlipEnergy0;%#ok<NASGU> 


%% 終了
VeloDelay=VeloDelay0;
mu_deff=mu_deff0;
SlipEnergy=SlipEnergy0;
close(ppt);

pptview(Title + ".pptx");
pptview(Title + ".pptx",'converttopdf');

warning


%% plot作成関数
function PlotScope(Path,GraphTitle,Torque,SlipRate,LongG)
    f=figure;
    f.Position(3:4) = [560 700]; %[width heigt]
    tiledlayout(3,1)
    nexttile
    title(GraphTitle)
    PlotScope2("Torque", Torque)
    nexttile
    PlotScope2("SlipRate", SlipRate)
    nexttile
    PlotScope1("LongG",LongG)
    
    saveas(f,Path);
    close(f)
end
function PlotScope2(Name, Data)
    hold on
    grid on
    plot(Data{1}.Values)
    plot(Data{2}.Values)
    set(gca,'FontSize',16);
    ylabel(Name)
    legend(Data{1}.Name,Data{2}.Name);
    hold off
end
function PlotScope1(Name, Data)
    hold on
    grid on
    plot(Data{1}.Values)
    set(gca,'FontSize',16);
    ylabel(Name)
    hold off
end