% How to Make MATLAB Report
% https://jp.mathworks.com/help/rptgen/ug/create-a-presentation-programmatically.html

import mlreportgen.ppt.*;

Model="AccSim_base.slx";
Param="VehicleParamsFEM20.m";

%% タイトルスライド
Title="Robustness Check of Traction Controll"; %タイトル
SubTille="Tsuyoshi SOGA" + newline + "2023/07/28"; %サブタイトル

ppt = Presentation(Title + ".pptx","MATLAB_Report.potx");
open(ppt);
slide1 = add(ppt,"Title Slide");
replace(slide1,"Title",Title);
replace(slide1,"Subtitle",SubTille);

%% スライド2
VeloDelay=80;
sim(Model)
PathPlot1='Pictures\Fig1.png';
PlotScope(PathPlot1,"VeloDelay="+VeloDelay,Torque,SlipRate,LongG);
Plot1=Picture(PathPlot1);

VeloDelay=20;
sim(Model)
PathPlot2='Pictures\Fig2.png';
PlotScope(PathPlot2,"VeloDelay="+VeloDelay,Torque,SlipRate,LongG);
Plot2=Picture(PathPlot2);

slide2 = add(ppt,"Title and 2Content");
replace(slide2,"Title","VeloDelay");
replace(slide2,"Content1",Plot1);
replace(slide2,"Content2",Plot2);

%% スライド3
VeloDelay=40;
sim(Model)
PathPlot2='Pictures\Fig2.png';
PlotScope(PathPlot2,"VeloDelay="+VeloDelay,Torque,SlipRate,LongG);
Plot2=Picture(PathPlot2);

slide3 = add(ppt,"Title and 2Content");
replace(slide3,"Title","VeloDelay");
replace(slide3,"Content1",Plot1);
replace(slide3,"Content2",Plot2);

%% スライド4
VeloDelay=160;
sim(Model)
PathPlot2='Pictures\Fig2.png';
PlotScope(PathPlot2,"VeloDelay="+VeloDelay,Torque,SlipRate,LongG);
Plot2=Picture(PathPlot2);

slide4 = add(ppt,"Title and 2Content");
replace(slide4,"Title","VeloDelay");
replace(slide4,"Content1",Plot1);
replace(slide4,"Content2",Plot2);

%% スライド5
VeloDelay=320;
sim(Model)
PathPlot2='Pictures\Fig2.png';
PlotScope(PathPlot2,"VeloDelay="+VeloDelay,Torque,SlipRate,LongG);
Plot2=Picture(PathPlot2);

slide5 = add(ppt,"Title and 2Content");
replace(slide5,"Title","VeloDelay");
replace(slide5,"Content1",Plot1);
replace(slide5,"Content2",Plot2);

%% スライド6
VeloDelay=80;
sim(Model)
PathPlot1='Pictures\Fig1.png';
PlotScope(PathPlot1,"VeloDelay="+VeloDelay,Torque,SlipRate,LongG);
Plot1=Picture(PathPlot1);

VeloDelay=20;
sim(Model)
PathPlot2='Pictures\Fig2.png';
PlotScope(PathPlot2,"VeloDelay="+VeloDelay,Torque,SlipRate,LongG);
Plot2=Picture(PathPlot2);

slide2 = add(ppt,"Title and 2Content");
replace(slide2,"Title","VeloDelay");
replace(slide2,"Content1",Plot1);
replace(slide2,"Content2",Plot2);
%% 終了
close(ppt);

%% 関数
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
    plot(Data{1}.Values)
    plot(Data{2}.Values)
    set(gca,'FontSize',16);
    ylabel(Name)
    legend(Data{1}.Name,Data{2}.Name);
    hold off
end
function PlotScope1(Name, Data)
    hold on
    plot(Data{1}.Values)
    set(gca,'FontSize',16);
    ylabel(Name)
    hold off
end