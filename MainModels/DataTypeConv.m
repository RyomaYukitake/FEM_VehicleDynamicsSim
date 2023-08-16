freq=100;
start=191*freq+1;
fin=200*freq;

time=1/freq*(0:fin-start);

Actual_speed_FL=L001_Actual_speed_FL(start:fin);
Actual_speed_FR=L001_Actual_speed_FR(start:fin);
Actual_speed_RL=L001_Actual_speed_RL(start:fin);
Actual_speed_RR=L001_Actual_speed_RR(start:fin);

Actual_torque_FL=L001_Actual_torque_FL(start:fin);
Actual_torque_FR=L001_Actual_torque_FR(start:fin);
Actual_torque_RL=L001_Actual_torque_RL(start:fin);
Actual_torque_RR=L001_Actual_torque_RR(start:fin);

VELOCITY_X=L001_VELOCITY_X(start:fin);
ACCEL_X=L001_ACCEL_X(start:fin);

AS_FL = timeseries(Actual_speed_FL,time);
AS_FR = timeseries(Actual_speed_FR,time);
AS_RL = timeseries(Actual_speed_RL,time);
AS_RR = timeseries(Actual_speed_RR,time);

AT_FL = timeseries(Actual_torque_FL,time);
AT_FR = timeseries(Actual_torque_FR,time);
AT_RL = timeseries(Actual_torque_RL,time);
AT_RR = timeseries(Actual_torque_RR,time);

VELO=timeseries(VELOCITY_X,time);
ACC=timeseries(ACCEL_X,time);