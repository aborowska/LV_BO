HV = '25';
path = 'Results/';
% load([path,'HG_all_steps_check_HV',HV,'.mat']);

Step1 = [path,'Sim_FS_LV_Ca_Cb_HV',HV,'.mat'];
Step2 = [path,'Sim_FS_LV_Ca_Cb_Klotz_HV',HV,'.mat'];
Step3 = [path,'Sim_FS_LV_af_bf_HV',HV,'.mat'];
Step4 = [path,'Sim_FS_LV_a_b_HV',HV,'.mat'];


IterNo = zeros(4,1);
IterMsg = zeros(4,1);
y_HG = zeros(4,1);

IterNo(1,1) = 100;
IterMsg(1,1) = NaN;

load(Step2,'output');
IterNo(2,1) = 2*output.funcCount;
load(Step2,'exitflag');
IterMsg(2,1) = exitflag;

load(Step3,'output');
IterNo(3,1) = output.funcCount;
load(Step3,'exitflag');
IterMsg(3,1) = exitflag;
load(Step3,'fval');
y_HG(3,1) = fval;


load(Step4,'output');
IterNo(4,1) = output.funcCount;
load(Step4,'exitflag');
IterMsg(4,1) = exitflag;
load(Step4,'fval');
y_HG(4,1) = fval;