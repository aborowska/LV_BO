function feval_total_fmincon = LVPassive_model_optimization_fmincon_ca_cb_v2(x)

workingDir = pwd();
global options;
opt_log_filename = options.opt_log_filename ;
abaqusDir = options.abaqusDir;
strainData = options.strainData;
strainDataFlag = options.strainDataFlag;
strain_index_1 = options.strain_index_1;
strain_index_24 = options.strain_index_24;

% mpara = options.mpara;

%% all parameters will be optimized
A = options.mpara.A*x(1);
B = options.mpara.B*x(2);
Af = options.mpara.Af*x(1);
Bf = options.mpara.Bf*x(2);
As = options.mpara.As*x(1);
Bs = options.mpara.Bs*x(2);
Afs = options.mpara.Afs*x(1);
Bfs = options.mpara.Bfs*x(2);

%%this is the constrain 
if A < 0.1
    A= 0.1;
end
if Af < 0.1
    Af = 0.1;
end
if As < 0.1
    As = 0.1;
end
if Afs < 0.1
    Afs = 0.1;
end

cd(abaqusDir);
fid_log = fopen(opt_log_filename ,'a');
cd(workingDir);
timestr =  datestr(clock());
fprintf(fid_log, '\n');
fprintf(fid_log, 'Step 1 optimization for Ca Cb %s\n', timestr);



%% calling one forward simulation with provided 8 unknown parameters 
pres = options.LVEDP_norm; 
[LVVolumeAba, ~,strainComparison, SuccessB] = ...
                    Sim_LVPassiveForwardSimulationMatPres(A,   B, ...
                                                      Af,  Bf, ...
                                                      As,  Bs, ...
                                                      Afs, Bfs, ...
                                                      pres, ...
                                                      false); 

fprintf(fid_log, 'finish one simulation with 8 mmHg\n\n');                                             

%%the objective function
strainComparison.strainAbaTotalOnEachSlice;
feval_total = [];
for i = strain_index_1:1:strain_index_24
    if strainDataFlag(i)>0.0 
        feval_total = [feval_total, strainComparison.strainAbaTotalOnEachSlice(i) - strainData(i)];
    end   
end
feval_total = [feval_total, (LVVolumeAba-options.LVEDVMRI)]; %relative volume difference
feval_total_fmincon = sum(feval_total.^2); 


fprintf(fid_log, 'abaqus running success for step 1: %d\n', SuccessB);
fprintf(fid_log, 'x updated:          %f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f\n', x(1), x(2), x(1), x(2), x(1), x(2), x(1), x(2));
fprintf(fid_log, 'parameters updated: %f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f\n', A, B, Af, Bf, As, Bs, Afs, Bfs);
fprintf(fid_log, 'strain: %f (target: %f)\n', mean(strainComparison.strainAbaTotalOnEachSlice(strain_index_1:strain_index_24)),   ...
                                              mean(strainData(strain_index_1: strain_index_24) ) );
fprintf(fid_log, 'Difference (total): %f\n', feval_total_fmincon);

fprintf(fid_log, 'iteration ends\n');
fprintf(fid_log, '\n');
fclose(fid_log);

assert(SuccessB==1);