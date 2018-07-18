Sim_DirConfigForwardComputation_HV02;


pres = 8.0; %%mmHg, end-diastolic pressure loading
A = 0.100000;
B = 7.545964;
Af = 3.886494;
Bf = 9.880257;
As = 0.722426;
Bs = 3.605100;
Afs = 0.100000;
Bfs = 7.888011; 

postprocess_only = false;

[LVVolumeAba, strainAbaTotal,strainComparison, SuccessB] = ...
                    Sim_LVPassiveForwardSimulationMatPres(A,   B, ...
                                                      Af,  Bf, ...
                                                      As,  Bs, ...
                                                      Afs, Bfs, ...
                                                      pres, ...
                                                      postprocess_only); 