function [fval, ResLVVol, ResLVStrain] = LV_objective2(Theta, LVEDVMRI, strainData, ...
    strain_index_1, strain_index_24,...
    options)

    pres = 8;
    postprocess_only = false;
%     if (nargin < 4)
%         strain_index_1 = 1;
%         strain_index_24 = 24;
%     end
    ii = 1;
    A = Theta(ii,1);
    B = Theta(ii,2);
    Af = Theta(ii,3); 
    Bf = Theta(ii,4); 
    As = Theta(ii,5);
    Bs = Theta(ii,6); 
    Afs = Theta(ii,7);
    Bfs = Theta(ii,8);

    [LVVolumeAba, strainAbaTotal, strainComparison, SuccessB] = ...
                    Sim_LVPassiveForwardSimulationMatPres_NoGlobal(A,   B, ...
                                                      Af,  Bf, ...
                                                      As,  Bs, ...
                                                      Afs, Bfs, ...
                                                      pres, ...
                                                      postprocess_only,...
                                                      options); 
                                                  
    ResLVVol = LVVolumeAba;
    ResLVStrain = strainComparison.strainAbaTotalOnEachSlice(strain_index_1:strain_index_24)';

    fval = [bsxfun(@minus, ResLVStrain(:,(strain_index_1:1:strain_index_24)), strainData(strain_index_1:1:strain_index_24)'),...
        (ResLVVol-LVEDVMRI)/LVEDVMRI];
    fval = sum(fval.^2,2);    

end