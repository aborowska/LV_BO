function val = Obj_fun(version,ResLVStrain,ResLVVol,strainData,LVEDVMRI,...
    strain_index_1,strain_index_24)
% Compute the value of the objective function (version 1 for fO1, version 2
% for fO2) for the ABAQUS outputs for strains (ResLVStrain) and volume
% (ResLVVol) given the MRI data for strains (strainData) and volume (LVEDVMRI)

% fO1 uses unnormalised volume (and thus is volume determined) <-- STEP 1
% fO2 uses normalised volume to match the strains <-- STEP 2 and further

    if (nargin == 5)
        strain_index_1 = 1;
        strain_index_24 = 24;
    end

    if (version == 1)
        val = [bsxfun(@minus, ResLVStrain(:,(strain_index_1:1:strain_index_24)), strainData(strain_index_1:1:strain_index_24)'),...
            (ResLVVol-LVEDVMRI)];
        val = sum(val.^2,2);
    elseif (version == 2)
        val = [bsxfun(@minus, ResLVStrain(:,(strain_index_1:1:strain_index_24)), strainData(strain_index_1:1:strain_index_24)'),...
            (ResLVVol-LVEDVMRI)/LVEDVMRI];
        val = sum(val.^2,2);
%     elseif (version == 3)
%         val = [bsxfun(@minus, ResLVStrain(:,(strain_index_1:1:strain_index_24)), strainData(strain_index_1:1:strain_index_24)'),...
%             (ResLVVol-LVEDVMRI)];
%         val = bsxfun(@divide, val, [ strainData(strain_index_1:1:strain_index_24)', LVEDVMRI]);        
%         val = sum(val.^2,2);
    else
        error('Incorrect version of the objective function!')
    end

end

%% the objective function
% feval_total_fmincon = zeros(N,1);
% for ii = 1:N
%     strainComparison = SimRes(ii);
% 
%     strainComparison.strainAbaTotalOnEachSlice;
%     feval_total = strainComparison.strainAbaTotalOnEachSlice(strain_index_1:1:strain_index_24) - strainData(strain_index_1:1:strain_index_24);
%     feval_total = feval_total(strainDataFlag(strain_index_1:1:strain_index_24)>0.0);
%     % feval_total = [];
%     % for i = strain_index_1:1:strain_index_24
%     %     if strainDataFlag(i)>0.0 
%     %         feval_total = [feval_total, strainComparison.strainAbaTotalOnEachSlice(i) - strainData(i)];
%     %     end   
%     % end
%     feval_total = [feval_total; (strainComparison.LVVolumeAba-options.LVEDVMRI)/options.LVEDVMRI]; %relative volume difference
%     feval_total_fmincon(ii) = sum(feval_total.^2); 
% end
 
%     LVVolumeAba = strainComparison.LVVolumeAba;
%     ResLVVol(ii,1) = LVVolumeAba;
%     ResLVStrain(ii,:) = strainComparison.strainAbaTotalOnEachSlice(strain_index_1:strain_index_24);
  