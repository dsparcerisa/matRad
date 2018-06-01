% Estructuracion del programa

% 1 - (Deshabilitado) Calculo y optimizacion para dosis fisica
% 2 - Calculo y optimizacion de dosis considerando el RBE = 1.1
% 3 - Calculo y optimizacion de dosis para el modelo McNamara
% 4 - Calculo y optimizacion de dosis para el modelo de planificacion UCM
% 5 - Graficas de perfil de dosis
% 6 - Graficas de perfil de dosis vs RBE
% 7 - Graficas de dosis 2D
% 8 - Representacion de DVH para todas las VOI
% 9 - (Deshabilitado) Calculo de RBE medio para RBEMCN con las diversas optimizaciones
% 10 - Comparacion de DVH para VOIs especificas
% 11 - Calculo de las estadisticas generales de dosis
% 12 - Renormalizacion y recalculo de estadisticas
% 13 - Calculo de NTCP
% 14 - "Resumen" de valores de NTCP


function  [ResultPhysical, ResultConstRBE, ResultRBEMCN, ResultRBEUCM, DoseStatistics, NTCP, meanNTCP, NTCPMCNall, Renorm] = ...
    prueba_NTCP(cst, pln, ct, phantomtype, DoseStatistics, GraphSel, DoseRecalc, DoseResults, StatsRef, CompDVH)

if ~isempty(DoseResults)
    ResultConstRBE = DoseResults{1,1};
    ResultRBEMCN = DoseResults{1,2} ;
    ResultRBEUCM = DoseResults{1,3} ;
end



%% 1 - Calculo y optimizacion para dosis fisica
%
% quantityOpt         = 'physicalDose';     % options: physicalDose, constRBE, effect, RBExD
% modelName           = 'none';             % none: for photons, protons, carbon                                    MCN: McNamara-variable RBE model for protons
%                                           % WED: Wedenberg-variable RBE model for protons
%
% scenGenType = 'nomScen';            % scenario creation type 'nomScen' 'wcScen' 'impScen' 'rndScen'
%
% % retrieve model parameters
% pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
% % retrieve scenarios for dose calculation and optimziation
% pln.multScen = matRad_multScen(ct,scenGenType);
% % generate steering file
% stf = matRad_generateStf(ct,cst,pln);
%
% % dose calculation
% if strcmp(pln.radiationMode,'photons')
%     dij = matRad_calcPhotonDose(ct,stf,pln,cst);
% elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'helium') || strcmp(pln.radiationMode,'carbon')
%     dij = matRad_calcParticleDose(ct,stf,pln,cst);
% end
%
% resultGUI = matRad_fluenceOptimization(dij,cst,pln);
%
% pause(1.5);
% close all
% clc
%
% ResultPhysical.Optimized.dij = dij;
% ResultPhysical.Optimized.resultGUI = resultGUI;
%
%
% %% Recalculo de dosis para ConstRBE Y RBEMCN
%
% [~ ,ResultPhysical.ConstRBEreCalc.resultGUI] = prueba_RecalcDose(ResultPhysical.Optimized.resultGUI, ct, stf, pln, cst, 'RBExD','constRBE');
%
% [~ ,ResultPhysical.RBEMCNreCalc.resultGUI] = prueba_RecalcDose(ResultPhysical.Optimized.resultGUI, ct, stf, pln, cst, 'effect', 'MCN' );
%
% clear dij resultGUI quantityOpt modelName
ResultPhysical = 'Not evaluated';

%% 2 - Calculo y optimizacion de dosis considerando el RBE = 1.1

if DoseRecalc{1,1}{1,1} > 0
    
        quantityOpt         = 'RBExD';     % options: physicalDose, constRBE, effect, RBExD
        modelName           = 'constRBE';             % none: for photons, protons, carbon                                    MCN: McNamara-variable RBE model for protons
        % WED: Wedenberg-variable RBE model for protons
        
        
        scenGenType = 'nomScen';            % scenario creation type 'nomScen' 'wcScen' 'impScen' 'rndScen'
        
        % retrieve model parameters
        pln.bioParam = matRad_bioModel(pln.radiationMode, quantityOpt, modelName);
        % retrieve scenarios for dose calculation and optimziation
        pln.multScen = matRad_multScen(ct,scenGenType);
        % generate steering file
        stf = matRad_generateStf(ct,cst,pln);
        
    if DoseRecalc{1,1}{1,2} == 0
        clear ResultConstRBE
                
        % dose calculation
        if strcmp(pln.radiationMode,'photons')
            dij = matRad_calcPhotonDose(ct,stf,pln,cst);
        elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'helium') || strcmp(pln.radiationMode,'carbon')
            dij = matRad_calcParticleDose(ct,stf,pln,cst);
        end
        
        resultGUI = matRad_fluenceOptimization(dij,cst,pln);
        
        pause(1.5);
        close all
        clc
        
        ResultConstRBE.Optimized.dij = dij;
        ResultConstRBE.Optimized.resultGUI = resultGUI;
        
        % Recalculo de dosis para RBEMCN 
        
        [~ ,ResultConstRBE.RBEMCNreCalc.resultGUI] = prueba_RecalcDose(resultGUI, ct, stf, pln, cst, 'effect', 'MCN');
        clear dij resultGUI quantityOpt modelName
        
    else
        dij = ResultConstRBE.Optimized.dij;
        resultGUI = matRad_fluenceOptimization(dij,cst,pln);
        
        pause(1.5);
        close all
        clc
        
        ResultConstRBE.Optimized.dij = dij;
        ResultConstRBE.Optimized.resultGUI = resultGUI;
        clear dij resultGUI quantityOpt modelName
    end
end

%% 3 - Calculo y optimizacion de dosis para el modelo McNamara

if DoseRecalc{2,1}{1,1} > 0
        
    quantityOpt         = 'RBExD';      % options: physicalDose, constRBE, effect, RBExD
    modelName           = 'MCN';             % none: for photons, protons, carbon                                    MCN: McNamara-variable RBE model for protons
    % WED: Wedenberg-variable RBE model for protons
    
    scenGenType = 'nomScen';            % scenario creation type 'nomScen' 'wcScen' 'impScen' 'rndScen'
    
    % retrieve model parameters
    pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
    pln.bioOpt = 1;
    
    % retrieve scenarios for dose calculation and optimziation
    pln.multScen = matRad_multScen(ct,scenGenType);
    
    stf = matRad_generateStf(ct,cst,pln);
    
    if DoseRecalc{2,1}{1,2} == 0
    clear ResultRBEMCN

    % dose calculation
    if strcmp(pln.radiationMode,'photons')
        dij = matRad_calcPhotonDose(ct,stf,pln,cst);
    elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'helium') || strcmp(pln.radiationMode,'carbon')
        dij = matRad_calcParticleDose(ct,stf,pln,cst);
    end
    
    resultGUI = matRad_fluenceOptimization(dij,cst,pln);
    pause(1.5);
    close all
    clc
    
    ResultRBEMCN.Optimized.dij = dij;
    ResultRBEMCN.Optimized.resultGUI = resultGUI;
    
    % Recalculo de dosis para ConstRBE
    
    [~ ,ResultRBEMCN.ConstRBEreCalc.resultGUI] = prueba_RecalcDose(ResultRBEMCN.Optimized.resultGUI, ct, stf, pln, cst,'RBExD','constRBE');
    clear dij resultGUI quantityOpt modelName
    
    else
        dij = ResultRBEMCN.Optimized.dij;
        resultGUI = matRad_fluenceOptimization(dij,cst,pln);
        
        pause(1.5);
        close all
        clc
        
        ResultConstRBE.Optimized.dij = dij;
        ResultConstRBE.Optimized.resultGUI = resultGUI;
        clear dij resultGUI quantityOpt modelName
    end
end

%% 4 - Calculo y optimizacion de dosis para el modelo de planificacion UCM

if DoseRecalc{3,1}{1,1} > 0
        
    quantityOpt         = 'RBExD';
    modelName           = 'UCM';             % none: for photons, protons, carbon                                    MCN: McNamara-variable RBE model for protons
    % WED: Wedenberg-variable RBE model for protons
    
    scenGenType = 'nomScen';            % scenario creation type 'nomScen' 'wcScen' 'impScen' 'rndScen'
    
    % retrieve model parameters
    pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
    
    % retrieve scenarios for dose calculation and optimziation
    pln.multScen = matRad_multScen(ct,scenGenType);
    
    stf = matRad_generateStf(ct,cst,pln);
    
    if DoseRecalc{3,1}{1,2} == 0
    clear ResultRBEUCM

    % dose calculation
    if strcmp(pln.radiationMode,'photons')
        dij = matRad_calcPhotonDose(ct,stf,pln,cst);
    elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'helium') || strcmp(pln.radiationMode,'carbon')
        dij = matRad_calcParticleDose(ct,stf,pln,cst);
    end
    
    resultGUI = matRad_fluenceOptimization(dij,cst,pln);
    pause(1.5);
    close all
    clc
    
    ResultRBEUCM.Optimized.dij = dij;
    ResultRBEUCM.Optimized.resultGUI = resultGUI;
    
    % Recalculo de dosis para ConstRBE Y RBEMCN
    
    [~ ,ResultRBEUCM.ConstRBEreCalc.resultGUI] = prueba_RecalcDose(ResultRBEUCM.Optimized.resultGUI, ct, stf, pln, cst, 'RBExD','constRBE');
    [~ ,ResultRBEUCM.RBEMCNreCalc.resultGUI] = prueba_RecalcDose(ResultRBEUCM.Optimized.resultGUI, ct, stf, pln, cst, 'effect', 'MCN');
    clear dij resultGUI quantityOpt modelName
    
    else
        dij = ResultConstRBE.Optimized.dij;
        resultGUI = matRad_fluenceOptimization(dij,cst,pln);
        
        pause(1.5);
        close all
        clc
        
        ResultConstRBE.Optimized.dij = dij;
        ResultConstRBE.Optimized.resultGUI = resultGUI;
        clear dij resultGUI quantityOpt modelName
    end

end

% 4.5 probar a renormalizar

%% 5 - Graficas de perfil de dosis

if GraphSel(1) > 0
    % ProfileType = longitudinal // lateral
    % DisplayOption = physicalDose // RBExD // physical_vs_RBExD
    % Para hacer comparaciones entre modelos DisplayOption == RBExD // physical_vs_RBExD
    % prueba_DoseGraphs (ct, pln, cst, NumBeam, ProfileType, DisplayOption, Result, Model1, Result2, Model2)
    
    % prueba_DoseGraphs (ct, pln, cst,1, 'longitudinal', 'physical_vs_RBExD', ResultPhysical.ConstRBEreCalc.resultGUI, 'ConstRBE',ResultPhysical.RBEMCNreCalc.resultGUI,'RBEMCN')
    % title('Physical dose optimized')
    %
    prueba_DoseGraphs (ct, pln, cst,1, 'longitudinal', 'physical_vs_RBExD', ResultConstRBE.Optimized.resultGUI, 'ConstRBE', ResultConstRBE.RBEMCNreCalc.resultGUI, 'RBEMCN')
    title('ConstRBE optimized')
    
    prueba_DoseGraphs (ct, pln, cst,1, 'longitudinal', 'physical_vs_RBExD', ResultRBEMCN.Optimized.resultGUI, 'RBEMCN',ResultRBEMCN.ConstRBEreCalc.resultGUI,'ConstRBE')
    title('RBEMCN optimized')
    
    prueba_DoseGraphs (ct, pln, cst,1, 'longitudinal', 'physical_vs_RBExD', ResultRBEUCM.RBEMCNreCalc.resultGUI, 'RBEMCN',ResultRBEUCM.ConstRBEreCalc.resultGUI,'ConstRBE')
    title('RBEUCM optimized')
    %
    % %Const RBE optimized
    % fig = figure();
    % hold off
    % plot(ct.x, squeeze(ResultConstRBE.Optimized.resultGUI.physicalDose(:, 92, 37)),'k','LineWidth',3)
    % xlabel('X coordinate (mm)');
    % ylabel('Dose (Gy) // Biological Dose (RBEGy)');
    % hold on
    % plot(ct.x, 1.1*squeeze(ResultConstRBE.Optimized.resultGUI.RBExD(:, 92, 37)),'b','LineWidth',3)
    % plot(ct.x, squeeze(ResultConstRBE.RBEMCNreCalc.resultGUI.RBExD(:, 92, 37)),'r','LineWidth',3)
    % axis([-150 50 0 3]);
    % title('Uniform RBE optimized')
    % legend('Physical dose', 'Uniform RBE weighted dose', 'Variable RBE weighted dose', 'Location', 'Northwest')
    % print(fig, 'prueba', '-dpng')
    %
    %
    clearvars -except ct phantomtype cst pln stf ResultRBEMCN ResultRBEUCM ResultConstRBE ResultPhysical GraphSel DoseRecalc StatsRef CompDVH
end


%% 6 - Graficas de perfil de dosis vs RBE

if GraphSel(2) > 0
    % ProfileType = longitudinal // lateral
    % DisplayOption = RBE // physicalDose // RBExD // physical_vs_RBExD
    % Para hacer comparaciones entre modelos DisplayOption == RBExD // physical_vs_RBExD
    % prueba_DoseGraphs (ct, pln, cst, NumBeam, ProfileType, DisplayOption, Result, Model1, Result2, Model2)
    if strcmp(version('-release'),'2014b') == 0
        
        prueba_DoseGraphs (ct, pln, cst,1, 'longitudinal', 'RBE', ResultConstRBE.RBEMCNreCalc.resultGUI, 'RBE', ResultConstRBE.RBEMCNreCalc.resultGUI, 'RBEMCN')
        title('RBE vs RBExD (ConstRBE optimized)')
        
        prueba_DoseGraphs (ct, pln, cst,1, 'longitudinal', 'RBE', ResultRBEMCN.Optimized.resultGUI, 'RBE', ResultRBEMCN.Optimized.resultGUI, 'RBEMCN')
        title('RBE vs RBExD (RBEMCN optimized)')
        
        prueba_DoseGraphs (ct, pln, cst,1, 'longitudinal', 'RBE', ResultRBEUCM.RBEMCNreCalc.resultGUI, 'RBE', ResultRBEUCM.RBEMCNreCalc.resultGUI, 'RBEUCM')
        title('RBE vs RBExD (RBEUCM optimized)')
        
    end
end


%% 7 - Graficas de dosis 2D

if GraphSel(3) > 0
    % prueba_DoseIntens (ct, pln, Dose, IsoDose_Levels, z_cut, TypeDose, Model)
    
    IsoDose_Levels = [30 40 64.6 66.84]; %(Gy)
    corte = [];
    
    % physicalDose Optimized
    % prueba_DoseIntens (ct, pln, ResultPhysical.Optimized.resultGUI.physicalDose, IsoDose_Levels, corte, 'Physical', 'Physical');
    % prueba_DoseIntens (ct, pln, ResultPhysical.ConstRBEreCalc.resultGUI.RBExD, IsoDose_Levels, corte, 'RBExD', 'ConstRBE');
    % prueba_DoseIntens (ct, pln, ResultPhysical.RBEMCNreCalc.resultGUI.RBExD, IsoDose_Levels, corte, 'RBExD', 'RBEMCN');
    
    % ConstRBE Optimized
    %prueba_DoseIntens (ct, pln, ResultConstRBE.Optimized.resultGUI.physicalDose, IsoDose_Levels, corte, 'Physical', 'Physical');
    prueba_DoseIntens (ct, pln, ResultConstRBE.Optimized.resultGUI.RBExD, IsoDose_Levels, corte, 'RBExD', 'ConstRBE');
    prueba_DoseIntens (ct, pln, ResultConstRBE.RBEMCNreCalc.resultGUI.RBExD, IsoDose_Levels, corte, 'RBExD', 'RBEMCN');
    
    % MCN Optimized
    %prueba_DoseIntens (ct, pln, ResultRBEMCN.Optimized.resultGUI.physicalDose, IsoDose_Levels, corte, 'Physical', 'Physical');
    prueba_DoseIntens (ct, pln, ResultRBEMCN.ConstRBEreCalc.resultGUI.RBExD, IsoDose_Levels, corte, 'RBExD', 'ConstRBE');
    prueba_DoseIntens (ct, pln, ResultRBEMCN.Optimized.resultGUI.RBExD, IsoDose_Levels, corte, 'RBExD', 'RBEMCN');
    
    % UCM Optimized
    %prueba_DoseIntens (ct, pln, ResultRBEUCM.Optimized.resultGUI.physicalDose,IsoDose_Levels, corte, 'Physical', 'Physical');
    prueba_DoseIntens (ct, pln, ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD, IsoDose_Levels, corte, 'RBExD', 'ConstRBE');
    prueba_DoseIntens (ct, pln, ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD, IsoDose_Levels, corte, 'RBExD', 'RBEMCN');
    
    clearvars -except ct phantomtype cst pln stf ResultRBEMCN ResultRBEUCM ResultConstRBE ResultPhysical GraphSel DoseRecalc StatsRef CompDVH
end

%% 8 - Representacion de DVH para todas las VOI

if GraphSel(4) == 1
    %Comparacion para RBE constante optimizado
    OptModel = cell(3,2);
    OptModel(1,1) = {'ConstRBE Optimized'};
    OptModel(1,2) = {'ConstRBE'};
    OptModel(2,2) = {'RBEMCN'};
    
    Dose{1,1} = ResultConstRBE.Optimized.resultGUI.RBExD;
    Dose{2,1} = ResultConstRBE.RBEMCNreCalc.resultGUI.RBExD;
    
    prueba_compDVH ( pln , cst, Dose, [], OptModel);
    clear OptModel
    
    %Comparacion para el modelo de McNamara optimizado
    OptModel = cell(3,2);
    OptModel(1,1) = {'RBEMCN Optimized'};
    OptModel(1,2) = {'ConstRBE'};
    OptModel(2,2) = {'RBEMCN'};
    
    Dose{1,1} = ResultRBEMCN.ConstRBEreCalc.resultGUI.RBExD;
    Dose{2,1} = ResultRBEMCN.Optimized.resultGUI.RBExD;
    
    prueba_compDVH ( pln , cst, Dose, [], OptModel);
    clear OptModel
    
    %Comparacion para el modelo UCM optimizado
    OptModel = cell(3,2);
    OptModel(1,1) = {'RBEUCM Optimized'};
    OptModel(1,2) = {'ConstRBE'};
    OptModel(2,2) = {'RBEMCN'};
    
    Dose{1,1} = ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD;
    Dose{2,1} = ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD;
    
    prueba_compDVH ( pln , cst, Dose, [], OptModel);
    clear OptModel
    
    
    clearvars -except ct phantomtype cst pln stf ResultRBEMCN ResultRBEUCM ResultConstRBE ResultPhysical GraphSel DoseRecalc StatsRef CompDVH
end 

%% 9 - Calculo de RBE medio para RBEMCN con las diversas optimizaciones

% clear midRBE
% if strcmp(phantomtype, 'Prostate') >0
%     Regions{1,1} = 'Rectum';
%     Regions{2,1} = 'PTV_68';
%     
% elseif strcmp(phantomtype, 'Head and Neck') >0
%     Regions = [];
% end
% 
% 
% Dose{1,1} = ResultConstRBE.RBEMCNreCalc.resultGUI;
% Dose{2,1} = ResultRBEMCN.Optimized.resultGUI;
% Dose{3,1} = ResultRBEUCM.RBEMCNreCalc.resultGUI;
% midRBE(1).OptModel = 'ConstRBE';
% midRBE(2).OptModel = 'RBEMCN';
% midRBE(3).OptModel = 'RBEUCM';
% 
% for i = 1:size(Regions,1)
%     for j = 1: size (Dose,1)
%         for k = 1:size(cst,1)
%             if strcmp (Regions{i,1},cst{k,2}) > 0
%                 indices     = cst{k,4}{1};
%                 % midRBE point by point
%                 RBE_pbp = Dose{j,1}.RBExD(indices) ./ Dose{j,1}.physicalDose(indices);
%                 RBE_pbp(isnan(RBE_pbp)>0) = 1.1;
%                 
%                 % midRBE sum points
%                 RBE_sum = sum(Dose{j,1}.RBExD(indices)) / sum(Dose{j,1}.physicalDose(indices));
%                 
%                 % midRBE over threshold
%                 RBE_thr = Dose{j,1}.RBExD(indices) ./ Dose{j,1}.physicalDose(indices);
%                 phyDose = Dose{j,1}.physicalDose(indices);
%                 thr_dose = 0.5;
%                 thr_index = phyDose < thr_dose;
%                 RBE_thr(thr_index) = 1.1;
%                 RBE_thr(isnan(RBE_thr)>0) = 1.1;
%                 
%                 if strcmpi (Regions{i,1}(1:3),'PTV') > 0
%                     midRBE(j).(strcat('PTV',Regions{i,1}(5:6),'_pbp')) = sum(RBE_pbp) / numel(indices);
%                     midRBE(j).(strcat('PTV',Regions{i,1}(5:6),'_sum')) = RBE_sum;
%                     midRBE(j).(strcat('PTV',Regions{i,1}(5:6),'_thr')) = sum(RBE_thr) / numel(indices);
%                     
%                     
%                 else
%                     midRBE(j).(strcat(Regions{i,1},'_pbp')) = sum(RBE_pbp) / numel(indices);
%                     midRBE(j).(strcat(Regions{i,1},'_sum')) = RBE_sum;
%                     midRBE(j).(strcat(Regions{i,1},'_thr')) = sum(RBE_thr) / numel(indices);
%                 end
%                 
%             end
%         end
%     end
% end
% 
% clearvars -except ct phantomtype cst pln stf ResultRBEMCN ResultRBEUCM ResultConstRBE midRBE ResultPhysical GraphSel DoseRecalc StatsRef
% ResultConstRBE ResultPhysical midRBE GraphSel DoseRecalc


%% 10 - Comparacion de DVH para VOIs especificas

if GraphSel(4) == 2 && ~isempty(CompDVH{1,1}) > 0
    
    Regions = CompDVH{1,1};
    VOIType = CompDVH{2,1};
       
    Dose{1,1} = ResultConstRBE.Optimized.resultGUI.RBExD;
    Dose{2,1} = ResultConstRBE.RBEMCNreCalc.resultGUI.RBExD;
    Dose{3,1} = ResultRBEMCN.Optimized.resultGUI.RBExD;
    Dose{4,1} = ResultRBEMCN.ConstRBEreCalc.resultGUI.RBExD;
    Dose{5,1} = ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD;
    Dose{6,1} = ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD;
    
    OptModel = cell(size(Dose,1),2);
    
    
    for i = 1:size(Dose,1)
        if Dose{i,1} == ResultConstRBE.Optimized.resultGUI.RBExD
            OptModel(i,1) = {' ConstRBEOpt '};
            OptModel(i,2) = {' ConstRBE '};
        elseif Dose{i,1} == ResultConstRBE.RBEMCNreCalc.resultGUI.RBExD
            OptModel(i,1) = {' ConstRBEOpt '};
            OptModel(i,2) = {' RBEMCN '};
        elseif Dose{i,1} == ResultRBEMCN.ConstRBEreCalc.resultGUI.RBExD
            OptModel(i,1) = {' RBEMCNOpt '};
            OptModel(i,2) = {' ConstRBE '};
        elseif  Dose{i,1} == ResultRBEMCN.Optimized.resultGUI.RBExD
            OptModel(i,1) = {' RBEMCNOpt '};
            OptModel(i,2) = {' RBEMCN '};
        elseif Dose{i,1} == ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD
            OptModel(i,1) = {' RBEUCMOpt '};
            OptModel(i,2) = {' ConstRBE '};
        elseif Dose{i,1} == ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD
            OptModel(i,1) = {' RBEUCMOpt '};
            OptModel(i,2) = {' RBEMCN '};
        end
    end
    
    prueba_compDVH (pln , cst, Dose, Regions, OptModel, VOIType);
    %saveas(image_All, 'DVHComp_All.png'); close;
    clear OptModel Dose
else
    fprintf('Warning: There are no selected regions. Skipping DVH comparison for specific VOIs');
end

clearvars -except ct phantomtype cst pln stf ResultRBEMCN ResultRBEUCM ResultConstRBE midRBE ResultPhysical GraphSel DoseRecalc StatsRef


%% 11 - Calculo de las estadisticas generales de dosis
DoseStatistics = 'Not Evaluated';

if exist('Renorm','var')>0
    if strcmp(Renorm,'Renormalized')>0
        Renorm = 'Renormalized';
    elseif DoseRecalc{1,1}{1,1} > 0 && DoseRecalc{2,1}{1,1} > 0 && DoseRecalc{3,1}{1,1} > 0
        Renorm = 'Not renormalized';
    elseif strcmp(Renorm,'Not renormalized')>0
        Renorm = 'Not renormalized';
    else
        Renorm = 'Not all result are renormalized';
    end
else
    Renorm = 'Not renormalized';
end

if GraphSel(5) > 0  
        
    if isempty(StatsRef)
        refVol = [];
        refGy = [];
    else
        refVol = StatsRef{1,1};
        refGy = StatsRef{2,1};
    end
    
    %prueba_DVHstatsComp (pln, cst, Dose, refVol, refGy, Models, Name, FigRem)
    
    clear DoseStatistics
    
    % Estadisticas de dosis para el caso de ConstRBE optimizado
    Dose{1,1} = ResultConstRBE.Optimized.resultGUI.RBExD;
    Dose{2,1} = ResultConstRBE.RBEMCNreCalc.resultGUI.RBExD;
    Models{1,1} = 'ConstRBE'; Models{1,2} = 'Optimized';
    Models{2,1} = 'RBEMCN';Models{2,2} = 'RBEMCNreCalc';
    
    DoseStatistics.ConstRBEOpt = prueba_DVHstatsComp(pln, cst, Dose,refVol,refGy, Models, 'ConstRBEOpt',0);
    clear Dose Models
    
    % Estadisticas de dosis para el modelo de McNamara optimizado
    Dose{1,1} = ResultRBEMCN.ConstRBEreCalc.resultGUI.RBExD;
    Dose{2,1} = ResultRBEMCN.Optimized.resultGUI.RBExD;
    Models{1,1} = 'ConstRBE'; Models{1,2} = 'ConstRBEreCalc';
    Models{2,1} = 'RBEMCN';Models{2,2} = 'Optimized';
    
    DoseStatistics.RBEMCNOpt = prueba_DVHstatsComp(pln, cst, Dose,refVol,refGy, Models, 'RBEMCNOpt',0);
    clear Dose Models
    
    % Estadisticas de dosis para el modelo UCM optimizado
    Dose{1,1} = ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD;
    Dose{2,1} = ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD;
    %Dose{3,1} = ResultRBEUCM.Optimized.resultGUI.RBExD;
    Models{1,1} = 'ConstRBE'; Models{1,2} = 'ConstRBEreCalc';
    Models{2,1} = 'RBEMCN';Models{2,2} = 'RBEMCNreCalc';
    %Models{3,1} = 'RBEUCM';Models{3,2} = 'Optimized';
    
    DoseStatistics.RBEUCMOpt = prueba_DVHstatsComp(pln, cst, Dose,refVol,refGy, Models, 'RBEUCMOpt',0);
    clear Dose Models    
    
    %% 12 - Renormalizacion y recalculo de estadisticas
    
  if GraphSel(5) == 2 && exist('DoseStatistics', 'var')
      
    Renorm = 'Renormalized';
    
    meanDoseTargetUCM = DoseStatistics.RBEUCMOpt.ConstRBEreCalc(6).mean;
    meanDoseTargetConst = DoseStatistics.ConstRBEOpt.Optimized(6).mean;
    meanDoseTargetMCN = DoseStatistics.RBEMCNOpt.ConstRBEreCalc(6).mean;
    
    if strcmpi(phantomtype, 'Head and Neck') > 0
        meanDoseTargetPrescription = 70;
    elseif strcmpi(phantomtype, 'Prostate') > 0
        meanDoseTargetPrescription = 78;
    elseif strcmpi(phantomtype, 'TG119') > 0
        meanDoseTargetPrescription = 55;
    end 
    
    FnormUCM = meanDoseTargetPrescription / meanDoseTargetUCM;
    FnormConst = meanDoseTargetPrescription / meanDoseTargetConst;
    FnormMCN = meanDoseTargetPrescription / meanDoseTargetMCN;
    
    ResultConstRBE.Optimized.resultGUI.RBExD = FnormConst * ResultConstRBE.Optimized.resultGUI.RBExD;
    ResultConstRBE.RBEMCNreCalc.resultGUI.RBExD = FnormConst * ResultConstRBE.RBEMCNreCalc.resultGUI.RBExD;
    ResultRBEMCN.ConstRBEreCalc.resultGUI.RBExD = FnormMCN * ResultRBEMCN.ConstRBEreCalc.resultGUI.RBExD;
    ResultRBEMCN.Optimized.resultGUI.RBExD = FnormMCN * ResultRBEMCN.Optimized.resultGUI.RBExD;
    ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD = FnormUCM * ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD;
    ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD = FnormUCM * ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD;

     % Estadisticas de dosis para el caso de ConstRBE optimizado
    Dose{1,1} = ResultConstRBE.Optimized.resultGUI.RBExD;
    Dose{2,1} = ResultConstRBE.RBEMCNreCalc.resultGUI.RBExD;
    Models{1,1} = 'ConstRBE'; Models{1,2} = 'Optimized';
    Models{2,1} = 'RBEMCN';Models{2,2} = 'RBEMCNreCalc';
    
    DoseStatistics.ConstRBEOpt = prueba_DVHstatsComp(pln, cst, Dose,refVol,refGy, Models, 'ConstRBEOpt',0);
    clear Dose Models
    
    % Estadisticas de dosis para el modelo de McNamara optimizado
    Dose{1,1} = ResultRBEMCN.ConstRBEreCalc.resultGUI.RBExD;
    Dose{2,1} = ResultRBEMCN.Optimized.resultGUI.RBExD;
    Models{1,1} = 'ConstRBE'; Models{1,2} = 'ConstRBEreCalc';
    Models{2,1} = 'RBEMCN';Models{2,2} = 'Optimized';
    
    DoseStatistics.RBEMCNOpt = prueba_DVHstatsComp(pln, cst, Dose,refVol,refGy, Models, 'RBEMCNOpt',0);
    clear Dose Models
    
    % Estadisticas de dosis para el modelo UCM optimizado
    Dose{1,1} = ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD;
    Dose{2,1} = ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD;
    %Dose{3,1} = ResultRBEUCM.Optimized.resultGUI.RBExD;
    Models{1,1} = 'ConstRBE'; Models{1,2} = 'ConstRBEreCalc';
    Models{2,1} = 'RBEMCN';Models{2,2} = 'RBEMCNreCalc';
    %Models{3,1} = 'RBEUCM';Models{3,2} = 'Optimized';
    
    DoseStatistics.RBEUCMOpt = prueba_DVHstatsComp(pln, cst, Dose,refVol,refGy, Models, 'RBEUCMOpt',0);
    clear Dose Models 
  end
    clearvars -except ct cst phantomtype pln stf ResultRBEMCN ResultRBEUCM ResultConstRBE ResultPhysical midRBE  DoseStatistics GraphSel DoseRecalc StatsRef Renorm
    
end


%% 13 - Calculo de NTCP

if strcmpi(phantomtype, 'TG119') > 0
    NTCP = 'No models';
else
    clear NTCP
    
    % % Calculos de NTCP para PhysicalDose
    % NTCP.PhysicalOpt.Optimized = prueba_NTCPcalc(pln, cst, phantomtype, ResultPhysical.Optimized.resultGUI.physicalDose);
    % NTCP.PhysicalOpt.ConstRBEreCalc = prueba_NTCPcalc(pln, cst, phantomtype, ResultPhysical.ConstRBEreCalc.resultGUI.RBExD);
    % NTCP.PhysicalOpt.RBEMCNreCalc = prueba_NTCPcalc(pln, cst, phantomtype, ResultPhysical.RBEMCNreCalc.resultGUI.RBExD);
    
    % Calculos NTCP para ConstRBEOpt
    NTCP.ConstRBEOpt.Optimized = prueba_NTCPcalc(pln, cst, phantomtype, ResultConstRBE.Optimized.resultGUI.RBExD);
    NTCP.ConstRBEOpt.RBEMCNreCalc = prueba_NTCPcalc(pln, cst, phantomtype, ResultConstRBE.RBEMCNreCalc.resultGUI.RBExD);
    
    % Calculos NTCP para RBEMCNOpt
    NTCP.RBEMCNOpt.ConstRBEreCalc = prueba_NTCPcalc(pln, cst, phantomtype, ResultRBEMCN.ConstRBEreCalc.resultGUI.RBExD);
    NTCP.RBEMCNOpt.Optimized = prueba_NTCPcalc(pln, cst, phantomtype, ResultRBEMCN.Optimized.resultGUI.RBExD);
    
    % Calculos NTCP para RBEUCMOpt
    NTCP.RBEUCMOpt.ConstRBEreCalc = prueba_NTCPcalc(pln, cst, phantomtype, ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD);
    NTCP.RBEUCMOpt.RBEMCNreCalc = prueba_NTCPcalc(pln, cst, phantomtype, ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD);
    
    clearvars -except ct cst phantomtype pln stf ResultRBEMCN ResultRBEUCM ResultConstRBE ResultPhysical midRBE  DoseStatistics NTCP GraphSel DoseRecalc StatsRef Renorm
    
end


%% 14 - "Resumen" de valores de NTCP

clear meanNTCP

if strcmp (phantomtype, 'Prostate') > 0
    NTCP_bio_phys = nan(9,1);
    NTCP_bio_const = nan(9,1);
    NTCP_bio_MCN = nan(9,1);
    NTCP_bio_UCM = nan(9,1);
    
    %     NTCP_bio_phys(1) = NTCP.PhysicalOpt.RBEMCNreCalc.Fukahori.NTCP(2).NTCP;
    %     NTCP_bio_phys(2) = NTCP.PhysicalOpt.RBEMCNreCalc.Burman.Rectum.NTCP;
    %     NTCP_bio_phys(3) = NTCP.PhysicalOpt.RBEMCNreCalc.Cheung.Rectum.NTCP(2).NTCP;
    %     NTCP_bio_phys(4) = NTCP.PhysicalOpt.RBEMCNreCalc.Liu.NTCP;
    %     NTCP_bio_phys(5) = NTCP.PhysicalOpt.RBEMCNreCalc.Tucker.NTCP;
    %     NTCP_bio_phys(6) = NTCP.PhysicalOpt.RBEMCNreCalc.Peeters.Bleeding.NTCP;
    %     NTCP_bio_phys(7) = NTCP.PhysicalOpt.RBEMCNreCalc.Peeters.Dep_Freq_Incr.NTCP;
    %     NTCP_bio_phys(8) = NTCP.PhysicalOpt.RBEMCNreCalc.Peeters.Fecal_Inc.NTCP;
    %     NTCP_bio_phys(9) = NTCP.PhysicalOpt.RBEMCNreCalc.Schaake.NTCP(2).NTCP;
    
    NTCP_bio_const(1) = NTCP.ConstRBEOpt.RBEMCNreCalc.Fukahori.NTCP(2).NTCP;
    NTCP_bio_const(2) = NTCP.ConstRBEOpt.RBEMCNreCalc.Burman.Rectum.NTCP;
    NTCP_bio_const(3) = NTCP.ConstRBEOpt.RBEMCNreCalc.Cheung.Rectum.NTCP(2).NTCP;
    NTCP_bio_const(4) = NTCP.ConstRBEOpt.RBEMCNreCalc.Liu.NTCP;
    NTCP_bio_const(5) = NTCP.ConstRBEOpt.RBEMCNreCalc.Tucker.NTCP;
    NTCP_bio_const(6) = NTCP.ConstRBEOpt.RBEMCNreCalc.Peeters.Bleeding.NTCP;
    NTCP_bio_const(7) = NTCP.ConstRBEOpt.RBEMCNreCalc.Peeters.Dep_Freq_Incr.NTCP;
    NTCP_bio_const(8) = NTCP.ConstRBEOpt.RBEMCNreCalc.Peeters.Fecal_Inc.NTCP;
    NTCP_bio_const(9) = NTCP.ConstRBEOpt.RBEMCNreCalc.Schaake.NTCP(2).NTCP;
    
    NTCP_bio_MCN(1) = NTCP.RBEMCNOpt.Optimized.Fukahori.NTCP(2).NTCP;
    NTCP_bio_MCN(2) = NTCP.RBEMCNOpt.Optimized.Burman.Rectum.NTCP;
    NTCP_bio_MCN(3) = NTCP.RBEMCNOpt.Optimized.Cheung.Rectum.NTCP(2).NTCP;
    NTCP_bio_MCN(4) = NTCP.RBEMCNOpt.Optimized.Liu.NTCP;
    NTCP_bio_MCN(5) = NTCP.RBEMCNOpt.Optimized.Tucker.NTCP;
    NTCP_bio_MCN(6) = NTCP.RBEMCNOpt.Optimized.Peeters.Bleeding.NTCP;
    NTCP_bio_MCN(7) = NTCP.RBEMCNOpt.Optimized.Peeters.Dep_Freq_Incr.NTCP;
    NTCP_bio_MCN(8) = NTCP.RBEMCNOpt.Optimized.Peeters.Fecal_Inc.NTCP;
    NTCP_bio_MCN(9) = NTCP.RBEMCNOpt.Optimized.Schaake.NTCP(2).NTCP;
    
    NTCP_bio_UCM(1) = NTCP.RBEUCMOpt.RBEMCNreCalc.Fukahori.NTCP(2).NTCP;
    NTCP_bio_UCM(2) = NTCP.RBEUCMOpt.RBEMCNreCalc.Burman.Rectum.NTCP;
    NTCP_bio_UCM(3) = NTCP.RBEUCMOpt.RBEMCNreCalc.Cheung.Rectum.NTCP(2).NTCP;
    NTCP_bio_UCM(4) = NTCP.RBEUCMOpt.RBEMCNreCalc.Liu.NTCP;
    NTCP_bio_UCM(5) = NTCP.RBEUCMOpt.RBEMCNreCalc.Tucker.NTCP;
    NTCP_bio_UCM(6) = NTCP.RBEUCMOpt.RBEMCNreCalc.Peeters.Bleeding.NTCP;
    NTCP_bio_UCM(7) = NTCP.RBEUCMOpt.RBEMCNreCalc.Peeters.Dep_Freq_Incr.NTCP;
    NTCP_bio_UCM(8) = NTCP.RBEUCMOpt.RBEMCNreCalc.Peeters.Fecal_Inc.NTCP;
    NTCP_bio_UCM(9) = NTCP.RBEUCMOpt.RBEMCNreCalc.Schaake.NTCP(2).NTCP;
    
    % Mostrar valores quitando Cheung y Fukahori
    mascara = [4 5 6 7 8 9];
    
    NTCPMCNall = [ NTCP_bio_phys(mascara), NTCP_bio_const(mascara), NTCP_bio_MCN(mascara), NTCP_bio_UCM(mascara)]
    meanNTCP = [mean(NTCP_bio_phys(mascara)),mean(NTCP_bio_const(mascara)), mean(NTCP_bio_MCN(mascara)), mean(NTCP_bio_UCM(mascara))]
    
elseif strcmpi(phantomtype, 'Head and Neck') > 0
    
    NTCP_bio_phys = nan(14,1);
    NTCP_bio_const = nan(14,1);
    NTCP_bio_MCN = nan(14,1);
    NTCP_bio_UCM = nan(14,1);
    
    %     NTCP_bio_phys(1) = NTCP.PhysicalOpt.RBEMCNreCalc.Semenenko.NTCP_Left;
    %     NTCP_bio_phys(2) = NTCP.PhysicalOpt.RBEMCNreCalc.Burman.NTCP_Left;
    %     NTCP_bio_phys(3) = NTCP.PhysicalOpt.RBEMCNreCalc.Eisbruch.NTCP_Left;
    %     NTCP_bio_phys(4) = NTCP.PhysicalOpt.RBEMCNreCalc.Luxton.Larynx.Necro.NTCP;
    %     NTCP_bio_phys(5) = NTCP.PhysicalOpt.RBEMCNreCalc.Luxton.Larynx.Edema.NTCP;
    %     NTCP_bio_phys(6) = NTCP.PhysicalOpt.RBEMCNreCalc.Luxton.OcularLens.NTCP_Left;
    %     NTCP_bio_phys(7) = NTCP.PhysicalOpt.RBEMCNreCalc.Luxton.Parotid.NTCP_Left;
    %     NTCP_bio_phys(8) = NTCP.PhysicalOpt.RBEMCNreCalc.Luxton.Skin.NTCP;
    %     NTCP_bio_phys(9) = NTCP.PhysicalOpt.RBEMCNreCalc.Luxton.Spinal_Cord.NTCP;
    %     NTCP_bio_phys(10) = NTCP.PhysicalOpt.RBEMCNreCalc.Luxton.Tm_Joint.NTCP_Left;
    %     NTCP_bio_phys(11) = NTCP.PhysicalOpt.RBEMCNreCalc.Roesink.Left.NTCP_SEF25;
    %     NTCP_bio_phys(12) = NTCP.PhysicalOpt.RBEMCNreCalc.Roesink.Left.NTCP_SEF35;
    %     NTCP_bio_phys(13) = NTCP.PhysicalOpt.RBEMCNreCalc.Roesink.Left.NTCP_SEF45;
    %     NTCP_bio_phys(14) = NTCP.PhysicalOpt.RBEMCNreCalc.Roesink.Left.NTCP_SEF55;
    
    NTCP_bio_const(1) = NTCP.ConstRBEOpt.RBEMCNreCalc.Semenenko.NTCP_Left;
    NTCP_bio_const(2) = NTCP.ConstRBEOpt.RBEMCNreCalc.Burman.NTCP_Left;
    NTCP_bio_const(3) = NTCP.ConstRBEOpt.RBEMCNreCalc.Eisbruch.NTCP_Left;
    NTCP_bio_const(4) = NTCP.ConstRBEOpt.RBEMCNreCalc.Luxton.Larynx.Necro.NTCP;
    NTCP_bio_const(5) = NTCP.ConstRBEOpt.RBEMCNreCalc.Luxton.Larynx.Edema.NTCP;
    NTCP_bio_const(6) = NTCP.ConstRBEOpt.RBEMCNreCalc.Luxton.OcularLens.NTCP_Left;
    NTCP_bio_const(7) = NTCP.ConstRBEOpt.RBEMCNreCalc.Luxton.Parotid.NTCP_Left;
    NTCP_bio_const(8) = NTCP.ConstRBEOpt.RBEMCNreCalc.Luxton.Skin.NTCP;
    NTCP_bio_const(9) = NTCP.ConstRBEOpt.RBEMCNreCalc.Luxton.Spinal_Cord.NTCP;
    NTCP_bio_const(10) = NTCP.ConstRBEOpt.RBEMCNreCalc.Luxton.Tm_Joint.NTCP_Left;
    NTCP_bio_const(11) = NTCP.ConstRBEOpt.RBEMCNreCalc.Roesink.Left.NTCP_SEF25;
    NTCP_bio_const(12) = NTCP.ConstRBEOpt.RBEMCNreCalc.Roesink.Left.NTCP_SEF35;
    NTCP_bio_const(13) = NTCP.ConstRBEOpt.RBEMCNreCalc.Roesink.Left.NTCP_SEF45;
    NTCP_bio_const(14) = NTCP.ConstRBEOpt.RBEMCNreCalc.Roesink.Left.NTCP_SEF55;
    
    NTCP_bio_MCN(1) = NTCP.RBEMCNOpt.Optimized.Semenenko.NTCP_Left;
    NTCP_bio_MCN(2) = NTCP.RBEMCNOpt.Optimized.Burman.NTCP_Left;
    NTCP_bio_MCN(3) = NTCP.RBEMCNOpt.Optimized.Eisbruch.NTCP_Left;
    NTCP_bio_MCN(4) = NTCP.RBEMCNOpt.Optimized.Luxton.Larynx.Necro.NTCP;
    NTCP_bio_MCN(5) = NTCP.RBEMCNOpt.Optimized.Luxton.Larynx.Edema.NTCP;
    NTCP_bio_MCN(6) = NTCP.RBEMCNOpt.Optimized.Luxton.OcularLens.NTCP_Left;
    NTCP_bio_MCN(7) = NTCP.RBEMCNOpt.Optimized.Luxton.Parotid.NTCP_Left;
    NTCP_bio_MCN(8) = NTCP.RBEMCNOpt.Optimized.Luxton.Skin.NTCP;
    NTCP_bio_MCN(9) = NTCP.RBEMCNOpt.Optimized.Luxton.Spinal_Cord.NTCP;
    NTCP_bio_MCN(10) = NTCP.RBEMCNOpt.Optimized.Luxton.Tm_Joint.NTCP_Left;
    NTCP_bio_MCN(11) = NTCP.RBEMCNOpt.Optimized.Roesink.Left.NTCP_SEF25;
    NTCP_bio_MCN(12) = NTCP.RBEMCNOpt.Optimized.Roesink.Left.NTCP_SEF35;
    NTCP_bio_MCN(13) = NTCP.RBEMCNOpt.Optimized.Roesink.Left.NTCP_SEF45;
    NTCP_bio_MCN(14) = NTCP.RBEMCNOpt.Optimized.Roesink.Left.NTCP_SEF55;
    
    
    NTCP_bio_UCM(1) = NTCP.RBEUCMOpt.RBEMCNreCalc.Semenenko.NTCP_Left;
    NTCP_bio_UCM(2) = NTCP.RBEUCMOpt.RBEMCNreCalc.Burman.NTCP_Left;
    NTCP_bio_UCM(3) = NTCP.RBEUCMOpt.RBEMCNreCalc.Eisbruch.NTCP_Left;
    NTCP_bio_UCM(4) = NTCP.RBEUCMOpt.RBEMCNreCalc.Luxton.Larynx.Necro.NTCP;
    NTCP_bio_UCM(5) = NTCP.RBEUCMOpt.RBEMCNreCalc.Luxton.Larynx.Edema.NTCP;
    NTCP_bio_UCM(6) = NTCP.RBEUCMOpt.RBEMCNreCalc.Luxton.OcularLens.NTCP_Left;
    NTCP_bio_UCM(7) = NTCP.RBEUCMOpt.RBEMCNreCalc.Luxton.Parotid.NTCP_Left;
    NTCP_bio_UCM(8) = NTCP.RBEUCMOpt.RBEMCNreCalc.Luxton.Skin.NTCP;
    NTCP_bio_UCM(9) = NTCP.RBEUCMOpt.RBEMCNreCalc.Luxton.Spinal_Cord.NTCP;
    NTCP_bio_UCM(10) = NTCP.RBEUCMOpt.RBEMCNreCalc.Luxton.Tm_Joint.NTCP_Left;
    NTCP_bio_UCM(11) = NTCP.RBEUCMOpt.RBEMCNreCalc.Roesink.Left.NTCP_SEF25;
    NTCP_bio_UCM(12) = NTCP.RBEUCMOpt.RBEMCNreCalc.Roesink.Left.NTCP_SEF35;
    NTCP_bio_UCM(13) = NTCP.RBEUCMOpt.RBEMCNreCalc.Roesink.Left.NTCP_SEF45;
    NTCP_bio_UCM(14) = NTCP.RBEUCMOpt.RBEMCNreCalc.Roesink.Left.NTCP_SEF55;
    
    
    mascara = [1 2 3 4 5 6 7 8 9 10 11 12 13 14];
    
    NTCPMCNall = [ NTCP_bio_phys(mascara), NTCP_bio_const(mascara), NTCP_bio_MCN(mascara), NTCP_bio_UCM(mascara)]
    meanNTCP = [mean(NTCP_bio_phys(mascara)),mean(NTCP_bio_const(mascara)), mean(NTCP_bio_MCN(mascara)), mean(NTCP_bio_UCM(mascara))]
    
    elseif strcmpi(phantomtype, 'TG119') > 0
         meanNTCP = 'No models';  
end


end