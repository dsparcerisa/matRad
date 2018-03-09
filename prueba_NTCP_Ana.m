%% Estructuracion del programa

% 1 - Carga del phantom/paciente
% 2 - Carga de parametros alpha y beta
% 3 - Definicion de restricciones "distintas"
% 4 - Introduccion de los datos basicos
% 5 - Calculo y optimizacion para dosis fisica
% 6 - Calculo y optimizacion de dosis considerando el RBE = 1.1
% 7 - Calculo y optimizacion de dosis para el modelo McNamara
% 8 - Calculo y optimizacion de dosis para el modelo de planificacion UCM
% 9 - Graficas de perfil de dosis
% 10 - Graficas de perfil de dosis vs RBE
% 11 - Graficas de dosis 2D
% 12 - Representacion de DVH para todas las VOI
% 13 - Calculo de RBE medio para RBEMCN con las diversas optimizaciones
% 14 - Comparacion de DVH para VOIs especificas
% 15 - Calculo de las estadisticas generales de dosis
% 16 - Calculo de NTCP
% 17 - "Resumen" de valores de NTCP
% 18 - Resultados de dij y resultGUI para ver en la GUI

%% 1 - Carga del phantom/paciente

close all
% clear
if ~ismac
    opengl software
end
load PROSTATE.mat; phantomtype = 'Prostate';
% %load HEAD_AND_neck.mat; phantomtype = 'Head and Neck';
% %load BOXPHANTOM.mat; phantomtype = 'Test';
% %load TG119.mat; phantomtype = 'Test';

%phantomtype = 'Prostate';

%% 2 - Carga de parametros alpha y beta

cst = prueba_abLoader (cst, phantomtype);

%% 3 - Definicion de restricciones "distintas"

%       .type                   || NaN parameters

% 'square underdosing'          ||.EUD - .volume
% 'square overdosing'           ||.EUD - .volume
% 'square deviation'            ||.EUD - .volume
% 'mean'                        ||.dose - .EUD - .volume
% 'EUD'                         ||.dose - .volume
% 'min dose constraint'         ||.penalty - .EUD - .volume
% 'max dose constraint'         ||.penalty - .EUD - .volume
% 'min mean dose constraint'    ||.EUD - .volume
% 'max mean dose constraint'    ||.EUD - .volume
% 'min EUD constraint'          ||.penalty - .dose - .volume
% 'max EUD constraint'          ||.penalty - .dose - .volume
% 'min DVH constraint'          ||.penalty - .EUD 
% 'max DVH constraint'          ||.penalty - .EUD
% 'min DVH objective'           ||.EUD
% 'max DVH objective'           ||.EUD


% EUD para el Rectum
cst{1,6}.type = 'EUD';
cst{1,6}.dose = NaN;
cst{1,6}.EUD = 85; % Con 70 se vuelve to loco
cst{1,6}.penalty = 50;%jugar
cst{1,6}.volume = 50;
cst{1,6}.robustness = 'none';

% PTV_68
cst{6,6}(2).type = 'min DVH objective';
cst{6,6}(2).dose = 68 / 30; % por fraccion
cst{6,6}(2).penalty = 1; % probando
cst{6,6}(2).EUD = NaN;
cst{6,6}(2).volume = 1; % tanto por uno, parece
cst{6,6}(2).robustness = 'none';

% % PTV_68
cst{6,6}(1).type = 'square deviation';
cst{6,6}(1).dose = 68;
cst{6,6}(1).penalty = 1000;
cst{6,6}(1).EUD = NaN;
cst{6,6}(1).volume = NaN;
cst{6,6}(1).robustness = 'none';

% PTV_56
cst{7,6}.type = 'square deviation';
cst{7,6}.dose = 56;
cst{7,6}.penalty = 1000;
cst{7,6}.EUD = NaN;
cst{7,6}.volume = NaN;
cst{7,6}.robustness = 'none';

% Bladder
cst{8,6}.type = 'square overdosing';
cst{8,6}.dose = 50;
cst{8,6}.penalty = 300;
cst{8,6}.EUD = NaN;
cst{8,6}.volume = NaN;
cst{8,6}.robustness = 'none';


% Body 
cst{9,6}.type = 'square overdosing';
cst{9,6}.dose = 30;
cst{9,6}.penalty = 100;
cst{9,6}.EUD = NaN;
cst{9,6}.volume = NaN;
cst{9,6}.robustness = 'none';



%% 4 - Introduccion de los datos basicos

pln.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = [0]; % [??]
pln.couchAngles     = [0]; % [??]
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = prod(ct.cubeDim);
pln.isoCenter       = ones(pln.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.voxelDimensions = ct.cubeDim;
pln.radiationMode   = 'protons';             % either photons / protons / helium / carbon
pln.scenGenType     = 'nomScen';             % scenario creation type'nomScen'  'wcScen' 'impScen' 'rndScen'

pln.numOfFractions  = 30;
pln.runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.machine         = 'Generic';
pln.robOpt          = false;
%pln.calcLET = true;


%% 5 - Calculo y optimizacion para dosis fisica

quantityOpt         = 'physicalDose';     % options: physicalDose, constRBE, effect, RBExD
modelName           = 'none';             % none: for photons, protons, carbon                                    MCN: McNamara-variable RBE model for protons
                                          % WED: Wedenberg-variable RBE model for protons 

% retrieve model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,pln.scenGenType); 
% generate steering file
stf = matRad_generateStf(ct,cst,pln);

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

ResultPhysical.Optimized.dij = dij;
ResultPhysical.Optimized.resultGUI = resultGUI;


%% Recalculo de dosis para ConstRBE Y RBEMCN

[~ ,ResultPhysical.ConstRBEreCalc.resultGUI] = prueba_RecalcDose(ResultPhysical.Optimized.resultGUI, ct, stf, pln, cst, 'RBExD','constRBE');

[~ ,ResultPhysical.RBEMCNreCalc.resultGUI] = prueba_RecalcDose(ResultPhysical.Optimized.resultGUI, ct, stf, pln, cst, 'effect', 'MCN' );

clear dij resultGUI quantityOpt modelName


%% 6 - Calculo y optimizacion de dosis considerando el RBE = 1.1

quantityOpt         = 'RBExD';     % options: physicalDose, constRBE, effect, RBExD
modelName           = 'constRBE';             % none: for photons, protons, carbon                                    MCN: McNamara-variable RBE model for protons
                                          % WED: Wedenberg-variable RBE model for protons 
% retrieve model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode, quantityOpt, modelName);
% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,pln.scenGenType); 
% generate steering file
stf = matRad_generateStf(ct,cst,pln);

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


%% Recalculo de dosis para RBEMCN

[~ ,ResultConstRBE.RBEMCNreCalc.resultGUI] = prueba_RecalcDose(resultGUI, ct, stf, pln, cst, 'effect', 'MCN');

clear dij resultGUI quantityOpt modelName


%% 7 - Calculo y optimizacion de dosis para el modelo McNamara

quantityOpt         = 'RBExD';      % options: physicalDose, constRBE, effect, RBExD
modelName           = 'MCN';             % none: for photons, protons, carbon                                    MCN: McNamara-variable RBE model for protons
                                         % WED: Wedenberg-variable RBE model for protons
                                                
% retrieve model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);
pln.bioOpt = 1;

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,pln.scenGenType); 

stf = matRad_generateStf(ct,cst,pln);

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


%% Recalculo de dosis para ConstRBE

[~ ,ResultRBEMCN.ConstRBEreCalc.resultGUI] = prueba_RecalcDose(ResultRBEMCN.Optimized.resultGUI, ct, stf, pln, cst,'RBExD','constRBE');

clear dij resultGUI quantityOpt modelName


%% 8 - Calculo y optimizacion de dosis para el modelo de planificacion UCM
quantityOpt         = 'RBExD';
modelName           = 'UCM';             % none: for photons, protons, carbon                                    MCN: McNamara-variable RBE model for protons
                                          % WED: Wedenberg-variable RBE model for protons   
                                             
% retrieve model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,pln.scenGenType); 

stf = matRad_generateStf(ct,cst,pln);

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


%% Recalculo de dosis para ConstRBE Y RBEMCN

[~ ,ResultRBEUCM.ConstRBEreCalc.resultGUI] = prueba_RecalcDose(ResultRBEUCM.Optimized.resultGUI, ct, stf, pln, cst, 'RBExD','constRBE');

[~ ,ResultRBEUCM.RBEMCNreCalc.resultGUI] = prueba_RecalcDose(ResultRBEUCM.Optimized.resultGUI, ct, stf, pln, cst, 'effect', 'MCN');

clear dij resultGUI quantityOpt modelName


%% 9 - Graficas de perfil de dosis
hkhk % para que de error aqui y poder parar.

% ProfileType = longitudinal // lateral
% DisplayOption = physicalDose // RBExD // physical_vs_RBExD
% Para hacer comparaciones entre modelos DisplayOption == RBExD // physical_vs_RBExD
% prueba_DoseGraphs (ct, pln, cst, NumBeam, ProfileType, DisplayOption, Result, Model1, Result2, Model2)

prueba_DoseGraphs (ct, pln, cst,1, 'longitudinal', 'physical_vs_RBExD', ResultPhysical.ConstRBEreCalc.resultGUI, 'ConstRBE',ResultPhysical.RBEMCNreCalc.resultGUI,'RBEMCN')
title('Physical dose optimized')

prueba_DoseGraphs (ct, pln, cst,1, 'longitudinal', 'physical_vs_RBExD', ResultConstRBE.Optimized.resultGUI, 'ConstRBE', ResultConstRBE.RBEMCNreCalc.resultGUI, 'RBEMCN')
title('ConstRBE optimized')

prueba_DoseGraphs (ct, pln, cst,1, 'longitudinal', 'physical_vs_RBExD', ResultRBEMCN.Optimized.resultGUI, 'RBEMCN',ResultRBEMCN.ConstRBEreCalc.resultGUI,'ConstRBE')
title('RBEMCN optimized')

prueba_DoseGraphs (ct, pln, cst,1, 'longitudinal', 'physical_vs_RBExD', ResultRBEUCM.RBEMCNreCalc.resultGUI, 'RBEMCN',ResultRBEUCM.ConstRBEreCalc.resultGUI,'ConstRBE')
title('RBEUCM optimized')

clearvars -except ct phantomtype cst pln ResultRBEMCN ResultRBEUCM ResultConstRBE ResultPhysical


%% 10 - Graficas de perfil de dosis vs RBE
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
%% 11 - Graficas de dosis 2D 
% prueba_DoseIntens (ct, pln, Dose, IsoDose_Levels, z_cut, TypeDose, Model)

IsoDose_Levels = [10 20 30 40 50 60 68]; %(Gy)
corte = 30;
% physicalDose Optimized
prueba_DoseIntens (ct, pln, ResultPhysical.Optimized.resultGUI.physicalDose, IsoDose_Levels, corte, 'Physical', 'Physical');
prueba_DoseIntens (ct, pln, ResultPhysical.ConstRBEreCalc.resultGUI.RBExD, IsoDose_Levels, corte, 'RBExD', 'ConstRBE');
prueba_DoseIntens (ct, pln, ResultPhysical.RBEMCNreCalc.resultGUI.RBExD, IsoDose_Levels, corte, 'RBExD', 'RBEMCN');

% ConstRBE Optimized
prueba_DoseIntens (ct, pln, ResultConstRBE.Optimized.resultGUI.physicalDose, IsoDose_Levels, corte, 'Physical', 'Physical');
prueba_DoseIntens (ct, pln, ResultConstRBE.Optimized.resultGUI.RBExD, IsoDose_Levels, corte, 'RBExD', 'ConstRBE');
prueba_DoseIntens (ct, pln, ResultConstRBE.RBEMCNreCalc.resultGUI.RBExD, IsoDose_Levels, corte, 'RBExD', 'RBEMCN');

% MCN Optimized
prueba_DoseIntens (ct, pln, ResultRBEMCN.Optimized.resultGUI.physicalDose, IsoDose_Levels, corte, 'Physical', 'Physical');
prueba_DoseIntens (ct, pln, ResultRBEMCN.ConstRBEreCalc.resultGUI.RBExD, IsoDose_Levels, corte, 'RBExD', 'ConstRBE');
prueba_DoseIntens (ct, pln, ResultRBEMCN.Optimized.resultGUI.RBExD, IsoDose_Levels, corte, 'RBExD', 'RBEMCN');

% UCM Optimized
prueba_DoseIntens (ct, pln, ResultRBEUCM.Optimized.resultGUI.physicalDose,IsoDose_Levels, corte, 'Physical', 'Physical');
prueba_DoseIntens (ct, pln, ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD, IsoDose_Levels, corte, 'RBExD', 'ConstRBE');
prueba_DoseIntens (ct, pln, ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD, IsoDose_Levels, corte, 'RBExD', 'RBEMCN');

clearvars -except ct phantomtype cst pln ResultRBEMCN ResultRBEUCM ResultConstRBE ResultPhysical


%% 12 - Representacion de DVH para todas las VOI

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


clearvars -except ct phantomtype cst pln ResultRBEMCN ResultRBEUCM ResultConstRBE ResultPhysical


%% 13 - Calculo de RBE medio para RBEMCN con las diversas optimizaciones
if strcmp(phantomtype, 'Prostate') >0
    Regions{1,1} = 'Rectum';
    Regions{2,1} = 'PTV_68';
    
elseif strcmp(phantomtype, 'Head and Neck') >0
    Regions = [];
end


Dose{1,1} = ResultConstRBE.RBEMCNreCalc.resultGUI;
Dose{2,1} = ResultRBEMCN.Optimized.resultGUI;
Dose{3,1} = ResultRBEUCM.RBEMCNreCalc.resultGUI;
midRBE(1).OptModel = 'ConstRBE';
midRBE(2).OptModel = 'RBEMCN';
midRBE(3).OptModel = 'RBEUCM';

for i = 1:size(Regions,1)
    for j = 1: size (Dose,1)
        for k = 1:size(cst,1)
            if strcmp (Regions{i,1},cst{k,2}) > 0
                indices     = cst{k,4}{1};
                % midRBE point by point
                RBE_pbp = Dose{j,1}.RBExD(indices) ./ Dose{j,1}.physicalDose(indices);
                RBE_pbp(isnan(RBE_pbp)>0) = 1.1;
                
                % midRBE sum points
                RBE_sum = sum(Dose{j,1}.RBExD(indices)) / sum(Dose{j,1}.physicalDose(indices));
                
                % midRBE over threshold              
                RBE_thr = Dose{j,1}.RBExD(indices) ./ Dose{j,1}.physicalDose(indices);
                phyDose = Dose{j,1}.physicalDose(indices);
                thr_dose = 0.5;
                thr_index = phyDose < thr_dose;  
                RBE_thr(thr_index) = 1.1;
                RBE_thr(isnan(RBE_thr)>0) = 1.1;
                  
                if strcmpi (Regions{i,1}(1:3),'PTV') > 0 
                    midRBE(j).(strcat('PTV',Regions{i,1}(5:6),'_pbp')) = sum(RBE_pbp) / numel(indices);
                    midRBE(j).(strcat('PTV',Regions{i,1}(5:6),'_sum')) = RBE_sum;
                    midRBE(j).(strcat('PTV',Regions{i,1}(5:6),'_thr')) = sum(RBE_thr) / numel(indices);
                    

                else
                    midRBE(j).(strcat(Regions{i,1},'_pbp')) = sum(RBE_pbp) / numel(indices);
                    midRBE(j).(strcat(Regions{i,1},'_sum')) = RBE_sum;
                    midRBE(j).(strcat(Regions{i,1},'_thr')) = sum(RBE_thr) / numel(indices);
                end               
                
            end      
        end
    end
end

clearvars -except ct phantomtype cst pln ResultRBEMCN ResultRBEUCM ResultConstRBE ResultPhysical midRBE


%% 14 - Comparacion de DVH para VOIs especificas

% Todos los modelos juntos
if strcmp(phantomtype, 'Prostate') >0
    Regions{1,1} = 'Rectum';
    Regions{2,1} = 'PTV_68';    
elseif strcmp(phantomtype, 'Head and Neck') >0
    Regions = [];
end

OptModel = cell(6,2);
OptModel([1,2],1) = {' ConstRBEOpt '};
OptModel([3,4],1) = {' RBEMCNOpt '};
OptModel([5,6],1) = {' RBEUCMOpt '};
OptModel([1,3,5],2) = {' ConstRBE '};
OptModel([2,4,6],2) = {' RBEMCN '};

Dose{1,1} = ResultConstRBE.Optimized.resultGUI.RBExD;
Dose{2,1} = ResultConstRBE.RBEMCNreCalc.resultGUI.RBExD;
Dose{3,1} = ResultRBEMCN.ConstRBEreCalc.resultGUI.RBExD;
Dose{4,1} = ResultRBEMCN.Optimized.resultGUI.RBExD;
Dose{5,1} = ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD;
Dose{6,1}= ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD;

image_All = prueba_compDVH ( pln , cst, Dose, Regions, OptModel);
%saveas(image_All, 'DVHComp_All.png'); close;
clear OptModel Dose

% % Solo ConstRBE
% OptModel = cell(6,2);
% OptModel(1,1) = {' ConstRBEOpt '};
% OptModel(2,1) = {' RBEMCNOpt '};
% OptModel(3,1) = {' RBEUCMOpt '};
% OptModel([1,2,3],2) = {' ConstRBE '};
% 
% Dose{1,1} = ResultConstRBE.Optimized.resultGUI.RBExD;
% Dose{2,1} = ResultRBEMCN.ConstRBEreCalc.resultGUI.RBExD;
% Dose{3,1} = ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD;
% 
% image_ConstRBE = prueba_compDVH ( pln , cst, Dose, Regions, OptModel);
% %saveas(image_ConstRBE, 'DVHComp_ConstRBE.png'); close;
% clear OptModel Dose
% 
% % Solo RBEMCN
% OptModel = cell(6,2);
% OptModel(1,1) = {' ConstRBEOpt '};
% OptModel(2,1) = {' RBEMCNOpt '};
% OptModel(3,1) = {' RBEUCMOpt '};
% OptModel([1,2,3],2) = {' RBEMCN '};
% 
% Dose{1,1} = ResultConstRBE.RBEMCNreCalc.resultGUI.RBExD;
% Dose{2,1} = ResultRBEMCN.Optimized.resultGUI.RBExD;
% Dose{3,1}= ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD;
% 
% image_RBEMCN = prueba_compDVH ( pln , cst, Dose, Regions, OptModel);
% 
% %saveas(image_RBEMCN, 'DVHComp_RBEMCN.png')

clearvars -except ct phantomtype cst pln ResultRBEMCN ResultRBEUCM ResultConstRBE midRBE ResultPhysical


%% 15 - Calculo de las estadisticas generales de dosis

%prueba_DVHstatsComp (pln, cst, Dose, refVol, refGy, Models, Name, FigRem)

% Estadisticas de dosis para el caso de ConstRBE optimizado
Dose{1,1} = ResultConstRBE.Optimized.resultGUI.RBExD;
Dose{2,1} = ResultConstRBE.RBEMCNreCalc.resultGUI.RBExD;
Models{1,1} = 'ConstRBE'; Models{1,2} = 'Optimized';
Models{2,1} = 'RBEMCN';Models{2,2} = 'RBEMCNreCalc';

DoseStatistics.ConstRBEOpt = prueba_DVHstatsComp(pln, cst, Dose,[],[], Models, 'ConstRBEOpt',0);
clear Dose Models

% Estadisticas de dosis para el modelo de McNamara optimizado
Dose{1,1} = ResultRBEMCN.ConstRBEreCalc.resultGUI.RBExD;
Dose{2,1} = ResultRBEMCN.Optimized.resultGUI.RBExD;
Models{1,1} = 'ConstRBE'; Models{1,2} = 'ConstRBEreCalc';
Models{2,1} = 'RBEMCN';Models{2,2} = 'Optimized';

DoseStatistics.RBEMCNOpt = prueba_DVHstatsComp(pln, cst, Dose,[],[], Models, 'RBEMCNOpt',0);
clear Dose Models

% Estadisticas de dosis para el modelo UCM optimizado
Dose{1,1} = ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD;
Dose{2,1} = ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD;
%Dose{3,1} = ResultRBEUCM.Optimized.resultGUI.RBExD;
Models{1,1} = 'ConstRBE'; Models{1,2} = 'ConstRBEreCalc';
Models{2,1} = 'RBEMCN';Models{2,2} = 'RBEMCNreCalc';
%Models{3,1} = 'RBEUCM';Models{3,2} = 'Optimized';

DoseStatistics.RBEUCMOpt = prueba_DVHstatsComp(pln, cst, Dose,[],[], Models, 'RBEUCMOpt',0);
clear Dose Models

clearvars -except ct cst phantomtype pln ResultRBEMCN ResultRBEUCM ResultConstRBE ResultPhysical midRBE  DoseStatistics


%% 16 - Calculo de NTCP

% Calculos NTCP para ConstRBEOpt
NTCP.ConstRBEOpt.Optimized = prueba_NTCPcalc(pln, cst, phantomtype, ResultConstRBE.Optimized.resultGUI.RBExD);
NTCP.ConstRBEOpt.RBEMCNreCalc = prueba_NTCPcalc(pln, cst, phantomtype, ResultConstRBE.RBEMCNreCalc.resultGUI.RBExD);

% Calculos NTCP para RBEMCNOpt
NTCP.RBEMCNOpt.ConstRBEreCalc = prueba_NTCPcalc(pln, cst, phantomtype, ResultRBEMCN.ConstRBEreCalc.resultGUI.RBExD);
NTCP.RBEMCNOpt.Optimized = prueba_NTCPcalc(pln, cst, phantomtype, ResultRBEMCN.Optimized.resultGUI.RBExD);

% Calculos NTCP para RBEUCMOpt
NTCP.RBEUCMOpt.ConstRBEreCalc = prueba_NTCPcalc(pln, cst, phantomtype, ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD);
NTCP.RBEUCMOpt.RBEMCNreCalc = prueba_NTCPcalc(pln, cst, phantomtype, ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD);


%% 17 - "Resumen" de valores de NTCP
NTCP_bio_const = nan(7,1);
NTCP_bio_MCN = nan(7,1);
NTCP_bio_UCM = nan(7,1);

NTCP_bio_const(1) = NTCP.ConstRBEOpt.RBEMCNreCalc.Fukahori.NTCP(2).NTCP;
NTCP_bio_const(2) = NTCP.ConstRBEOpt.RBEMCNreCalc.Burman.Rectum.NTCP;
NTCP_bio_const(3) = NTCP.ConstRBEOpt.RBEMCNreCalc.Cheung.Rectum.NTCP(2).NTCP;
NTCP_bio_const(4) = NTCP.ConstRBEOpt.RBEMCNreCalc.Liu.NTCP;
NTCP_bio_const(5) = NTCP.ConstRBEOpt.RBEMCNreCalc.Tucker.NTCP;
NTCP_bio_const(6) = NTCP.ConstRBEOpt.RBEMCNreCalc.Peeters.Bleeding.NTCP;
NTCP_bio_const(7) = NTCP.ConstRBEOpt.RBEMCNreCalc.Schaake.NTCP(2).NTCP;

NTCP_bio_MCN(1) = NTCP.RBEMCNOpt.Optimized.Fukahori.NTCP(2).NTCP;
NTCP_bio_MCN(2) = NTCP.RBEMCNOpt.Optimized.Burman.Rectum.NTCP;
NTCP_bio_MCN(3) = NTCP.RBEMCNOpt.Optimized.Cheung.Rectum.NTCP(2).NTCP;
NTCP_bio_MCN(4) = NTCP.RBEMCNOpt.Optimized.Liu.NTCP;
NTCP_bio_MCN(5) = NTCP.RBEMCNOpt.Optimized.Tucker.NTCP;
NTCP_bio_MCN(6) = NTCP.RBEMCNOpt.Optimized.Peeters.Bleeding.NTCP;
NTCP_bio_MCN(7) = NTCP.RBEMCNOpt.Optimized.Schaake.NTCP(2).NTCP;

NTCP_bio_UCM(1) = NTCP.RBEUCMOpt.RBEMCNreCalc.Fukahori.NTCP(2).NTCP;
NTCP_bio_UCM(2) = NTCP.RBEUCMOpt.RBEMCNreCalc.Burman.Rectum.NTCP;
NTCP_bio_UCM(3) = NTCP.RBEUCMOpt.RBEMCNreCalc.Cheung.Rectum.NTCP(2).NTCP;
NTCP_bio_UCM(4) = NTCP.RBEUCMOpt.RBEMCNreCalc.Liu.NTCP;
NTCP_bio_UCM(5) = NTCP.RBEUCMOpt.RBEMCNreCalc.Tucker.NTCP;
NTCP_bio_UCM(6) = NTCP.RBEUCMOpt.RBEMCNreCalc.Peeters.Bleeding.NTCP;
NTCP_bio_UCM(7) = NTCP.RBEUCMOpt.RBEMCNreCalc.Schaake.NTCP(2).NTCP;

% Mostrar valores quitando Cheung y Fukahori
mascara = [2 4 5 6 7];
[NTCP_bio_const(mascara), NTCP_bio_MCN(mascara), NTCP_bio_UCM(mascara)]
[mean(NTCP_bio_const(mascara)), mean(NTCP_bio_MCN(mascara)), mean(NTCP_bio_UCM(mascara))]


%% 18 - Resultados de dij y resultGUI para ver en la GUI

% NOTA IMPORTANTE: Solo se puede cargar uno cada vez

% Resultados para dosis física optimizada
     dij = ResultPhysical.Optimized.dij;
    % resultGUI = ResultPhysical.Optimized.resultGUI;
    % resultGUI = ResultPhysical.ConstRBEreCalc.resultGUI;
     resultGUI = ResultPhysical.RBEMCNreCalc.resultGUI;

% Resultados para RBExD para el modelo ConstRBE optimizado
    % dij = ResultConstRBE.Optimized.dij;
    % resultGUI = ResultConstRBE.Optimized.resultGUI;
    % resultGUI = ResultConstRBE.RBEMCNreCalc.resultGUI;


% Resultados para RBExD para el modelo RBEMCN optimizado
    % dij = ResultRBEMCN.Optimized.dij;
    % resultGUI = ResultRBEMCN.ConstRBEreCalc.resultGUI;
    % resultGUI = ResultRBEMCN.Optimized.resultGUI;

% Resultados para RBExD para el modelo RBEUCM optimizado
    % dij = ResultRBEUCM.Optimized.dij;
    % resultGUI = ResultRBEUCM.ConstRBEreCalc.resultGUI;
    % resultGUI = ResultRBEUCM.RBEMCNreCalc.resultGUI;


%% Cobertura del PTV68
DoseStatistics.ConstRBEOpt.Optimized(6)
DoseStatistics.RBEUCMOpt.ConstRBEreCalc(6)
DoseStatistics.RBEMCNOpt.ConstRBEreCalc(6) % hola