%% Hacer carga del phantom

close all
clear
if ispc
    opengl software
end
load PROSTATE.mat; phantomtype = 'Prostate';
%load HEAD_AND_neck.mat; phantomtype = 'Head and Neck';
%load BOXPHANTOM.mat; phantomtype = 'Test';
%load TG119.mat; phantomtype = 'Test';


%% Carga de parametros alpha y beta

cst = prueba_abLoader (cst, phantomtype);

%% Introduccion de los datos 

pln.bixelWidth      = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.gantryAngles    = [0]; % [Â°]
pln.couchAngles     = [0]; % [Â°]
pln.numOfBeams      = numel(pln.gantryAngles);
pln.numOfVoxels     = prod(ct.cubeDim);
pln.isoCenter       = ones(pln.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.voxelDimensions = ct.cubeDim;
pln.radiationMode   = 'protons';             % either photons / protons / carbonm
pln.scenGenType     = 'nomScen';             % scenario creation type'nomScen'  'wcScen' 'impScen' 'rndScen'

pln.numOfFractions  = 30;
pln.runSequencing   = false; % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
pln.runDAO          = false; % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.machine         = 'Generic';
pln.robOpt          = false;

%% Cálculo y optimización para dosis fisica

pln.bioOptimization = 'none_physicalDose';   % none_physicalDose: physical optimization;                              constRBE_RBExD; constant RBE of 1.1;  
                                             % MCN_effect; McNamara-variable RBE model for protons (effect based)     MCN_RBExD; McNamara-variable RBE model for protons (RBExD) based
                                             % WED_effect; Wedenberg-variable RBE model for protons (effect based)    WED_RBExD; Wedenberg-variable RBE model for protons (RBExD) based
                                             % LEM_effect: effect-based optimization;                                 LEM_RBExD: optimization of RBE-weighted dose

% retrieve model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,pln.bioOptimization);
% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,pln.scenGenType); 
% generate steering file
stf = matRad_generateStf(ct,cst,pln);

% dose calculation
if strcmp(pln.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,pln,cst);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    dij = matRad_calcParticleDose(ct,stf,pln,cst);
end

resultGUI = matRad_fluenceOptimization(dij,cst,pln);

pause(1.5);
close all
clc

ResultPhysical.Optimized.dij = dij;
ResultPhysical.Optimized.resultGUI = resultGUI;


%% Recalculo de dosis para ConstRBE Y RBEMCN

%[ResultRBEUCM.ConstRBEreCalc.dij,ResultRBEUCM.ConstRBEreCalc.resultGUI] = prueba_RecalcDose(ResultRBEUCM.Optimized.resultGUI, ct, stf, pln, cst, 'constRBE_RBExD');
[~ ,ResultPhysical.ConstRBEreCalc.resultGUI] = prueba_RecalcDose(ResultPhysical.Optimized.resultGUI, ct, stf, pln, cst, 'constRBE_RBExD');

%[ResultRBEUCM.RBEMCNreCalc.dij,ResultRBEUCM.RBEMCNreCalc.resultGUI] = prueba_RecalcDose(ResultRBEUCM.Optimized.resultGUI, ct, stf, pln, cst, 'MCN_RBExD');
[~ ,ResultPhysical.RBEMCNreCalc.resultGUI] = prueba_RecalcDose(ResultPhysical.Optimized.resultGUI, ct, stf, pln, cst, 'MCN_RBExD');

clear dij resultGUI 

%% Calculo y optimizacion de dosis considerando el RBE = 1.1

pln.bioOptimization = 'constRBE_RBExD';      % none_physicalDose: physical optimization;                              constRBE_RBExD; constant RBE of 1.1;  
                                             % MCN_effect; McNamara-variable RBE model for protons (effect based)     MCN_RBExD; McNamara-variable RBE model for protons (RBExD) based
                                             % WED_effect; Wedenberg-variable RBE model for protons (effect based)    WED_RBExD; Wedenberg-variable RBE model for protons (RBExD) based
                                             % LEM_effect: effect-based optimization;                                 LEM_RBExD: optimization of RBE-weighted dose

% retrieve model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,pln.bioOptimization);
% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,pln.scenGenType); 
% generate steering file
stf = matRad_generateStf(ct,cst,pln);

% dose calculation
if strcmp(pln.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,pln,cst);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    dij = matRad_calcParticleDose(ct,stf,pln,cst);
end

resultGUI = matRad_fluenceOptimization(dij,cst,pln);

pause(1.5);
close all
clc

ResultConstRBE.Optimized.dij = dij;
ResultConstRBE.Optimized.resultGUI = resultGUI;


%% Recalculo de dosis para RBEMCN

%[ResultConstRBE.RBEMCNreCalc.dij,ResultConstRBE.RBEMCNreCalc.resultGUI] = prueba_RecalcDose(resultGUI, ct, stf, pln, cst, 'MCN_RBExD');
[~ ,ResultConstRBE.RBEMCNreCalc.resultGUI] = prueba_RecalcDose(resultGUI, ct, stf, pln, cst, 'MCN_RBExD');

clear dij resultGUI 

%% Cambio a modelo McNamara, calculo de dosis y optimizacion de este


pln.bioOptimization = 'MCN_RBExD';           % none_physicalDose: physical optimization;                              constRBE_RBExD; constant RBE of 1.1;  
                                             % MCN_effect; McNamara-variable RBE model for protons (effect based)     MCN_RBExD; McNamara-variable RBE model for protons (RBExD) based
                                             % WED_effect; Wedenberg-variable RBE model for protons (effect based)    WED_RBExD; Wedenberg-variable RBE model for protons (RBExD) based
                                             % LEM_effect: effect-based optimization;      
                                             
% retrieve model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,pln.bioOptimization);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,pln.scenGenType); 

stf = matRad_generateStf(ct,cst,pln);

if strcmp(pln.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,pln,cst);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    dij = matRad_calcParticleDose(ct,stf,pln,cst);
end

resultGUI = matRad_fluenceOptimization(dij,cst,pln);
pause(1.5);
close all
clc

ResultRBEMCN.Optimized.dij = dij;
ResultRBEMCN.Optimized.resultGUI = resultGUI;


%% Recalculo de dosis para ConstRBE

%[ResultRBEMCN.ConstRBEreCalc.dij,ResultRBEMCN.ConstRBEreCalc.resultGUI] = prueba_RecalcDose(ResultRBEMCN.Optimized.resultGUI, ct, stf, pln, cst,'constRBE_RBExD');
[~ ,ResultRBEMCN.ConstRBEreCalc.resultGUI] = prueba_RecalcDose(ResultRBEMCN.Optimized.resultGUI, ct, stf, pln, cst,'constRBE_RBExD');

clear dij resultGUI 

%% Cambio a modelo UCM, calculo de dosis y optimizacion de este

pln.bioOptimization = 'UCM_RBExD';           % none_physicalDose: physical optimization;                              constRBE_RBExD; constant RBE of 1.1;  
                                             % MCN_effect; McNamara-variable RBE model for protons (effect based)     MCN_RBExD; McNamara-variable RBE model for protons (RBExD) based
                                             % WED_effect; Wedenberg-variable RBE model for protons (effect based)    WED_RBExD; Wedenberg-variable RBE model for protons (RBExD) based
                                             % LEM_effect: effect-based optimization;      
                                             
% retrieve model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,pln.bioOptimization);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,pln.scenGenType); 

stf = matRad_generateStf(ct,cst,pln);

if strcmp(pln.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,pln,cst);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'carbon')
    dij = matRad_calcParticleDose(ct,stf,pln,cst);
end

resultGUI = matRad_fluenceOptimization(dij,cst,pln);
pause(1.5);
close all
clc

ResultRBEUCM.Optimized.dij = dij;
ResultRBEUCM.Optimized.resultGUI = resultGUI;


%% Recalculo de dosis para ConstRBE Y RBEMCN

%[ResultRBEUCM.ConstRBEreCalc.dij,ResultRBEUCM.ConstRBEreCalc.resultGUI] = prueba_RecalcDose(ResultRBEUCM.Optimized.resultGUI, ct, stf, pln, cst, 'constRBE_RBExD');
[~ ,ResultRBEUCM.ConstRBEreCalc.resultGUI] = prueba_RecalcDose(ResultRBEUCM.Optimized.resultGUI, ct, stf, pln, cst, 'constRBE_RBExD');

%[ResultRBEUCM.RBEMCNreCalc.dij,ResultRBEUCM.RBEMCNreCalc.resultGUI] = prueba_RecalcDose(ResultRBEUCM.Optimized.resultGUI, ct, stf, pln, cst, 'MCN_RBExD');
[~ ,ResultRBEUCM.RBEMCNreCalc.resultGUI] = prueba_RecalcDose(ResultRBEUCM.Optimized.resultGUI, ct, stf, pln, cst, 'MCN_RBExD');

clear dij resultGUI 

%% Graficas de perfil de dosis

%ProfileType = longitudinal // lateral
%DisplayOption = physicalDose // RBExD // physical_vs_RBExD
%Para hacer comparaciones entre modelos DisplayOption == RBExD // physical_vs_RBExD
% prueba_DoseGraphs (ct, pln, cst, NumBeam, ProfileType, DisplayOption, Result, Model1, Result2, Model2)

prueba_DoseGraphs (ct, pln, cst,1, 'longitudinal', 'physical_vs_RBExD', ResultPhysical.ConstRBEreCalc.resultGUI, 'ConstRBE',ResultPhysical.RBEMCNreCalc.resultGUI,'RBEMCN')

prueba_DoseGraphs (ct, pln, cst,1, 'longitudinal', 'physical_vs_RBExD', ResultConstRBE.Optimized.resultGUI, 'ConstRBE', ResultConstRBE.RBEMCNreCalc, 'RBEMCN')

prueba_DoseGraphs (ct, pln, cst,1, 'longitudinal', 'physical_vs_RBExD', ResultRBEMCN.Optimized.resultGUI, 'RBEMCN',ResultRBEMCN.ConstRBEreCalc.resultGUI,'ConstRBE')

prueba_DoseGraphs (ct, pln, cst,1, 'longitudinal', 'physical_vs_RBExD', ResultRBEUCM.RBEMCNreCalc.resultGUI, 'RBEMCN',ResultRBEUCM.ConstRBEreCalc.resultGUI,'ConstRBE')


%% Graficas 2D de dosis
% prueba_DoseIntens (ct, pln, Dose, z_cut, TypeDose, Model)

% physicalDose Optimized
prueba_DoseIntens (ct, pln, ResultPhysical.Optimized.resultGUI.physicalDose, [], 'Physical', 'Physical');
prueba_DoseIntens (ct, pln, ResultPhysical.ConstRBEreCalc.resultGUI.RBExD, [], 'RBExD', 'ConstRBE');
prueba_DoseIntens (ct, pln, ResultPhysical.RBEMCNreCalc.resultGUI.RBExD, [], 'RBExD', 'MCN');

% ConstRBE Optimized
prueba_DoseIntens (ct, pln, ResultConstRBE.Optimized.resultGUI.physicalDose, [], 'Physical', 'Physical');
prueba_DoseIntens (ct, pln, ResultConstRBE.Optimized.resultGUI.RBExD, [], 'RBExD', 'ConstRBE');
prueba_DoseIntens (ct, pln, ResultConstRBE.RBEMCN.resultGUI.RBExD, [], 'RBExD', 'RBEMCN');

% MCN Optimized
prueba_DoseIntens (ct, pln, ResultRBEMCN.Optimized.resultGUI.physicalDose, [], 'Physical', 'Physical');
prueba_DoseIntens (ct, pln, ResultRBEMCN.ConstRBEreCalc.resultGUI.RBExD, [], 'RBExD', 'ConstRBE');
prueba_DoseIntens (ct, pln, ResultRBEMCN.Optimized.resultGUI.RBExD, [], 'RBExD', 'RBEMCN');

% UCM Optimized
prueba_DoseIntens (ct, pln, ResultRBEUCM.Optimized.resultGUI.physicalDose, [], 'Physical', 'Physical');
prueba_DoseIntens (ct, pln, ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD, [], 'RBExD', 'ConstRBE');
prueba_DoseIntens (ct, pln, ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD, [], 'RBExD', 'RBEMCN');


%% Representacion DVH todas las VOI

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
OptModel(1,1) = {'UCMMCN Optimized'};
OptModel(1,2) = {'ConstRBE'};
OptModel(2,2) = {'RBEMCN'};
OptModel(3,2) = {'UCMMCN'};

Dose{1,1} = ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD;
Dose{2,1} = ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD;
Dose{3,1} = ResultRBEUCM.Optimized.resultGUI.RBExD;

prueba_compDVH ( pln , cst, Dose, [], OptModel);
clear OptModel


%Comparacion entre los tres modelos optimizados
OptModel = cell(3,2);
OptModel(1,1) = {'3 Models Optimized'};
OptModel(1,2) = {'ConstRBE'};
OptModel(2,2) = {'RBEMCN'};
OptModel(3,2) = {'UCMMCN'};

Dose{1,1} = ResultConstRBE.Optimized.resultGUI.RBExD;
Dose{2,1} = ResultRBEMCN.Optimized.resultGUI.RBExD;
Dose{3,1} = ResultRBEUCM.Optimized.resultGUI.RBExD;

prueba_compDVH ( pln , cst, Dose, [], OptModel);


clearvars -except ct phantomtype cst pln ResultRBEMCN ResultRBEUCM ResultConstRBE 

%% Comparacion DVH para VOIs específicas 

if strcmp(phantomtype, 'Prostate') >0
    Regions{1,1} = 'Rectum';
    Regions{2,1} = 'PTV 68';
    
elseif strcmp(phantomtype, 'Head and Neck') >0
    Regions = [];
end 

OptModel = cell(6,2);
OptModel([1,2],1) = {' ConstRBE Model '};
OptModel([3,4],1) = {' RBEMCN Model '};
OptModel([5,6],1) = {' RBEUCM Plan '};
OptModel([1,3,5],2) = {' ConstRBE '};
OptModel([2,4,6],2) = {' RBEMCN '};

Dose{1,1} = ResultConstRBE.Optimized.resultGUI.RBExD;
Dose{2,1} = ResultConstRBE.RBEMCNreCalc.resultGUI.RBExD;
Dose{3,1} = ResultRBEMCN.Optimized.resultGUI.RBExD;
Dose{4,1} = ResultRBEMCN.ConstRBEreCalc.resultGUI.RBExD;
Dose{5,1} = ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD;
Dose{6,1}= ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD;

prueba_compDVH ( pln , cst, Dose, Regions, OptModel);


clearvars -except ct phantomtype cst pln ResultRBEMCN ResultRBEUCM ResultConstRBE 


%% Calculo de las estadisticas de dosis

% Estadisticas de dosis para el caso de ConstRBE optimizado
DoseStatistics.ConstRBEOpt = prueba_DVHstatsComp (ResultConstRBE.Optimized.resultGUI.RBExD,...
    ResultConstRBE.RBEMCNreCalc.resultGUI.RBExD, [],'ConstRBE',[],[],pln, cst,1,1);

% Estadisticas de dosis para el modelo de McNamara optimizado
DoseStatistics.RBEMCNOpt = prueba_DVHstatsComp (ResultRBEMCN.ConstRBEreCalc.resultGUI.RBExD,...
    ResultRBEMCN.Optimized.resultGUI.RBExD,[],'McNamara''s model',[],[],pln, cst,1,1);

% Estadisticas de dosis para el modelo UCM optimizado
DoseStatistics.RBEUCMOpt.NTCP_DVHStats = prueba_DVHstatsComp (ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD,...
    ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD,...
    ResultRBEUCM.Optimized.resultGUI.RBExD,...
    'UCM''s model',[],[],pln, cst,1,1);

% Estadisticas para los 3 modelos optimizados
DoseStatistics.AllOpt = prueba_DVHstatsComp (ResultConstRBE.Optimized.resultGUI.RBExD,...
    ResultRBEMCN.Optimized.resultGUI.RBExD,...
    ResultRBEUCM.Optimized.resultGUI.RBExD,...
    '3 models',[],[],pln, cst,1,1);

clearvars -except ct cst phantomtype pln ResultRBEMCN ResultRBEUCM ResultConstRBE DoseStatistics

%% Valores V_x y D_x necesarios para los calculos de los modelos NTCP
if strcmp(phantomtype, 'Prostate') > 0
    refGy = [];
    refVol = [70 95 98];
    EasyStats = 0;
    
    
elseif strcmp(phantomtype, 'Head and Neck') >0
    refGy = [];
    refVol = [95 98];
    EasyStats = 1;
end

% Estadisticas de dosis para el caso de ConstRBE optimizado
NTCP_DVHStats.ConstRBEOpt = prueba_DVHstatsComp (ResultConstRBE.Optimized.resultGUI.RBExD,...
    ResultConstRBE.RBEMCNreCalc.resultGUI.RBExD, [],'ConstRBE', refGy, refVol, pln, cst,0, EasyStats);

% Estadisticas de dosis para el modelo de McNamara optimizado
NTCP_DVHStats.RBEMCNOpt = prueba_DVHstatsComp (ResultRBEMCN.ConstRBEreCalc.resultGUI.RBExD,...
    ResultRBEMCN.Optimized.resultGUI.RBExD, [], 'McNamara''s model', refGy, refVol, pln, cst,0, EasyStats);

% Estadisticas de dosis para el modelo UCM optimizado
NTCP_DVHStats.RBEUCMOpt = prueba_DVHstatsComp (ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD,...
    ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD,...
    ResultRBEUCM.Optimized.resultGUI.RBExD,...
    'UCM''s model', refGy, refVol, pln, cst,0, EasyStats);

% Estadisticas para los 3 modelos optimizados
NTCP_DVHStats.AllOpt = prueba_DVHstatsComp (ResultConstRBE.Optimized.resultGUI.RBExD,...
    ResultRBEMCN.Optimized.resultGUI.RBExD,...
    ResultRBEUCM.Optimized.resultGUI.RBExD,...
    '3 models', refGy, refVol, pln, cst,0, EasyStats);

clearvars -except ct phantomtype cst pln ResultRBEMCN ResultRBEUCM ResultConstRBE DoseStatistics NTCP_DVHStats

%% Calculo y comparacion de ambos NTCP

NTCP = cell(3,2);
NTCP{1,1} = 'ConstRBEOpt';
NTCP{2,1} = 'MCNOpt';
NTCP{3,1} = 'UCMOpt model'; 
NTCP{4,1} = '3 Models optimized';


NTCP{1,2} = prueba_NTCPcalc (cst, phantomtype, NTCP_DVHStats.ConstRBEOpt, NTCP_DVHStats.ConstRBEOpt, NTCP_DVHStats.ConstRBEOpt, 0);

NTCP{2,2} = prueba_NTCPcalc (cst, phantomtype, NTCP_DVHStats.RBEMCNOpt, NTCP_DVHStats.RBEMCNOpt, NTCP_DVHStats.RBEMCNOpt, 0);

NTCP{3,2} = prueba_NTCPcalc (cst, phantomtype, NTCP_DVHStats.RBEUCMOpt, NTCP_DVHStats.RBEUCMOpt, NTCP_DVHStats.RBEUCMOpt, 1);

NTCP{4,2} = prueba_NTCPcalc (cst, phantomtype, NTCP_DVHStats.ConstRBEOpt, NTCP_DVHStats.RBEMCNOpt, NTCP_DVHStats.RBEUCMOpt, 1);


clearvars -except ct phantomtype cst pln ResultRBEMCN ResultRBEUCM ResultConstRBE DoseStatistics NTCP_DVHStats NTCP