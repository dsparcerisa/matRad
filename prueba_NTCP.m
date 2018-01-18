%% Hacer carga del phantom

close all
clear
if ispc
    opengl software
end
%load PROSTATE.mat; phantomtype = 'Prostate';
%load HEAD_AND_neck.mat; phantomtype = 'Head and Neck';
%load BOXPHANTOM.mat; phantomtype = 'Test';
load TG119.mat; phantomtype = 'Test';


%% Carga de parámetros alpha y beta

cst = prueba_abLoader (cst, phantomtype);

%% Introducción de los datos 

pln.bixelWidth      = 5; % [mm] / also correspondsto lateral spot spacing for particles
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



%% Cálculo y optimización de dosis considerando el RBE = 1.1

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


%% Recálculo de dosis para RBEMCN

%[ResultConstRBE.RBEMCNreCalc.dij,ResultConstRBE.RBEMCNreCalc.resultGUI] = prueba_RecalcDose(resultGUI, ct, stf, pln, cst, 'MCN_RBExD');
[~ ,ResultConstRBE.RBEMCNreCalc.resultGUI] = prueba_RecalcDose(resultGUI, ct, stf, pln, cst, 'MCN_RBExD');

clear dij resultGUI 

%% Cambio a modelo McNamara, cálculo de dosis y optimización de este


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


%% Recálculo de dosis para ConstRBE

%[ResultRBEMCN.ConstRBEreCalc.dij,ResultRBEMCN.ConstRBEreCalc.resultGUI] = prueba_RecalcDose(ResultRBEMCN.Optimized.resultGUI, ct, stf, pln, cst,'constRBE_RBExD');
[~ ,ResultRBEMCN.ConstRBEreCalc.resultGUI] = prueba_RecalcDose(ResultRBEMCN.Optimized.resultGUI, ct, stf, pln, cst,'constRBE_RBExD');

clear dij resultGUI 

%% Cambio a modelo UCM, cálculo de dosis y optimización de este

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


%% Recálculo de dosis para ConstRBE Y RBEMCN

%[ResultRBEUCM.ConstRBEreCalc.dij,ResultRBEUCM.ConstRBEreCalc.resultGUI] = prueba_RecalcDose(ResultRBEUCM.Optimized.resultGUI, ct, stf, pln, cst, 'constRBE_RBExD');
[~ ,ResultRBEUCM.ConstRBEreCalc.resultGUI] = prueba_RecalcDose(ResultRBEUCM.Optimized.resultGUI, ct, stf, pln, cst, 'constRBE_RBExD');

%[ResultRBEUCM.RBEMCNreCalc.dij,ResultRBEUCM.RBEMCNreCalc.resultGUI] = prueba_RecalcDose(ResultRBEUCM.Optimized.resultGUI, ct, stf, pln, cst, 'MCN_RBExD');
[~ ,ResultRBEUCM.RBEMCNreCalc.resultGUI] = prueba_RecalcDose(ResultRBEUCM.Optimized.resultGUI, ct, stf, pln, cst, 'MCN_RBExD');

clear dij resultGUI 

%% Gráficas de perfil e intensidad de dosis


%ProfileType = longitudinal // lateral
%DisplayOption = physicalDose // RBExD // physical_vs_RBExD
%Para hacer comparaciones entre modelos DisplayOption == RBExD // physical_vs_RBExD
% prueba_DoseGraphs (ct, pln, cst, NumBeam, ProfileType, DisplayOption, Result, Model1, Result2, Model2)

prueba_DoseGraphs (ct, pln, cst,1, 'longitudinal', 'physical_vs_RBExD', ResultConstRBE.Optimized.resultGUI, 'ConstRBE',[],[])
%%
% prueba_DoseIntens (ct, pln, Dose, z_cut, TypeDose, Model)
prueba_DoseIntens (ct, pln, ResultConstRBE.Optimized.resultGUI.RBExD, [], 'RBExD', 'ConstRBE')


%% Representación de las comparaciones de los DVH

%Comparación para RBE constante optimizado
prueba_compDVH (ResultConstRBE.RBEMCNreCalc.resultGUI.RBExD, [],...
    ResultConstRBE.Optimized.resultGUI.RBExD, pln, cst,'Constant RBE','McNamara''s','','Constant RBE');

%Comparación para el modelo de McNamara optimizado
prueba_compDVH (ResultRBEMCN.ConstRBEreCalc.resultGUI.RBExD, [],...
    ResultRBEMCN.Optimized.resultGUI.RBExD, pln, cst,'McNamara''s model','Constant RBE','','McNamara''s');

%Comparación para el modelo UCM optimizado
prueba_compDVH (ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD,...
    ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD,...
    ResultRBEUCM.Optimized.resultGUI.RBExD, pln, cst,'UCM''s model','Constant RBE','McNamara''s','UCM''s');

%Comparación entre los tres modelos optimizados
prueba_compDVH (ResultConstRBE.Optimized.resultGUI.RBExD,...
    ResultRBEMCN.Optimized.resultGUI.RBExD,...
    ResultRBEUCM.Optimized.resultGUI.RBExD, pln, cst,'3 models','Constant RBE','McNamara''s','UCM''s');

clearvars -except ct phantomtype cst pln ResultRBEMCN ResultRBEUCM ResultConstRBE 

%% Cálculo de las estadísticas de dosis

% Estadisticas de dosis para el caso de ConstRBE optimizado
DoseStatistics.ConstRBEOpt = prueba_DVHstatsComp (ResultConstRBE.Optimized.resultGUI.RBExD,...
    ResultConstRBE.RBEMCNreCalc.resultGUI.RBExD, [],'ConstRBE',[],[],pln, cst,1,1);

% Estadísticas de dosis para el modelo de McNamara optimizado
DoseStatistics.RBEMCNOpt = prueba_DVHstatsComp (ResultRBEMCN.ConstRBEreCalc.resultGUI.RBExD,...
    ResultRBEMCN.Optimized.resultGUI.RBExD,[],'McNamara''s model',[],[],pln, cst,1,1);

% Estadísticas de dosis para el modelo UCM optimizado
DoseStatistics.RBEUCMOpt.NTCP_DVHStats = prueba_DVHstatsComp (ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD,...
    ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD,...
    ResultRBEUCM.Optimized.resultGUI.RBExD,...
    'UCM''s model',[],[],pln, cst,1,1);

% Estadísticas para los 3 modelos optimizados
DoseStatistics.AllOpt = prueba_DVHstatsComp (ResultConstRBE.Optimized.resultGUI.RBExD,...
    ResultRBEMCN.Optimized.resultGUI.RBExD,...
    ResultRBEUCM.Optimized.resultGUI.RBExD,...
    '3 models',[],[],pln, cst,1,1);

clearvars -except ct cst phantomtype pln ResultRBEMCN ResultRBEUCM ResultConstRBE DoseStatistics

%% Valores V_x y D_x necesarios para los cálculos de los modelos NTCP
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

% Estadísticas de dosis para el modelo de McNamara optimizado
NTCP_DVHStats.RBEMCNOpt = prueba_DVHstatsComp (ResultRBEMCN.ConstRBEreCalc.resultGUI.RBExD,...
    ResultRBEMCN.Optimized.resultGUI.RBExD, [], 'McNamara''s model', refGy, refVol, pln, cst,0, EasyStats);

% Estadísticas de dosis para el modelo UCM optimizado
NTCP_DVHStats.RBEUCMOpt = prueba_DVHstatsComp (ResultRBEUCM.ConstRBEreCalc.resultGUI.RBExD,...
    ResultRBEUCM.RBEMCNreCalc.resultGUI.RBExD,...
    ResultRBEUCM.Optimized.resultGUI.RBExD,...
    'UCM''s model', refGy, refVol, pln, cst,0, EasyStats);

% Estadísticas para los 3 modelos optimizados
NTCP_DVHStats.AllOpt = prueba_DVHstatsComp (ResultConstRBE.Optimized.resultGUI.RBExD,...
    ResultRBEMCN.Optimized.resultGUI.RBExD,...
    ResultRBEUCM.Optimized.resultGUI.RBExD,...
    '3 models', refGy, refVol, pln, cst,0, EasyStats);

clearvars -except ct phantomtype cst pln ResultRBEMCN ResultRBEUCM ResultConstRBE DoseStatistics NTCP_DVHStats

%% Cálculo y comparacion de ambos NTCP

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