%% Estructuracion del programa

% 1 - Carga del phantom/paciente
% 2 - Carga de parametros alpha y beta
% 3 - Definicion de restricciones "distintas"
% 4 - Introduccion de los datos basicos
% 5 - Seleccion de calculos
% 6 - Calculos
% 7 - Exportacion de resultados a la GUI
% 8 - Apertura de la GUI

%% 1 - Carga del phantom/paciente

close all
clear
if ~ismac
    opengl software
end

load LIVER_test.mat; phantomtype = 'Liver';

%% 2 - Carga de parametros alpha y beta

cst = prueba_abLoader (cst, phantomtype);

%% 3 - Definicion de restricciones "distintas"

%       .type                   || NaN parameters               ||OAR or Target         ||Total or dose/frac 
 
% 'square underdosing'          ||.EUD - .volume                ||OAR/Target            ||Total
% 'square overdosing'           ||.EUD - .volume                ||OAR/Target            ||Total                
% 'square deviation'            ||.EUD - .volume                ||Target                ||Total
% 'mean'                        ||.dose - .EUD - .volume        ||OAR                   ||Total
% 'EUD'                         ||.dose - .volume               ||OAR/Target            ||Total        
% 'min dose constraint'         ||.penalty - .EUD - .volume     ||Target                ||
% 'max dose constraint'         ||.penalty - .EUD - .volume     ||OAR/Target            ||
% 'min mean dose constraint'    ||.EUD - .volume                ||Target                ||
% 'max mean dose constraint'    ||.EUD - .volume                ||OAR/Target            ||
% 'min EUD constraint'          ||.penalty - .dose - .volume    ||OAR/Target            ||Total
% 'max EUD constraint'          ||.penalty - .dose - .volume    ||OAR/Target            ||Total
% 'min DVH constraint'          ||.penalty - .EUD               ||Target                ||
% 'max DVH constraint'          ||.penalty - .EUD               ||OAR/Target            ||
% 'min DVH objective'           ||.EUD                          ||Target                ||Total
% 'max DVH objective'           ||.EUD                          ||OAR/Target            ||Total

% VOLUMEN EN TANTO POR CIENTO SIEMPRE!

%% Objetivos

% PTV
cst{15,5}.Priority = 1;
cst{15,6}.type = 'square deviation';
cst{15,6}.dose = 48;
cst{15,6}.penalty = 1000;
cst{15,6}.EUD = NaN;
cst{15,6}.volume = NaN;
cst{15,6}.robustness = 'none';

% Skin
cst{14,5}.Priority = 2;
cst{14,6}.type = 'square overdosing';
cst{14,6}.dose = 21;
cst{14,6}.penalty = 300;
cst{14,6}.EUD = NaN;
cst{14,6}.volume = NaN;
cst{14,6}.robustness = 'none';

% Stomach
cst{4,5}.Priority = 2;
cst{4,6}.type = 'square overdosing';
cst{4,6}.dose = 64;
cst{4,6}.penalty = 300;
cst{4,6}.EUD = NaN;
cst{4,6}.volume = NaN;
cst{4,6}.robustness = 'none';

% Small Bowel
cst{5,5}.Priority = 2;
cst{5,6}.type = 'square overdosing';
cst{5,6}.dose = 64;
cst{5,6}.penalty = 300;
cst{5,6}.EUD = NaN;
cst{5,6}.volume = NaN;
cst{5,6}.robustness = 'none';

% Large Bowel
cst{6,5}.Priority = 2;
cst{6,6}.type = 'square overdosing';
cst{6,6}.dose = 64;
cst{6,6}.penalty = 300;
cst{6,6}.EUD = NaN;
cst{6,6}.volume = NaN;
cst{6,6}.robustness = 'none';


% ----------------------------------------

%% 4 - Introduccion de los datos basicos

% meta information for treatment plan 
pln.numOfFractions = 8;
pln.radiationMode = 'protons';           % either photons / protons / helium / carbon
pln.machine = 'Generic';

% beam geometry settings
pln.propStf.bixelWidth     = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles   = [0]; % [?];
pln.propStf.couchAngles    = [0]; % [?];
pln.propStf.numOfBeams     = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter      = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln.propOpt.runDAO         = false;   % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing  = false;   % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
%pln.calcLET = true;

%% 5 - Seleccion de calculos

% Seleccion de graficas y estadisticas que mostrar (0 = Desactivado // 1 = Activado)
perfGraphs = 0;       % Graficas de perfil de dosis
perfRBEGraphs = 0;    % Graficas de perfil de dosis vs RBE
DGraphs = 0;          % Graficas de dosis 2D en z = z(dij max)
DVHGraphs = 0;        % Representacion de DVH (1 = Generales // 2 = Especificas)
DVHStats = 0;         % Calculo de las estadisticas generales de dosis (1 = C�lculo de estad�sticas // 2 = Renormalizacion de dosis )

GraphSel = [perfGraphs perfRBEGraphs DGraphs DVHGraphs DVHStats];


% Seleccion de estadisticas V_X y D_X
refVol = []; % Valores D_X
refGy = []; % Valores V_X

StatsRef{1,1} = refVol;
StatsRef{2,1} = refGy;

% Seleccion de modelos de dosis a calcular 

ConstRBE{1,1} = 0;  ConstRBE{1,2} = 0;
RBEMCN{1,1} = 0;    RBEMCN{1,2} = 0;
RBEUCM{1,1} = 0;    RBEUCM{1,2} = 0;

DoseRecalc{1,1} = ConstRBE;
DoseRecalc{2,1} = RBEMCN;
DoseRecalc{3,1} = RBEUCM;
 

DVHRegions{1,1} = 'PTV';
DVHRegions{2,1} = 'Skin';
VOIType{1,1}    = 'TARGET';
VOIType{2,1}    = 'OAR';


%% 6 - Calculos

graphlaunch = 0;

while graphlaunch < 1
    
    CompDVH{1,1} = DVHRegions;
    CompDVH{2,1} = VOIType;
    
    if exist('ResultConstRBE','var') > 0 && exist('ResultRBEMCN', 'var') > 0 && exist('ResultRBEUCM', 'var') > 0
        % Si ya se ha realizado un calculo de todas las matrices de dosis y solo se quiere reevaluar alguna de ellas
        clear DoseResults
        DoseResults{1,1} = ResultConstRBE;
        DoseResults{1,2} = ResultRBEMCN;
        DoseResults{1,3} = ResultRBEUCM;
        
        [~, ResultConstRBE, ResultRBEMCN, ResultRBEUCM, DoseStatistics, NTCP, meanNTCP, Renorm] = ...
            prueba_NTCP(cst, pln, ct, phantomtype, DoseStatistics, GraphSel, DoseRecalc, DoseResults, StatsRef, CompDVH);
        
        graphlaunch = 1;
    else
        % Si no se ha calculado ninguna vez los resultados, ignora DoseRecalc y calcula todas las matrices de dosis automaticamente
        clear DoseResults
        GraphSel_ignore = [0 0 0 0 0];
        DoseStatistics = 'Not evaluated';
        
        ConstRBE_ig{1,1} = 1;  ConstRBE_ig{1,2} = 0;
        RBEMCN_ig{1,1} = 1;    RBEMCN_ig{1,2} = 0;
        RBEUCM_ig{1,1} = 1;    RBEUCM_ig{1,2} = 0;
        
        DoseRecalc_ig{1,1} = ConstRBE_ig;
        DoseRecalc_ig{2,1} = RBEMCN_ig;
        DoseRecalc_ig{3,1} = RBEUCM_ig;
        
        DoseResults = [];
        [~, ResultConstRBE, ResultRBEMCN, ResultRBEUCM, DoseStatistics, NTCP, meanNTCP, Renorm] = ...
            prueba_NTCP(cst, pln, ct, phantomtype, DoseStatistics, GraphSel_ignore, DoseRecalc_ig, DoseResults, StatsRef, []);
        
    end
end

    clearvars -except ct cst CompDVH GraphSel phantomtype pln stf ResultRBEMCN ResultRBEUCM ResultConstRBE ResultPhysical midRBE  DoseStatistics NTCP GraphSel DVHRegions

    
%% 7 - Exportacion de resultados a la GUI

% NOTA IMPORTANTE: Solo se puede cargar uno cada vez
% NOTA 2: Segun el resultado que se quiera exportar hay que modificar el parametro modelname 

clear dij resultGUI quantityOpt modelName stf

quantityOpt         = 'RBExD';
modelName           = 'MCN';       % 'constRBE', 'MCN', 'UCM'
scenGenType = 'nomScen';
pln.bioParam = matRad_bioModel(pln.radiationMode, quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,scenGenType);
stf = matRad_generateStf(ct,cst,pln);

% Resultados para RBExD para el modelo ConstRBE optimizado
% dij = ResultConstRBE.Optimized.dij;
% resultGUI = ResultConstRBE.Optimized.resultGUI;
% resultGUI = ResultConstRBE.RBEMCNreCalc.resultGUI;


% Resultados para RBExD para el modelo RBEMCN optimizado
dij = ResultRBEMCN.Optimized.dij;
%resultGUI = ResultRBEMCN.ConstRBEreCalc.resultGUI;
resultGUI = ResultRBEMCN.Optimized.resultGUI;


% Resultados para RBExD para el modelo RBEUCM optimizado
% dij = ResultRBEUCM.Optimized.dij;
% resultGUI = ResultRBEUCM.ConstRBEreCalc.resultGUI;
% resultGUI = ResultRBEUCM.RBEMCNreCalc.resultGUI;

%% 8 - Apertura de la GUI

matRadGUI
