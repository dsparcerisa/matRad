%% Estructuracion del programa

% 1 - Carga del phantom/paciente
% 2 - Carga de parametros alpha y beta
% 3 - Definicion de restricciones "distintas"
% 4 - Introduccion de los datos basicos
% 5 - Calculos
% 6 - Exportacion de resultados a la GUI
% 7 - Apertura de la GUI

%% 1 - Carga del phantom/paciente

close all
clear
if ~ismac
    opengl software
end
load TG119.mat; phantomtype = 'TG119';

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


% Core
cst{1,6}.type = 'square overdosing';
cst{1,6}.dose = 25;
cst{1,6}.EUD = NaN;
cst{1,6}.penalty = 300;
cst{1,6}.volume = NaN;
cst{1,6}.robustness = 'none';

% Target
cst{2,6}.type = 'square deviation';
cst{2,6}.dose = 50;
cst{2,6}.EUD = NaN;
cst{2,6}.penalty = 1000;
cst{2,6}.volume = NaN;
cst{2,6}.robustness = 'none';

% Body
cst{3,6}.type = 'square overdosing';
cst{3,6}.dose = 30;
cst{3,6}.EUD = NaN;
cst{3,6}.penalty = 100;
cst{3,6}.volume = NaN;
cst{3,6}.robustness = 'none';

    
%% 4 - Introduccion de los datos basicos

% meta information for treatment plan 
pln.numOfFractions = 30;
pln.radiationMode = 'protons';           % either photons / protons / helium / carbon
pln.machine = 'Generic';

% beam geometry settings
pln.propStf.bixelWidth     = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles   = [45 315]; % [?];
pln.propStf.couchAngles    = [0 0]; % [?];
pln.propStf.numOfBeams     = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter      = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln.propOpt.runDAO         = false;   % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing  = false;   % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
%pln.calcLET = true;

%% 5 - Calculos

% Seleccion de graficas y estadisticas que mostrar (0 = Desactivado // 1 = Activado)
perfGraphs = 0;       % Graficas de perfil de dosis
perfRBEGraphs = 0;    % Graficas de perfil de dosis vs RBE
DGraphs = 0;          % Graficas de dosis 2D en z = z(dij max)
DVHGraphs = 0;        % Representacion de DVH (1 = Generales // 2 = Especificas)
DVHStats = 0;         % Calculo de las estadisticas generales de dosis

GraphSel = [perfGraphs perfRBEGraphs DGraphs DVHGraphs DVHStats];

% Seleccion de modelos de dosis a calcular 
% Dosis{i,j}     j = 1 -> Activación del calculo (0 = Desactivado // 1 = Activado)
%                j = 2 -> Recalculo de dij y reoptimizacion (0 = Todo // 1 = Solo reoptimizar)

ConstRBE{1,1} = 0;  ConstRBE{1,2} = 0;
RBEMCN{1,1} = 0;    RBEMCN{1,2} = 0;
RBEUCM{1,1} = 0;    RBEUCM{1,2} = 0;

DoseRecalc{1,1} = ConstRBE;
DoseRecalc{2,1} = RBEMCN;
DoseRecalc{3,1} = RBEUCM;


% Calculos
if exist('ResultConstRBE','var') > 0 && exist('ResultRBEMCN', 'var') > 0 && exist('ResultRBEUCM', 'var') > 0
    % Si ya se ha realizado un calculo de todas las matrices de dosis y solo se quiere reevaluar alguna de ellas
    clear DoseResults
    DoseResults{1,1} = ResultConstRBE;
    DoseResults{1,2} = ResultRBEMCN;
    DoseResults{1,3} = ResultRBEUCM;
    
    [~, ResultConstRBE, ResultRBEMCN, ResultRBEUCM, DoseStatistics, NTCP, meanNTCP] = ...
        prueba_NTCP(cst, pln, ct, phantomtype, DoseStatistics, GraphSel, DoseRecalc, DoseResults);
else
    % Si no se ha calculado ninguna vez los resultados, ignora DoseRecalc y calcula todas las matrices de dosis automaticamente
    clear DoseResults
    DoseStatistics = 'Not evaluated';
    ConstRBE{1,1} = 1;  ConstRBE{1,2} = 0;
    RBEMCN{1,1} = 1;    RBEMCN{1,2} = 0;
    RBEUCM{1,1} = 1;    RBEUCM{1,2} = 0;
    
    DoseRecalc{1,1} = ConstRBE;
    DoseRecalc{2,1} = RBEMCN;
    DoseRecalc{3,1} = RBEUCM;
    DoseResults = [];
    [~, ResultConstRBE, ResultRBEMCN, ResultRBEUCM, DoseStatistics, NTCP, meanNTCP] = ...
        prueba_NTCP(cst, pln, ct, phantomtype, DoseStatistics, GraphSel, DoseRecalc, DoseResults);
end
    
    clearvars -except ct cst phantomtype pln stf ResultRBEMCN ResultRBEUCM ResultConstRBE ResultPhysical midRBE  DoseStatistics NTCP

    
%% 6 - Exportacion de resultados a la GUI

% NOTA IMPORTANTE: Solo se puede cargar uno cada vez
% NOTA 2: Segun el resultado que se quiera exportar hay que modificar el parametro modelname 

clear dij resultGUI quantityOpt modelName stf

quantityOpt         = 'RBExD';
modelName           = 'constRBE';       % 'constRBE', 'MCN', 'UCM'
scenGenType = 'nomScen';
pln.bioParam = matRad_bioModel(pln.radiationMode, quantityOpt, modelName);
pln.multScen = matRad_multScen(ct,scenGenType);
stf = matRad_generateStf(ct,cst,pln);

% Resultados para RBExD para el modelo ConstRBE optimizado
dij = ResultConstRBE.Optimized.dij;
resultGUI = ResultConstRBE.Optimized.resultGUI;
% resultGUI = ResultConstRBE.RBEMCNreCalc.resultGUI;


% Resultados para RBExD para el modelo RBEMCN optimizado
% dij = ResultRBEMCN.Optimized.dij;
% resultGUI = ResultRBEMCN.ConstRBEreCalc.resultGUI;
% resultGUI = ResultRBEMCN.Optimized.resultGUI;


% Resultados para RBExD para el modelo RBEUCM optimizado
% dij = ResultRBEUCM.Optimized.dij;
% resultGUI = ResultRBEUCM.ConstRBEreCalc.resultGUI;
% resultGUI = ResultRBEUCM.RBEMCNreCalc.resultGUI;

%% 7 - Apertura de la GUI

matRadGUI
