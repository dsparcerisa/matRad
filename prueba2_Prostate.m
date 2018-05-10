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
load PROSTATE.mat; phantomtype = 'Prostate';

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

% PTV_68
cst{6,5}.Priority = 1;
cst{6,6}(1).type = 'square deviation';
cst{6,6}(1).dose = 78;
cst{6,6}(1).penalty = 1000;
cst{6,6}(1).EUD = NaN;
cst{6,6}(1).volume = NaN;
cst{6,6}(1).robustness = 'none';

cst{6,6}(2).type = 'min DVH objective';
cst{6,6}(2).dose = 78;
cst{6,6}(2).penalty = 500;
cst{6,6}(2).EUD = NaN;
cst{6,6}(2).volume = 97;
cst{6,6}(2).robustness = 'none';

cst{6,6}(3).type = 'max DVH objective';
cst{6,6}(3).dose = 83;
cst{6,6}(3).penalty = 100;
cst{6,6}(3).EUD = NaN;
cst{6,6}(3).volume = 1;
cst{6,6}(3).robustness = 'none';

% PTV_56
cst{7,5}.Priority = 2;
cst{7,6}.type = 'square deviation';
cst{7,6}.dose = 56;
cst{7,6}.penalty = 1000;
cst{7,6}.EUD = NaN;
cst{7,6}.volume = NaN;
cst{7,6}.robustness = 'none';

% Bladder
cst{8,5}.Priority = 3;
cst{8,6}(1).type = 'max DVH objective';
cst{8,6}(1).dose = 70;
cst{8,6}(1).penalty = 300;
cst{8,6}(1).EUD = NaN;
cst{8,6}(1).volume = 35;
cst{8,6}(1).robustness = 'none';

cst{8,6}(2).type = 'max DVH objective';
cst{8,6}(2).dose = 50;
cst{8,6}(2).penalty = 300;
cst{8,6}(2).EUD = NaN;
cst{8,6}(2).volume = 60;
cst{8,6}(2).robustness = 'none';

% Rectum
cst{1,5}.Priority = 3;
cst{1,6}(1).type = 'max DVH objective';
cst{1,6}(1).dose = 70;
cst{1,6}(1).penalty = 300;
cst{1,6}(1).EUD = NaN;
cst{1,6}(1).volume = 30;
cst{1,6}(1).robustness = 'none';

cst{1,6}(2).type = 'max DVH objective';
cst{1,6}(2).dose = 50;
cst{1,6}(2).penalty = 300;
cst{1,6}(2).EUD = NaN;
cst{1,6}(2).volume = 50;
cst{1,6}(2).robustness = 'none';

% Femoral heads
cst{4,5}.Priority = 3;
cst{4,6}.type = 'max DVH objective';
cst{4,6}.dose = 50;
cst{4,6}.penalty = 1000;
cst{4,6}.EUD = NaN;
cst{4,6}.volume = 5;
cst{4,6}.robustness = 'none';

cst{10,5}.Priority = 3;
cst{10,6}.type = 'max DVH objective';
cst{10,6}.dose = 50;
cst{10,6}.penalty = 1000;
cst{10,6}.EUD = NaN;
cst{10,6}.volume = 5;
cst{10,6}.robustness = 'none';

% Body
cst{9,5}.Priority = 4;
cst{9,6}.type = 'square overdosing';
cst{9,6}.dose = 30;
cst{9,6}.penalty = 100;
cst{9,6}.EUD = NaN;
cst{9,6}.volume = NaN;
cst{9,6}.robustness = 'none';

% ----------------------------------------

%% 4 - Introduccion de los datos basicos

% meta information for treatment plan 
pln.numOfFractions = 39;
pln.radiationMode = 'protons';           % either photons / protons / helium / carbon
pln.machine = 'Generic';

% beam geometry settings
pln.propStf.bixelWidth     = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles   = [90 270]; % [?];
pln.propStf.couchAngles    = [0 0]; % [?];
pln.propStf.numOfBeams     = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter      = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln.propOpt.runDAO         = false;   % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing  = false;   % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below
%pln.calcLET = true;


%% 5 - Seleccion de calculos

% Seleccion de graficas y estadisticas que mostrar (0 = Desactivado // 1 = Activado)
perfGraphs = 1;       % Graficas de perfil de dosis
perfRBEGraphs = 1;    % Graficas de perfil de dosis vs RBE
DGraphs = 1;          % Graficas de dosis 2D en z = z(dij max)
DVHGraphs = 1;        % Representacion de DVH (1 = Generales // 2 = Especificas)
DVHStats = 1;         % Calculo de las estadisticas generales de dosis

GraphSel = [perfGraphs perfRBEGraphs DGraphs DVHGraphs DVHStats];


% Seleccion de estadisticas V_X y D_X
refVol = []; % Valores D_X
refGy = []; % Valores V_X

StatsRef{1,1} = refVol;
StatsRef{2,1} = refGy;

% Seleccion de modelos de dosis a calcular 
% Dosis{i,j}     j = 1 -> Activacion del calculo (0 = Desactivado // 1 = Activado)
%                j = 2 -> Recalculo de dij y reoptimizacion (0 = Todo // 1 = Solo reoptimizar)

ConstRBE{1,1} = 0;  ConstRBE{1,2} = 0;
RBEMCN{1,1} = 0;    RBEMCN{1,2} = 0;
RBEUCM{1,1} = 0;    RBEUCM{1,2} = 0;

DoseRecalc{1,1} = ConstRBE;
DoseRecalc{2,1} = RBEMCN;
DoseRecalc{3,1} = RBEUCM;

if DVHGraphs == 2 && exist('ResultConstRBE','var') > 0 && exist('ResultRBEMCN', 'var') > 0 && exist('ResultRBEUCM', 'var') > 0
    % Seleccion de regiones especificas para comparar DVH
    % NumReg actua como contador de VOIs Target y OAR.
    % Ej: Si tenemos DVHRegions{i,1} = 'OAR_i' --> NumReg{i,1} = i
    % Ej: Si tenemos DVHRegions{1,1} = 'OAR' y DVHRegions{2,1} = 'Target' --> NumReg{[1,2],1} = 1
    %       Es decir, el numero de OAR y Target se cuenta por separado (NumReg maximo = 3)
    
    DVHRegions{1,1} = 'PTV_68';
    NumReg{1,1} = 1; % Numero de Targets
    DVHRegions{2,1} = 'Rectum';
    NumReg{2,1} = 1; % Numero de OARs
    % DVHRegions{3,1} = 'Bladder';
    % NumReg{3,1} = 2; % Numero de OARs
    
    
    % Seleccion de modelos en comparacion de DVH
    % Dose{i,1} = a.b.resultGUI.RBExD --> Seleccion de la matriz de dosis
    %                               a = ResultConstRBE // ResultRBEMCN // ResultRBEUCM
    %                               b = Optimized // ConstRBEreCalc // RBEMCNreCalc
    
    
    Dose{1,1} = ResultConstRBE.Optimized.resultGUI.RBExD;
    Dose{2,1} = ResultConstRBE.RBEMCNreCalc.resultGUI.RBExD;
    
    % Variable a traves de la que entran todas las opciones al activar DVHGraphs = 2
    CompDVH{1,1} = DVHRegions;
    CompDVH{2,1} = NumReg;
    CompDVH{3,1} = Dose;

end

%% 6 - Calculos
if exist('ResultConstRBE','var') > 0 && exist('ResultRBEMCN', 'var') > 0 && exist('ResultRBEUCM', 'var') > 0
    % Si ya se ha realizado un calculo de todas las matrices de dosis y solo se quiere reevaluar alguna de ellas
    clear DoseResults
    DoseResults{1,1} = ResultConstRBE;
    DoseResults{1,2} = ResultRBEMCN;
    DoseResults{1,3} = ResultRBEUCM;
    
    [~, ResultConstRBE, ResultRBEMCN, ResultRBEUCM, DoseStatistics, NTCP, meanNTCP] = ...
        prueba_NTCP(cst, pln, ct, phantomtype, DoseStatistics, GraphSel, DoseRecalc, DoseResults, StatsRef, CompDVH);
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
    CompDVH{3,1} = [];
    [~, ResultConstRBE, ResultRBEMCN, ResultRBEUCM, DoseStatistics, NTCP, meanNTCP] = ...
        prueba_NTCP(cst, pln, ct, phantomtype, DoseStatistics, GraphSel, DoseRecalc, DoseResults, StatsRef, CompDVH);
end
    
    clearvars -except ct cst phantomtype pln stf ResultRBEMCN ResultRBEUCM ResultConstRBE ResultPhysical midRBE  DoseStatistics NTCP

    
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

% matRadGUI

