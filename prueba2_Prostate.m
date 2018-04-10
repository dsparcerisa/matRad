%% Estructuracion del programa

% 1 - Carga del phantom/paciente
% 2 - Carga de parametros alpha y beta
% 3 - Definicion de restricciones "distintas"
% 4 - Introduccion de los datos basicos
% 5 - Cálculos

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
cst{1,6}.penalty = 50;
cst{1,6}.volume = 50;
cst{1,6}.robustness = 'none';

% PTV_68
cst{6,6}.type = 'square deviation';
cst{6,6}.dose = 68;
cst{6,6}.penalty = 1000;
cst{6,6}.EUD = NaN;
cst{6,6}.volume = NaN;
cst{6,6}.robustness = 'none';
% prueba 9 de marzo
cst{6,6}.type = 'min DVH objective';
cst{6,6}.dose = 64;
cst{6,6}.penalty = 1000;
cst{6,6}.EUD = 0.95;
cst{6,6}.volume = NaN;
cst{6,6}.robustness = 'none';

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

% meta information for treatment plan 
pln.numOfFractions = 30;
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

%% 5 - Cálculos

[~, ResultConstRBE, ResultRBEMCN, ResultRBEUCM, DoseStatistics, NTCP] = prueba_NTCP(cst, pln, ct, phantomtype);