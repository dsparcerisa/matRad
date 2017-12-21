%% Función de recálculo de dosis en un modelo diferente respecto a uno ya optimizado.

function [dijreCalc,resultGUIreCalc] = prueba_RecalcDose (resultGUI, ct, stf, pln, cst, bioOpt)
% recalculation only makes sense if ...
evalin('base','exist(''pln'',''var'')') && ...
    evalin('base','exist(''stf'',''var'')') && ...
    evalin('base','exist(''ct'',''var'')') && ...
    evalin('base','exist(''cst'',''var'')') && ...
    evalin('base','exist(''resultGUI'',''var'')')

% get all data from workspace
pln       = evalin('base','pln');
stf       = evalin('base','stf');
ct        = evalin('base','ct');
cst       = evalin('base','cst');
resultGUI = evalin('base','resultGUI');

pln.bioOptimization = bioOpt;
% retrieve model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,pln.bioOptimization);
% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,pln.scenGenType);

SelectedCube = 'RBExD';
Suffix = strsplit(SelectedCube,'_');
if length(Suffix)>1
    Suffix = ['_' Suffix{2}];
else
    Suffix = '';
end

if sum([stf.totalNumOfBixels]) ~= length(resultGUI.(['w' Suffix]))
    warndlg('weight vector does not corresponding to current steering file');
    return
end

% change isocenter if that was changed and do _not_ recreate steering
% information
for i = 1:numel(pln.gantryAngles)
    stf(i).isoCenter = pln.isoCenter(i,:);
end

% recalculate influence matrix
dij = matRad_calcParticleDose(ct,stf,pln,cst);
dijreCalc = dij;

% recalculate cubes in resultGUI
resultGUIreCalc = matRad_calcCubes(resultGUI.(['w' Suffix]),dij,cst);
% overwrite the "standard" fields
sNames = fieldnames(resultGUIreCalc);
for j = 1:length(sNames)
    resultGUI.(sNames{j}) = resultGUIreCalc.(sNames{j});
end

end
