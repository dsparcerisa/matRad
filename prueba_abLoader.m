%% Función para cargar los parámetros alpha y beta

function cst = prueba_abLoader (cst, phantomtype)

%parse variables from base-workspace
cst = evalin('base','cst');

%Datos generales
AvailableAlphaXBetaX = {'BODY', [0.035 0.01],  'Bentzen (1988)(VdK-2009)', 'Arbitrario (musc-vascular)';
    'LIPS', [0.0474  0.0206], 'Emami(1985)', 'Skin(late reactions)';
    'Skin', [0.0474 0.0206], 'Emami(1985)', 'Skin(late reactions)';
    'Heart', [0.0702 0.0351], 'Emami (1985)', 'Heart';
    'Liver', [0.0683 0.0455], 'Emami(1985)+Dörr (2009', 'Liver';
    'Bladder', [0.0878 0.0146],'Emami (1985)', 'Bladder';
    'KidneyL', [0.0106 0.0036], 'Emami(1985)+Beyzadeoglu (2010)', 'Kidney';
    'Stomach', [0.0927 0.01112], 'Emami (1985)', 'Stomach';
    'LargeBowel', [0.0845 0.0141], 'Emami (1985)', 'Small intestine';
    'SmallBowel', [0.0845 0.0141], 'Emami (1985)', 'Small intestine';
    'duodenum', [0.0845 0.0141], 'Emami (1985)', 'Small intestine';
    'Rectum', [0.0484 0.0124], 'Emami(1985)+Fowler (2001)', 'Rectum';
    'PAROTID LT', [0.0628 0.0209], 'Emami(1985)', 'Parotid';
    'LARYNX', [0.0875 0.0208], 'Emami (1985)', 'Larynx';
    'Lt femoral head', [0.03 0.01], 'Sfjro (2010)', 'Bone';
    'Rt femoral head', [0.03 0.01], 'Sfjro (2010)', 'Bone';
    'SPINAL CORD', [0.0445 0.0136], 'Emami (1985) + Dische (1999)(VdK-2009)', 'Spinal Cord';
    'Lymph Nodes', [0.035 0.01], 'Bentzen (1988)(VdK-2009)', 'Arbitrario (musc-vascular)';
    'Penile bulb', [0.035 0.01], 'Bentzen (1988)(VdK-2009)', 'Arbitrario (musc-vascular)';
    'Celiac', [0.035 0.01], 'Bentzen (1988)(VdK-2009)', 'Arbitrario (musc-vascular)';
    'SMASMV', [0.035 0.01], 'Bentzen (1988)(VdK-2009)', 'Arbitrario (musc-vascular)';  
    'BRAIN STEM', [0.0491 0.0234], 'Emami (1985)', 'Brain stem';
    'CEREBELLUM', [0.0499 0.0238], 'Emami (1985)', 'Brain';
    'TEMP LOBE LT', [0.0499 0.0238], 'Emami (1985)', 'Brain';
    'TM JOINT LT', [0.0761 0.0217], 'Emami (1985)', 'Temporomandibular & mandible';
    'CHIASMA', [0.0586 0.0195], 'Emami(1985) +Jiang (1994) (VdK-2009)', 'Optic Chiasma (Eye)'; 
    'LENS LT', [0.0686 0.0572], 'Emami (1985)+Dörr (2009)', 'Eye lens';
    'OPTIC NRV LT', [0.0586 0.0195], 'Emami(1985)', 'Optic nerve'};

%Datos de tumores
if strcmp (phantomtype, 'Prostate')
    TumorAlphaXBetaX = {'PTV', [0.026 0.024], 'Bentzen and Ritter (2005)+Chapman(2015)', 'Prostate cancer';
        'prostate bed', [0.026 0.024], 'Bentzen and Ritter (2005)+Chapman(2015)', 'Prostate Cancer'};
    
else strcmp (phantomtype, 'Head and Neck')
    TumorAlphaXBetaX = {'PTV', [0.4 0.21], 'Björk-Eriksson (2000)', 'H&N carcinoma'};
    
end


for i = 1:size(cst,1)
    % Carga de datos para OAR
    for j=1:size(AvailableAlphaXBetaX,1)
        if strcmpi (cst{i,2},AvailableAlphaXBetaX{j,1}) > 0
            cst{i,5}.alphaX = AvailableAlphaXBetaX{j,2}(1);
            cst{i,5}.betaX = AvailableAlphaXBetaX{j,2}(2);
            
        elseif length (cst{i,2}) >= 5 && length (AvailableAlphaXBetaX{j,1}) >= 5 && strcmpi (cst{i,2}(1:5),AvailableAlphaXBetaX{j,1}(1:5)) > 0
            cst{i,5}.alphaX = AvailableAlphaXBetaX{j,2}(1);
            cst{i,5}.betaX = AvailableAlphaXBetaX{j,2}(2);
        end
    end
    % Carga de datos tumorales
    for k = 1:size(TumorAlphaXBetaX,1)
        if strcmpi (cst{i,2}(1:3),'PTV') > 0 || strcmpi (cst{i,2}(1:3),'GTV') > 0 ||...
                strcmpi (cst{i,2}(1:3),'CTV') > 0 || strcmpi (cst{i,2},TumorAlphaXBetaX{k,1}) > 0
            cst{i,5}.alphaX = TumorAlphaXBetaX{k,2}(1);
            cst{i,5}.betaX = TumorAlphaXBetaX{k,2}(2);
        end
    end
    % Si alguno de los parámetros está vacío elegir los parámetros alpha y
    % beta de una lista
    if isempty (cst{i,5}.alphaX) ||  isempty (cst{i,5}.betaX)
        while val < 1
            fprintf('Data not found for the %s region, please, select one of the next list',cst{i,2});
            fprintf('1 - Musc. Vascular\n 2 - Skin\n 3 - Brain\n 4 - Eye\ 5 - Salival gland\n 6 - Intestine\n ');
            abselect = input('Region selected: ');
            clc
            if abselect == 1;
                cst{i,5}.alphaX = 0.035;
                cst{i,5}.betaX = 0.01;
                val = 1;
            elseif abselect == 2;
                cst{i,5}.alphaX = 0.0474;
                cst{i,5}.betaX = 0.0206;
                val = 1;
            elseif abselect == 3;
                cst{i,5}.alphaX = 0.0499;
                cst{i,5}.betaX = 0.0238;
                val = 1;
            elseif abselect == 4;
                cst{i,5}.alphaX = 0.0586;
                cst{i,5}.betaX = 0.0195;
                val = 1;
            elseif abselect == 5;
                cst{i,5}.alphaX = 0.0628;
                cst{i,5}.betaX = 0.0209;
                val = 1;
            elseif abselect == 6;
                cst{i,5}.alphaX = 0.0845;
                cst{i,5}.betaX = 0.0141;
                val = 1;              

            else
                fprintf('Error!!! Select a correct number.\n\n');
            end  
        end
    end
end
end
