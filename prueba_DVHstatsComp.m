
%% Función para hacer comparaciones de dosis

% Explicación de las variables
% 1: RBExD - Dosis para calcular las estadísticas
% 2: name0 - Nombre de la tabla
% 3: refVol/refGy - Parametros de volumen y dosis para los que se realizan
%    los cálculos
% 4: HICIAct - Parámetro de activación del cálculo de los coeficientes HI 
%    y CI
% 5: EasyStat - Activación del cálculo de los valores básicos (max, min, media)

function DoseStatistics = prueba_DVHstatsComp (RBExD_1, RBExD_2, RBExD_3, name0, refVol, refGy, pln, cst, HICIAct, EasyStat)

if isempty(refVol) && isempty(refGy)
    refVol = [2 5 50 95 98];
    refGy = linspace(0,max(pln.numOfFractions.*RBExD_3(:)),6);
end

figure('Color',[0.5 0.5 0.5],'Position',([0 30 1500 650]));

% Create the column and row names in cell arrays 
cnames = {'dummy_a'};
rnames = cst(:,2);
% Create the uitable
table = uitable(gcf,'Data',zeros(length(rnames),length(cnames)),...
            'ColumnName',cnames,... 
            'RowName',rnames,'ColumnWidth',{70});
        
pos = get(subplot(3,1,3),'position');
xlabel('Dose x RBE statistics RBEUCM');
set(subplot(3,1,3),'yTick',[])
set(subplot(3,1,3),'xTick',[])

set(table,'units','normalized')
set(table,'position',pos)



% get quality indicators and fill table
% calculate QIs per VOI
for runVoi = 1:size(cst,1)
    
    indices     = cst{runVoi,4}{1};
    numOfVoxels = numel(indices); 
    voiPrint = sprintf('%3d %20s',cst{runVoi,1},cst{runVoi,2}); %String that will print quality indicators
    
    % get Dose, dose is sorted to simplify calculations
    doseInVoi    = sort(pln.numOfFractions.*RBExD_1(indices));
        
    if ~isempty(doseInVoi)
        if EasyStat > 0
            % easy stats
            QI(runVoi).mean = mean(doseInVoi);
            QI(runVoi).std  = std(doseInVoi);
            QI(runVoi).max  = doseInVoi(end);
            QI(runVoi).min  = doseInVoi(1);
            
            voiPrint = sprintf('%s - Mean dose = %5.2f Gy +/- %5.2f Gy (Max dose = %5.2f Gy, Min dose = %5.2f Gy)\n%27s', ...
                voiPrint,QI(runVoi).mean,QI(runVoi).std,QI(runVoi).max,QI(runVoi).min,' ');
        end
        
        DX = @(x) doseInVoi(ceil((100-x)*0.01*numOfVoxels));
        VX = @(x) numel(doseInVoi(doseInVoi >= x)) / numOfVoxels;

        % create VX and DX struct fieldnames at runtime and fill
        for runDX = 1:numel(refVol)
            QI(runVoi).(strcat('D',num2str(refVol(runDX)))) = DX(refVol(runDX));
            voiPrint = sprintf('%sD%d%% = %5.2f Gy, ',voiPrint,refVol(runDX),DX(refVol(runDX)));
        end
        voiPrint = sprintf('%s\n%27s',voiPrint,' ');
        for runVX = 1:numel(refGy)
            sRefGy = num2str(refGy(runVX),3);
            QI(runVoi).(['V' strrep(sRefGy,'.','_') 'Gy']) = VX(refGy(runVX));
            voiPrint = sprintf(['%sV' sRefGy 'Gy = %6.2f%%, '],voiPrint,VX(refGy(runVX))*100);
        end
        voiPrint = sprintf('%s\n%27s',voiPrint,' ');

        % if current voi is a target -> calculate homogeneity and conformity
        if strcmp(cst{runVoi,3},'TARGET') > 0 && HICIAct>0

            % loop over target objectives and get the lowest dose objective 
            referenceDose = inf;
            for runObjective = 1:numel(cst{runVoi,6})
               % check if this is an objective that penalizes underdosing 
               if strcmp(cst{runVoi,6}(runObjective).type,'square deviation') > 0 || strcmp(cst{runVoi,6}(runObjective).type,'square underdosing') > 0
                   referenceDose = (min(cst{runVoi,6}(runObjective).dose,referenceDose));
               end            
            end

            if referenceDose == inf 
                voiPrint = sprintf('%s%s',voiPrint,'Warning: target has no objective that penalizes underdosage, ');
            else
 
                StringReferenceDose = regexprep(num2str(round(referenceDose*100)/100),'\D','_');
                % Conformity Index, fieldname contains reference dose
                VTarget95 = sum(doseInVoi >= 0.95*referenceDose); % number of target voxels recieving dose >= 0.95 dPres
                VTreated95 = sum(pln.numOfFractions.*RBExD_1(:) >= 0.95*referenceDose);  %number of all voxels recieving dose >= 0.95 dPres ("treated volume")
                QI(runVoi).(['CI_' StringReferenceDose 'Gy']) = VTarget95^2/(numOfVoxels * VTreated95); 

                % Homogeneity Index (one out of many), fieldname contains reference dose        
                QI(runVoi).(['HI_' StringReferenceDose 'Gy']) = (DX(5) - DX(95))/referenceDose * 100;

                voiPrint = sprintf('%sCI = %6.4f, HI = %5.2f for reference dose of %3.1f Gy\n',voiPrint,...
                                   QI(runVoi).(['CI_' StringReferenceDose 'Gy']),QI(runVoi).(['HI_' StringReferenceDose 'Gy']),referenceDose);
            end
        end
        fprintf('%s\n',voiPrint);
    else    
    fprintf('%3d %20s - No dose information.\n',cst{runVoi,1},cst{runVoi,2});
    end
    
end

VOI = struct('VOI', cst(:,2))';
names = [fieldnames(VOI); fieldnames(QI)];
DoseStatistics.ConstRBE = cell2struct([struct2cell(VOI); struct2cell(QI)], names, 1);

set(table,'ColumnName',fieldnames(QI));
set(table,'Data',(squeeze(struct2cell(QI)))');

clearvars -except cst pln refVol refGy RBExD_1 RBExD_2 RBExD_3 name0 DoseStatistics HICIAct EasyStat


% Create the column and row names in cell arrays 
cnames = {'dummy_a'};
rnames = cst(:,2);
% Create the uitable
table = uitable(gcf,'Data',zeros(length(rnames),length(cnames)),...
            'ColumnName',cnames,... 
            'RowName',rnames,'ColumnWidth',{70});
        
pos = get(subplot(3,1,2),'position');
xlabel('Dose x RBE statistics RBEMcNamara');
set(subplot(3,1,2),'yTick',[])
set(subplot(3,1,2),'xTick',[])

set(table,'units','normalized')
set(table,'position',pos)

% get quality indicators and fill table
% calculate QIs per VOI
for runVoi = 1:size(cst,1)
    
    indices     = cst{runVoi,4}{1};
    numOfVoxels = numel(indices); 
    voiPrint = sprintf('%3d %20s',cst{runVoi,1},cst{runVoi,2}); %String that will print quality indicators
    
    % get Dose, dose is sorted to simplify calculations
    doseInVoi    = sort(pln.numOfFractions.*RBExD_2(indices));
        
    if ~isempty(doseInVoi)
        
        % easy stats
        if EasyStat > 0
            QI(runVoi).mean = mean(doseInVoi);
            QI(runVoi).std  = std(doseInVoi);
            QI(runVoi).max  = doseInVoi(end);
            QI(runVoi).min  = doseInVoi(1);
            
            voiPrint = sprintf('%s - Mean dose = %5.2f Gy +/- %5.2f Gy (Max dose = %5.2f Gy, Min dose = %5.2f Gy)\n%27s', ...
                voiPrint,QI(runVoi).mean,QI(runVoi).std,QI(runVoi).max,QI(runVoi).min,' ');
        end
        
        DX = @(x) doseInVoi(ceil((100-x)*0.01*numOfVoxels));
        VX = @(x) numel(doseInVoi(doseInVoi >= x)) / numOfVoxels;

        % create VX and DX struct fieldnames at runtime and fill
        for runDX = 1:numel(refVol)
            QI(runVoi).(strcat('D',num2str(refVol(runDX)))) = DX(refVol(runDX));
            voiPrint = sprintf('%sD%d%% = %5.2f Gy, ',voiPrint,refVol(runDX),DX(refVol(runDX)));
        end
        voiPrint = sprintf('%s\n%27s',voiPrint,' ');
        for runVX = 1:numel(refGy)
            sRefGy = num2str(refGy(runVX),3);
            QI(runVoi).(['V' strrep(sRefGy,'.','_') 'Gy']) = VX(refGy(runVX));
            voiPrint = sprintf(['%sV' sRefGy 'Gy = %6.2f%%, '],voiPrint,VX(refGy(runVX))*100);
        end
        voiPrint = sprintf('%s\n%27s',voiPrint,' ');

        % if current voi is a target -> calculate homogeneity and conformity
        if strcmp(cst{runVoi,3},'TARGET') > 0 && HICIAct>0

            % loop over target objectives and get the lowest dose objective 
            referenceDose = inf;
            for runObjective = 1:numel(cst{runVoi,6})
               % check if this is an objective that penalizes underdosing 
               if strcmp(cst{runVoi,6}(runObjective).type,'square deviation') > 0 || strcmp(cst{runVoi,6}(runObjective).type,'square underdosing') > 0
                   referenceDose = (min(cst{runVoi,6}(runObjective).dose,referenceDose));
               end            
            end

            if referenceDose == inf 
                voiPrint = sprintf('%s%s',voiPrint,'Warning: target has no objective that penalizes underdosage, ');
            else
 
                StringReferenceDose = regexprep(num2str(round(referenceDose*100)/100),'\D','_');
                % Conformity Index, fieldname contains reference dose
                VTarget95 = sum(doseInVoi >= 0.95*referenceDose); % number of target voxels recieving dose >= 0.95 dPres
                VTreated95 = sum(pln.numOfFractions.*RBExD_2(:) >= 0.95*referenceDose);  %number of all voxels recieving dose >= 0.95 dPres ("treated volume")
                QI(runVoi).(['CI_' StringReferenceDose 'Gy']) = VTarget95^2/(numOfVoxels * VTreated95); 

                % Homogeneity Index (one out of many), fieldname contains reference dose        
                QI(runVoi).(['HI_' StringReferenceDose 'Gy']) = (DX(5) - DX(95))/referenceDose * 100;

                voiPrint = sprintf('%sCI = %6.4f, HI = %5.2f for reference dose of %3.1f Gy\n',voiPrint,...
                                   QI(runVoi).(['CI_' StringReferenceDose 'Gy']),QI(runVoi).(['HI_' StringReferenceDose 'Gy']),referenceDose);
            end
        end
        fprintf('%s\n',voiPrint);
    else    
    fprintf('%3d %20s - No dose information.\n',cst{runVoi,1},cst{runVoi,2});
    end
    
end

VOI = struct('VOI', cst(:,2))';
names = [fieldnames(VOI); fieldnames(QI)];
DoseStatistics.RBEMCN = cell2struct([struct2cell(VOI); struct2cell(QI)], names, 1);

set(table,'ColumnName',fieldnames(QI));
set(table,'Data',(squeeze(struct2cell(QI)))');
clearvars -except cst pln refVol refGy RBExD_1 RBExD_2 RBExD_3 name0 DoseStatistics HICIAct EasyStat

% Create the column and row names in cell arrays 
cnames = {'dummy_a'};
rnames = cst(:,2);
% Create the uitable
table = uitable(gcf,'Data',zeros(length(rnames),length(cnames)),...
            'ColumnName',cnames,... 
            'RowName',rnames,'ColumnWidth',{70});
        
pos = get(subplot(3,1,1),'position');
xlabel('Dose x RBE statistics ConstRBE');
set(subplot(3,1,1),'yTick',[])
set(subplot(3,1,1),'xTick',[])

set(table,'units','normalized')
set(table,'position',pos)

% get quality indicators and fill table
% calculate QIs per VOI
for runVoi = 1:size(cst,1)
    
    indices     = cst{runVoi,4}{1};
    numOfVoxels = numel(indices); 
    voiPrint = sprintf('%3d %20s',cst{runVoi,1},cst{runVoi,2}); %String that will print quality indicators
    
    % get Dose, dose is sorted to simplify calculations
    doseInVoi    = sort(pln.numOfFractions.*RBExD_3(indices));
        
    if ~isempty(doseInVoi)
        if EasyStat > 0
            % easy stats
            QI(runVoi).mean = mean(doseInVoi);
            QI(runVoi).std  = std(doseInVoi);
            QI(runVoi).max  = doseInVoi(end);
            QI(runVoi).min  = doseInVoi(1);
            
            voiPrint = sprintf('%s - Mean dose = %5.2f Gy +/- %5.2f Gy (Max dose = %5.2f Gy, Min dose = %5.2f Gy)\n%27s', ...
                voiPrint,QI(runVoi).mean,QI(runVoi).std,QI(runVoi).max,QI(runVoi).min,' ');
        end
        
        DX = @(x) doseInVoi(ceil((100-x)*0.01*numOfVoxels));
        VX = @(x) numel(doseInVoi(doseInVoi >= x)) / numOfVoxels;

        % create VX and DX struct fieldnames at runtime and fill
        for runDX = 1:numel(refVol)
            QI(runVoi).(strcat('D',num2str(refVol(runDX)))) = DX(refVol(runDX));
            voiPrint = sprintf('%sD%d%% = %5.2f Gy, ',voiPrint,refVol(runDX),DX(refVol(runDX)));
        end
        voiPrint = sprintf('%s\n%27s',voiPrint,' ');
        for runVX = 1:numel(refGy)
            sRefGy = num2str(refGy(runVX),3);
            QI(runVoi).(['V' strrep(sRefGy,'.','_') 'Gy']) = VX(refGy(runVX));
            voiPrint = sprintf(['%sV' sRefGy 'Gy = %6.2f%%, '],voiPrint,VX(refGy(runVX))*100);
        end
        voiPrint = sprintf('%s\n%27s',voiPrint,' ');

        % if current voi is a target -> calculate homogeneity and conformity
        if strcmp(cst{runVoi,3},'TARGET') > 0 && HICIAct>0

            % loop over target objectives and get the lowest dose objective 
            referenceDose = inf;
            for runObjective = 1:numel(cst{runVoi,6})
               % check if this is an objective that penalizes underdosing 
               if strcmp(cst{runVoi,6}(runObjective).type,'square deviation') > 0 || strcmp(cst{runVoi,6}(runObjective).type,'square underdosing') > 0
                   referenceDose = (min(cst{runVoi,6}(runObjective).dose,referenceDose));
               end            
            end

            if referenceDose == inf 
                voiPrint = sprintf('%s%s',voiPrint,'Warning: target has no objective that penalizes underdosage, ');
            else
 
                StringReferenceDose = regexprep(num2str(round(referenceDose*100)/100),'\D','_');
                % Conformity Index, fieldname contains reference dose
                VTarget95 = sum(doseInVoi >= 0.95*referenceDose); % number of target voxels recieving dose >= 0.95 dPres
                VTreated95 = sum(pln.numOfFractions.*RBExD_3(:) >= 0.95*referenceDose);  %number of all voxels recieving dose >= 0.95 dPres ("treated volume")
                QI(runVoi).(['CI_' StringReferenceDose 'Gy']) = VTarget95^2/(numOfVoxels * VTreated95); 

                % Homogeneity Index (one out of many), fieldname contains reference dose        
                QI(runVoi).(['HI_' StringReferenceDose 'Gy']) = (DX(5) - DX(95))/referenceDose * 100;

                voiPrint = sprintf('%sCI = %6.4f, HI = %5.2f for reference dose of %3.1f Gy\n',voiPrint,...
                                   QI(runVoi).(['CI_' StringReferenceDose 'Gy']),QI(runVoi).(['HI_' StringReferenceDose 'Gy']),referenceDose);
            end
        end
        fprintf('%s\n',voiPrint);
    else    
    fprintf('%3d %20s - No dose information.\n',cst{runVoi,1},cst{runVoi,2});
    end
    
end

VOI = struct('VOI', cst(:,2))';
names = [fieldnames(VOI); fieldnames(QI)];
DoseStatistics.RBEUCM = cell2struct([struct2cell(VOI); struct2cell(QI)], names, 1);

set(table,'ColumnName',fieldnames(QI));
set(table,'Data',(squeeze(struct2cell(QI)))');

title(sprintf('DVH comparation for optimized %s', name0));
hold off

end