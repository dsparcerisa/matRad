
%%Funcion para hacer comparaciones entre DVHs para los diferentes calculos
%%de dosis

%recomendado que la tercera sea la optimizada

function image = prueba_compDVH ( pln , cst, Dose, VOI, OptModel, VOIType)


figure('Name','DVH','Color',[0.5 0.5 0.5],'Position',([0 30 1500 650]));
hold on
numOfVois = size(cst,1);
LineStyle(1:3,1) = {'-' ':' '--'};

% calculate and print the dvh

if isempty(VOI) > 0
    colorMx    = colorcube;
    colorMx    = colorMx(1:floor(64/(numOfVois)):64,:);
else
    ColorLine    = colorcube;
    ColorLine    = ColorLine(1:floor(64/(size(VOI,1))):64,:);
end

n = 1000;

for i = 1:size(Dose,1)
    maxDose{i,1} = max(Dose{i,1}(:));
    dvh{i,1} = nan(1,n);
end

dvhPoints = linspace (0, pln.numOfFractions.* max([maxDose{i,1}])*1.01, n);


if ~isempty (VOI)
    for i = 1:size(Dose,1)
        for j = 1:numOfVois
            for l = 1:size(VOI,1)
                if cst{j,5}.Visible && strcmpi(VOI(l,1),cst{j,2}) > 0
                    indices     = cst{j,4}{1};
                    numOfVoxels = numel(indices);
                    doseInVoi   = pln.numOfFractions.*Dose{i,1}(indices);
                    for k = 1:n
                        dvh{i,1}(k) = sum(doseInVoi > dvhPoints(k));
                    end
                    dvh{i,1} = dvh{i,1} ./ numOfVoxels * 100;
                                     
                    
                    % Line type
                    if strcmp (OptModel{i,1}, ' ConstRBEOpt ')
                        LineStyle = '-';
                        LegendName = 'Uniform RBE Optimization';
                    elseif strcmp (OptModel{i,1}, ' RBEMCNOpt ')
                        LineStyle = ':';
                        LegendName = 'Variable RBE Optimization';
                    elseif strcmp(OptModel{i,1}, ' RBEUCMOpt ')
                        LineStyle = '--';
                        LegendName = 'MultiRBE Optimization';
                    end
                    
                    if strcmp (OptModel{i,2}, ' ConstRBE ') && strcmpi (VOIType{l,1},'Target') > 0

                    plot(dvhPoints,dvh{i,1}, LineStyle,'MarkerSize',2,...
                            'LineWidth',2,'Color', ColorLine(l,:),'DisplayName',...
                            strcat(cst{j,2}, OptModel{i,2},' (', LegendName,')'));hold on
                        
                    elseif strcmp (OptModel{i,2}, ' RBEMCN ') && strcmpi (VOIType{l,1},'OAR') > 0
                       
                    plot(dvhPoints,dvh{i,1}, LineStyle,'MarkerSize',2,...
                            'LineWidth',2,'Color', ColorLine(l,:),'DisplayName',...
                            strcat(cst{j,2}, OptModel{i,2},' (', LegendName,')'));hold on
                    end
                    
                end
            end
        end
    end
    
else
    for i = 1:size(Dose,1)
        for j = 1:numOfVois
            if cst{j,5}.Visible
                indices     = cst{j,4}{1};
                numOfVoxels = numel(indices);
                doseInVoi   = pln.numOfFractions.*Dose{i,1}(indices);
                for k = 1:n
                    dvh{i,1}(k) = sum(doseInVoi > dvhPoints(k));
                end
                if i>1
                    dvh{i,1} = dvh{i,1} ./ numOfVoxels * 100;
                    plot(dvhPoints,dvh{i,1}, LineStyle{i,1},'MarkerSize',2,...
                        'LineWidth',2,'Color',colorMx(j,:),'DisplayName',cst{j,2},...
                        'HandleVisibility','off');hold on
                else
                    dvh{i,1} = dvh{i,1} ./ numOfVoxels * 100;
                    plot(dvhPoints,dvh{i,1}, LineStyle{i,1},'MarkerSize',2,...
                        'LineWidth',2,'Color',colorMx(j,:),'DisplayName',cst{j,2});hold on
                end
            end
        end
    end
    
end


fontSizeValue = 14;
myLegend = legend('show','location','NorthEast');
set(myLegend,'FontSize',fontSizeValue,'Interpreter','none');
%legend boxoff

ylim([0 105]);
xlim([0 max(dvhPoints)]);
set(gca,'YTick',0:20:120)

grid on,grid minor
box(gca,'on');
set(gca,'LineWidth',1.5,'FontSize',fontSizeValue);

ylabel('Volume [%]','FontSize',fontSizeValue);

if ~isempty (VOI)
    title(sprintf('DVH comparison'));
    xlabel('Dose [RBEGy]','FontSize',fontSizeValue);
else
    title(sprintf('DVH comparison for %s',OptModel{1,1}));
    if size(Dose,1) == 2
        xlabel(sprintf('Dose [RBEGy]  [%s (Solid) // %s (Dotted)]', OptModel{1,2}, OptModel{2,2}),'FontSize',fontSizeValue);
    elseif size(Dose,1) == 3
        xlabel(sprintf('Dose [RBEGy]  [%s (Solid) // %s (Dotted) // %s (Dashed)]', OptModel{1,2}, OptModel{2,2}, OptModel{3,2}),'FontSize',fontSizeValue);
    end
end

image = figure(1);
image.Position = [100 100 900 500];
end