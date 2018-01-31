
%%Funci�n para hacer comparaciones entre DVHs para los diferentes c�lculos
%%de dosis

%recomendado que la tercera sea la optimizada

function image = prueba_compDVH ( pln , cst, Dose, VOI, OptModel)

figure('Name','DVH','Color',[0.5 0.5 0.5],'Position',([0 30 1500 650]));
hold on
numOfVois = size(cst,1);
LineStyle(1:3,1) = {'-' ':' '--'};

% calculate and print the dvh
colorMx    = colorcube;
colorMx    = colorMx(1:floor(64/(numOfVois)):64,:);

n = 1000;

for i = 1:size(Dose,1)
    maxDose{i,1} = max(Dose{i,1}(:));
    dvh{i,1} = nan(1,n);
end

dvhPoints = linspace (0, pln.numOfFractions.* max([maxDose{i,1}])*1.05, n);




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
                    elseif strcmp (OptModel{i,1}, ' RBEMCNOpt ')
                        LineStyle = ':';
                    elseif strcmp(OptModel{i,1}, ' RBEUCMOpt ')
                        LineStyle = '--';
                    end
                    
                    % Line color
                    if strcmp (cst{j,2}, 'Rectum') && strcmp (OptModel{i,2}, ' ConstRBE ')
                        ColorLine = 'c';
                    elseif strcmp (cst{j,2}, 'Rectum') && strcmp (OptModel{i,2}, ' RBEMCN ')
                        ColorLine = 'b';
                    elseif strcmp (cst{j,2}, 'PTV 68') && strcmp (OptModel{i,2}, ' ConstRBE ')
                        ColorLine = 'm';
                    elseif strcmp (cst{j,2}, 'PTV 68') && strcmp (OptModel{i,2}, ' RBEMCN ')
                        ColorLine = 'r';
                    end
                    
                    plot(dvhPoints,dvh{i,1}, LineStyle,'MarkerSize',2,...
                            'LineWidth',2,'Color', ColorLine,'DisplayName',...
                            strcat(cst{j,2}, OptModel{i,2},' (', OptModel{i,1},')'));hold on
                        
                        
                    % Marcadores espec�ficos
                    
                    % Dosis del PTV
                    if strcmp (cst{j,2}, 'PTV 68')
                        
                        % Restricciones 
                        for m = 1:size(cst{j,6},1)
                            if i==1
                                if strcmp(cst{j,6}(m).type,'max DVH constraint') > 0 || strcmp(cst{j,6}(m).type,'max DVH objective') > 0
                                    line([cst{j,6}(m).dose cst{j,6}(m).dose], [0 110],'DisplayName', 'PTV 68 max DVH objective',...
                                        'LineStyle' , '--' ,'Color','k','LineWidth',2); hold on
                                elseif strcmp(cst{j,6}(m).type,'square deviation') > 0
                                    line([cst{j,6}(m).dose cst{j,6}(m).dose], [0 110],'DisplayName', 'PTV 68 square deviation',...
                                        'LineStyle' , '-.' ,'Color','k','LineWidth',2); hold on
                                elseif strcmp(cst{j,6}(m).type,'min DVH objective') > 0 
                                    line([cst{j,6}(m).dose cst{j,6}(m).dose], [0 110],'DisplayName', 'PTV 68 min DVH objective',...
                                        'LineStyle' , '--' ,'Color',[0.38,0.38,0.38],'LineWidth',2); hold on
                                end
                            end
                        end
                    
                    
                    elseif strcmp (cst{j,2}, 'Rectum')
                                              
                        % V_X en el recto
                        [~, V40] = min(abs(dvhPoints-40));
                        [~, V60] = min(abs(dvhPoints-60));
                        [~, V70] = min(abs(dvhPoints-70));
                        [~, V75] = min(abs(dvhPoints-75));
                        V_Points = [V40 V60 V70 V75];
                        if i==1
                            plot(dvhPoints(V_Points),dvh{i,1}(V_Points),'.', 'MarkerSize',15,'LineWidth',2,...
                                'Color', 'k', 'DisplayName', 'Rectum V40 V60 V70 V75');hold on
                        else
                            plot(dvhPoints(V_Points),dvh{i,1}(V_Points),'.', 'MarkerSize',15,'LineWidth',2,'Color', 'k','HandleVisibility','off');hold on
                        end
                        
                        % Restricciones maximas y m�nimas del DVH
                        for m = 1:size(cst{j,6},1)
                            if i==1
                                if strcmp(cst{j,6}(m).type,'max DVH constraint') > 0 || strcmp(cst{j,6}(m).type,'max DVH objective') > 0
                                    line([cst{j,6}(m).dose cst{j,6}(m).dose], [0 110],'DisplayName', 'Rectum max DVH objective',...
                                        'LineStyle' , '--' ,'Color',[0,0.4,0],'LineWidth',2); hold on
                                elseif strcmp(cst{j,6}(m).type,'min DVH constraint') > 0 || strcmp(cst{j,6}(m).type,'min DVH objective') > 0
                                    line([cst{j,6}(m).dose cst{j,6}(m).dose], [0 110],'DisplayName', 'Rectum min DVH objective',...
                                        'LineStyle' , '--' ,'Color',[0,1,0],'LineWidth',2); hold on
                                end
                            end
                        end
                        
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


fontSizeValue = 11;
myLegend = legend('show','location','NorthWestOutside');
set(myLegend,'FontSize',7,'Interpreter','none');
legend boxoff

ylim([0 110]);
xlim([0 1.2*max(dvhPoints)]);
set(gca,'YTick',0:20:120)

grid on,grid minor
box(gca,'on');
set(gca,'LineWidth',1.5,'FontSize',fontSizeValue);

ylabel('Volume [%]','FontSize',fontSizeValue);

if ~isempty (VOI)
    title(sprintf('DVH comparation'));
    xlabel('Dose [RBEGy]','FontSize',fontSizeValue);
else
    title(sprintf('DVH comparation for %s',OptModel{1,1}));
    if size(Dose,1) == 2
        xlabel(sprintf('Dose [RBEGy]  [%s (Solid) // %s (Dotted)]', OptModel{1,2}, OptModel{2,2}),'FontSize',fontSizeValue);
    elseif size(Dose,1) == 3
        xlabel(sprintf('Dose [RBEGy]  [%s (Solid) // %s (Dotted) // %s (Dashed)]', OptModel{1,2}, OptModel{2,2}, OptModel{3,2}),'FontSize',fontSizeValue);
    end
end

image = figure(1);

end