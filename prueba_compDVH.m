%%Función para hacer comparaciones entre DVHs para los diferentes cálculos
%%de dosis

%recomendado que la tercera sea la optimizada

function  prueba_compDVH ( pln , cst, Dose, VOI, OptModel)

figure('Name','DVH','Color',[0.5 0.5 0.5],'Position',([0 30 1500 650]));
hold on
numOfVois = size(cst,1);
LineStyle(1:6,1) = {'-' ':' '--' '-.' '-o' ':p'};

% calculate and print the dvh
colorMx    = colorcube;
colorMx    = colorMx(1:floor(64/numOfVois):64,:);

n = 1000;

for i = 1:size(Dose,1)
    maxDose{i,1} = max(Dose{i,1}(:));
    dvh{i,1} = nan(1,n);
end

dvhPoints = linspace(0,pln.numOfFractions.*max(maxDose{:})*1.05,n);



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
                    plot(dvhPoints,dvh{i,1}, LineStyle{i,1},'MarkerSize',2,...
                        'LineWidth',2,'Color',colorMx(j,:),'DisplayName',...
                        strcat(cst{j,2}, OptModel{i,2},' (Opt', OptModel{i,1},')'));hold on
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


end