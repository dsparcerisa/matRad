%%Función para hacer comparaciones entre DVHs para los diferentes cálculos
%%de dosis

%recomendado que la tercera sea la optimizada

function  prueba_compDVH (Dose1, Dose2, Dose3, pln , cst, name0 ,name1, name2, name3)

if isempty(Dose1)
    Dose1 = zeros (size(Dose3));
    emptyselec = 1;
elseif isempty(Dose2)
    Dose2 = zeros (size(Dose3));
    emptyselec = 2;
else
    emptyselec = 0;
end



figure('Name','DVH','Color',[0.5 0.5 0.5],'Position',([0 30 1500 650]));
hold on
lineStyleIndicator = 1;
numOfVois = size(cst,1);

% calculate and print the dvh
colorMx    = colorcube;
colorMx    = colorMx(1:floor(64/numOfVois):64,:);
lineStyles = {'-',':','--','-.'};

n = 1000;

maxDose1 = max(pln.numOfFractions.*Dose1(:))*1.05;
maxDose2 = max(pln.numOfFractions.*Dose2(:))*1.05;
maxDose3 = max(pln.numOfFractions.*Dose3(:))*1.05;

maxDose = [maxDose1 maxDose2 maxDose3];

dvhPoints = linspace(0,max(maxDose(:))*1.05,n);
dvh       = nan(1,n);
dvhaux    = nan(1,n);
dvhaux2   = nan(1,n);

for i = 1:numOfVois
    if cst{i,5}.Visible
        indices     = cst{i,4}{1};
        numOfVoxels = numel(indices);
        doseInVoi   = pln.numOfFractions.*Dose1(indices);
        doseInVoiaux  = pln.numOfFractions.*Dose2(indices);
        doseInVoiaux2 = pln.numOfFractions.*Dose3(indices);
        
        for j = 1:n
            dvh(j) = sum(doseInVoi > dvhPoints(j));
            dvhaux(j) = sum(doseInVoiaux > dvhPoints(j));
            dvhaux2(j) = sum(doseInVoiaux2 > dvhPoints(j));
        end
        
        dvh = dvh ./ numOfVoxels * 100;
        dvhaux = dvhaux ./ numOfVoxels * 100;
        dvhaux2 = dvhaux2 ./ numOfVoxels * 100;
        
        plot(dvhPoints,dvh,'LineWidth',2,'Color',colorMx(i,:), ...
            'LineStyle','-.','DisplayName',cst{i,2});hold on
        plot(dvhPoints,dvhaux,'LineWidth',2,'Color',colorMx(i,:), ...
            'LineStyle',':','DisplayName',cst{i,2},'HandleVisibility','off');hold on
        plot(dvhPoints,dvhaux2,'LineWidth',2,'Color',colorMx(i,:), ...
            'LineStyle','-','DisplayName',cst{i,2},'HandleVisibility','off');hold on
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

if emptyselec == 1
    title(sprintf('DVH comparation for optimized %s', name0));
    ylabel('Volume [%]','FontSize',fontSizeValue);
    xlabel(sprintf('%s model [RBEGy] (Dashed)  // %s model [RBEGy] (Solid)', name2, name3),'FontSize',fontSizeValue)
    hold off
    
elseif emptyselec == 2
    title(sprintf('DVH comparation for optimized %s', name0));
    ylabel('Volume [%]','FontSize',fontSizeValue);
    xlabel(sprintf('%s model [RBEGy] (Dotted) // %s model [RBEGy] (Solid)', name1, name3),'FontSize',fontSizeValue)
    hold off
    
else
    title(sprintf('DVH comparation for optimized %s', name0));
    ylabel('Volume [%]','FontSize',fontSizeValue);
    xlabel(sprintf('%s model [RBEGy] (Dotted) // %s model [RBEGy] (Dashed)  // %s model [RBEGy] (Solid)', name1, name2, name3),'FontSize',fontSizeValue)
    hold off
end

end