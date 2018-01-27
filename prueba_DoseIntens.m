function prueba_DoseIntens (ct, pln, Dose, IsoDose_Levels, z_cut, TypeDose, Model)
figure
% Selector de corte de CT para dosis maxima
if isempty(z_cut)
    [~,I] = max(Dose(:));
    [~,~,z_dijmax] = ind2sub(size(Dose),I);
    fprintf('Selected cut off for maximun radiation dose (z = %.f)\n',z_dijmax);
    z_cut = z_dijmax;
end

if isunix > 1
    cd plotting
end

% Grafica de intesidad de dosis
axis off
title(sprintf('%s map for %s model', TypeDose, Model));
axesHandle = axes;
set(gca,'YDir','Reverse'); %Y axis inversion

matRad_plotCtSlice(axesHandle,ct.cube,1,3,z_cut);

Total_dose = pln.numOfFractions .* Dose;
cmin = min(Total_dose(:));
cmax = max(Total_dose(:));
colormap(axesHandle, jet(64));
caxis([cmin cmax]);
matRad_plotDoseSlice(axesHandle,Total_dose,3,z_cut,0.01,0.8,jet(64), [cmin cmax+0.05]);
axis off

% Lineas de Isodosis
if isempty(IsoDose_Levels)
    IsoDose_Levels = pln.numOfFractions .* [0.18 0.37 0.55 0.74 0.92 1.1 1.3 1.5 1.7 1.8 1.9 2 2.1 2.2];
end

IsoDose_Contours = matRad_computeIsoDoseContours(Total_dose ,IsoDose_Levels);

matRad_plotIsoDoseLines(axesHandle,Total_dose,IsoDose_Contours,IsoDose_Levels,...
   0,3,z_cut,jet(64),[cmin cmax+0.05],'LineWidth',1.5);


if isunix > 1
    cd ..
end

end
