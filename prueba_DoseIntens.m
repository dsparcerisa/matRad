function prueba_DoseIntens (ct, pln, Dose, z_cut, TypeDose, Model)
 figure
% Selector de corte de CT para dosis máxima
if isempty(z_cut)
    [~,I] = max(Dose(:));
    [~,~,z_dijmax] = ind2sub(size(Dose),I);
    fprintf('Selected cut off for maximun radiation dose (z = %.f)\n',z_dijmax);
    z_cut = z_dijmax;
end

cd plotting\
% Gráfica de intesidad de dosis
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
cd ..
end
