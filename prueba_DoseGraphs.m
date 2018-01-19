function prueba_DoseGraphs (ct, pln, cst, NumBeam, ProfileType, DisplayOption, Result, Model1, Result2, Model2)

figure 

handles.selectedBeam = NumBeam;
handles.ProfileType = ProfileType;
handles.popupDisplayOption = DisplayOption; 
handles.profileOffset = 0;

fileName = [pln.radiationMode '_' pln.machine];
try
    load(fileName);
    SAD = machine.meta.SAD;
catch
    error(['Could not find the following machine file: ' fileName ]);
end


set(gca,'YDir','normal');
ylabel('{\color{black}dose [Gy]}')
cColor ={'black','green','magenta','cyan','yellow','red','blue', 'orange', 'limegreen', 'lavender', 'banana', 'beige', 'wheat'};

% Rotate the system into the beam.
% passive rotation & row vector multiplication & inverted rotation requires triple matrix transpose
rotMat_system_T = transpose(matRad_getRotationMatrix(pln.gantryAngles(handles.selectedBeam),pln.couchAngles(handles.selectedBeam)));

if strcmp(handles.ProfileType,'longitudinal')
    sourcePointBEV = [handles.profileOffset -SAD   0];
    targetPointBEV = [handles.profileOffset  SAD   0];
elseif strcmp(handles.ProfileType,'lateral')
    sourcePointBEV = [-SAD handles.profileOffset   0];
    targetPointBEV = [ SAD handles.profileOffset   0];
end

rotSourcePointBEV = sourcePointBEV * rotMat_system_T;
rotTargetPointBEV = targetPointBEV * rotMat_system_T;

% perform raytracing on the central axis of the selected beam
[~,l,rho,~,ix] = matRad_siddonRayTracer(pln.isoCenter(handles.selectedBeam,:),ct.resolution,rotSourcePointBEV,rotTargetPointBEV,{ct.cube{1}});
d = [0 l .* rho{1}];
% Calculate accumulated d sum.
vX = cumsum(d(1:end-1));


SelectedCube = handles.popupDisplayOption;
if sum(strcmp(SelectedCube,{'physicalDose','effect','RBExD','alpha','beta','RBE', 'physical_vs_RBExD'})) > 0
    Suffix = '';
else
    Idx    = find(SelectedCube == '_');
    Suffix = SelectedCube(Idx:end);
end

% plot counter
Cnt=2;

% plots physical, RBExD and physical vs RBExD
if strcmp(handles.popupDisplayOption, 'physicalDose')
    mPhysDose = Result.(['physicalDose' Suffix]);
    PlotHandles{1} = plot(vX,pln.numOfFractions .* mPhysDose(ix),'color',cColor{1,1},'LineWidth',3); hold on;
    PlotHandles{1,2} ='physicalDose';
    ylabel('dose in [Gy]');
    
elseif strcmp(handles.popupDisplayOption, 'RBExD')
    mRBExD = Result.(['RBExD' Suffix]);
    PlotHandles{2,1} = plot(vX,pln.numOfFractions .* mRBExD(ix),'color',cColor{1,6},'LineWidth',3); hold on;
    PlotHandles{2,2} = strcat( Model1,' RBExD');
    ylabel('Dose [Gy // Gy(RBE)]');
    
    if ~isempty(Model2) && ~isempty(Result2)
        mRBExD2 = Result2.(['RBExD' Suffix]);
        PlotHandles{3,1} = plot(vX,pln.numOfFractions .* mRBExD2(ix),'color',cColor{1,7},'LineWidth',3); hold on;
        PlotHandles{3,2} = strcat( Model2,' RBExD');
    end
    
elseif  strcmp(handles.popupDisplayOption, 'physical_vs_RBExD')
    
    mPhysDose1 = Result.(['physicalDose' Suffix]);
    PlotHandles{1,1} = plot(vX,pln.numOfFractions .* mPhysDose1(ix),'color',cColor{1,1},'LineWidth',3); hold on;
    PlotHandles{1,2} ='physicalDose';
    
    mRBExD = Result.(['RBExD' Suffix]);
    PlotHandles{2,1} = plot(vX,pln.numOfFractions .* mRBExD(ix),'color',cColor{1,6},'LineWidth',3); hold on;
    PlotHandles{2,2} = strcat( Model1,' RBExD');
    ylabel('Dose [Gy // Gy(RBE)]');
    
    if ~isempty(Model2) && ~isempty(Result2)
        mRBExD2 = Result2.(['RBExD' Suffix]);
        PlotHandles{3,1} = plot(vX,pln.numOfFractions .* mRBExD2(ix),'color',cColor{1,7},'LineWidth',3); hold on;
        PlotHandles{3,2} = strcat( Model2,' RBExD');
    end
    
    Cnt=Cnt+2;
end

%%
% asses target coordinates
tmpPrior = intmax;
tmpSize = 0;
for i=1:size(cst,1)
    if strcmp(cst{i,3},'TARGET') && tmpPrior >= cst{i,5}.Priority && tmpSize<numel(cst{i,4}{1})
        linIdxTarget = unique(cst{i,4}{1});
        tmpPrior=cst{i,5}.Priority;
        tmpSize=numel(cst{i,4}{1});
        VOI = cst{i,2};
    end
end

str = sprintf('profile plot - central axis of %d beam gantry angle %d? couch angle %d?',...
    handles.selectedBeam ,pln.gantryAngles(handles.selectedBeam),pln.couchAngles(handles.selectedBeam));

% plot target boundaries
mTargetCube = zeros(ct.cubeDim);
mTargetCube(linIdxTarget) = 1;
vProfile = mTargetCube(ix);
WEPL_Target_Entry = vX(find(vProfile,1,'first'));
WEPL_Target_Exit  = vX(find(vProfile,1,'last'));
PlotHandles{Cnt,2} =[VOI 'boundary'];

mdylim = Result.(['RBExD' Suffix]);
ysuplim = max(mdylim(:))*pln.numOfFractions;  % boundaries line limits

if ~isempty(WEPL_Target_Entry) && ~isempty(WEPL_Target_Exit)
    hold on
    PlotHandles{Cnt,1} = ...
        line([WEPL_Target_Entry WEPL_Target_Entry], [0 ysuplim+10], 'LineStyle','--','Linewidth',3,'color','k');hold on
    line([WEPL_Target_Exit WEPL_Target_Exit], [0 ysuplim+10], 'LineStyle','--','Linewidth',3,'color','k');hold on
    
else
    PlotHandles{Cnt,1} =[];
end

Lines1  = PlotHandles(~cellfun(@isempty,PlotHandles(:,1)),1);
Lines2  = PlotHandles(~cellfun(@isempty,PlotHandles(:,2)),1);
Labels = PlotHandles(~cellfun(@isempty,PlotHandles(:,1)),2);

legend([Lines1{:} Lines2{:}],Labels{:});

xlabel('radiological depth [mm]','FontSize',8);
grid on, grid minor

hold off

end
