%% Función para los cálculos de los daños NTCP

function [NTCP] = prueba_NTCPcalc (pln, cst, phantomtype, Dose)


% Modelo LKB
LKM = @(x)(exp (-0.5.*x.^2/2));

% Example LKB model
% Dose = ResultConstRBE.Optimized.resultGUI.RBExD;
% numOfVois = size(cst,1);
%
% for j = 1:numOfVois
%     if strcmpi('Rectum',cst{j,2}) > 0
%         indices     = cst{j,4}{1};
%         numOfVoxels = numel(indices);
%         doseInVoi   = pln.numOfFractions.*Dose(indices);
%         EUD = (1/numel(indices).*sum(Dose(indices).^(1/n))^n);
%         t =  (EUD - TD_50)/(m * TD_50);
%         NTCP = integral(LKM, -inf, t);
%     end
% end

%%

if strcmp (phantomtype, 'Prostate')>0
    
    % Fukahori (2016) -> NTCP Rectal bleeding   
    NTCP.Fukahori.Model = 'Fukahori'' Model';
    NTCP.Fukahori.Risk = 'Rectal Bleeding';
    
    numOfVois = size(cst,1);
    for Grade = 1:2
        NTCP.Fukahori.NTCP(Grade).Grade = Grade;
        if Grade == 1;
            n = 0.035;
            m = 0.10;
            TD_50 = 63.6; % Gy(RBE)
            
        elseif Grade == 2;
            n = 0.012;
            m = 0.046;
            TD_50 = 69.1; % Gy(RBE)
        end
        for j = 1:numOfVois
            if strcmpi('Rectum',cst{j,2}) > 0
                indices     = cst{j,4}{1};
                numOfVoxels = numel(indices);
                doseInVoi   = pln.numOfFractions.*Dose(indices);
                EUD = (1/numel(indices).*sum(Dose(indices).^(1/n))^n);
                t =  (EUD - TD_50)/(m * TD_50);
                NTCP.Fukahori.NTCP(Grade).NTCP = integral(LKM, -inf, t);
            end
        end
    end 
    
%     % Schaake (2016) -> NTCP rectal bleeding
%     refVol = [70];
%     refGy = [];
%     
%      for i = 1:size(cst,1)
%             if strcmpi(cst{i,2},'Rectum')
%                 NTCP_final.Model = 'Shaake''s model';
%                 NTCP_final.Risk = 'Rectal bleeding';
%                 DoseStatistics = prueba_DVHstatsComp(pln, cst, Dose,[],[], Models, 'ConstRBEOpt',1);
%                 
%                 V70_anorectum = DoseStats.ConstRBE.V70Gy;% Aquí faltan los calculos salidos de prueba_DVHstatsComp
%             end
%      end
%       for delta = 1:2
%             S = -8.09 + 0.32 .* V70_anorectum .*100 + 1.19 .* (delta-1);
%             NTCP(delta).NTCP = 1./(1+exp(-S));
%             NTCP_final.Calc(delta).AnticUse = delta-1;
%             NTCP_final.Calc(delta).NTCP_ConstRBE = NTCP(delta).NTCP;
%         end
%     
% end

end