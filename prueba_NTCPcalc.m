%% Funci�n para los c�lculos de los da�os NTCP

function [NTCP_final] = prueba_NTCPcalc (cst, phantomtype, DoseStats_1, DoseStats_2, DoseStats_3, UCMcalc)

if  UCMcalc == 1
    if strcmp (phantomtype, 'Prostate')>0
        
        % Schaake (2016) -> NTCP rectal bleeding
        for i = 1:size(cst,1)
            if strcmpi(cst{i,2},'Rectum')
                NTCP_final.Model = 'Shaake''s model';
                NTCP_final.Risk = 'Rectal bleeding';
                V70_anorectum(1) = DoseStats_1.ConstRBE.V70Gy;
                V70_anorectum(2) = DoseStats_2.RBEMCN.V70Gy;
                V70_anorectum(3)= DoseStats_3.RBEUCM.V70Gy;
                
            end
        end
        
        for delta = 1:2
            S = -8.09 + 0.32 .* V70_anorectum .*100 + 1.19 .* (delta-1);
            NTCP(delta).NTCP = 1./(1+exp(-S));
            NTCP_final.Calc(delta).AnticUse = delta-1;
            NTCP_final.Calc(delta).NTCP_ConstRBE = NTCP(delta).NTCP(1);
            NTCP_final.Calc(delta).NTCP_RBEMCN = NTCP(delta).NTCP(2);
            NTCP_final.Calc(delta).NTCP_RBEUCM = NTCP(delta).NTCP(3);
        end
        
    end  
else
     if strcmp (phantomtype, 'Prostate')>0
        
        % Schaake (2016) -> NTCP rectal bleeding
        for i = 1:size(cst,1)
            if strcmpi(cst{i,2},'Rectum')
                NTCP_final.Model = 'Shaake''s model';
                NTCP_final.Risk = 'Rectal bleeding';
                V70_anorectum(1) = DoseStats_1.ConstRBE.V70Gy;
                V70_anorectum(2) = DoseStats_2.RBEMCN.V70Gy;                
            end
        end
        
        for delta = 1:2
            S = -8.09 + 0.32 .* V70_anorectum + 1.19 .* (delta-1);
            NTCP(delta).NTCP = 1./(1+exp(-S));
            NTCP_final.Calc(delta).AnticUse = delta-1;
            NTCP_final.Calc(delta).NTCP_ConstRBE = NTCP(delta).NTCP(1);
            NTCP_final.Calc(delta).NTCP_RBEMCN = NTCP(delta).NTCP(2);
        end
        
    end     
        
end



%     
% elseif strcmp (phantomtype, 'Head and Neck')>0
%     
%     % Beetz (2012)-> NTCP xerostomia
%     for xerostomia_score = 1:5;
%         for i = 1:size(cst,1)
%             if strcmp(cst{i,2},'PAROTID LT')
%                 mean_dolse_parotidLT(1) =  DoseStatistics.ConstRBE.mean;
%                 mean_dolse_parotidLT(2) =  DoseStatistics.RBEMCN.mean;
%                 mean_dolse_parotidLT(3) =  DoseStatistics.RBEUCM.mean;
%             end
%             if strcmp(cst{i,2},'PAROTID RT')
%                 mean_dolse_parotidRT(1) =  DoseStatistics.ConstRBE.mean;
%                 mean_dolse_parotidRT(2) =  DoseStatistics.RBEMCN.mean;
%                 mean_dolse_parotidRT(3) =  DoseStatistics.RBEUCM.mean;
%             end
%         end
%         S.LT = -1.443 + mean_dolse_parotidLT.*0.047 + xerostomia_score.*0.720;
%         S.RT = -1.443 + mean_dolse_parotidRT.*0.047 + xerostomia_score.*0.720;
%         
%         NTCP(xerostomia_score).xeroscore = xerostomia_score;
%         NTCP(xerostomia_score).LT = 1./(1+exp(-S.LT));
%         NTCP(xerostomia_score).RT = 1./(1+exp(-S.RT));
%         
%         %         fprintf('According to Beetz''s model, the NTCP for xerostomia with score = %.f is:\n',NTCP(xerostomia_score).xeroscore);
%         %         fprintf('In the LT Parotid:\n For a constant RBE: %.8f %%\n For McNamara''s RBE: %.8f %%\n For the UCM''s RBE: %.8f %%\n',...
%         %             NTCP(xerostomia_score).LT(1)*100,NTCP(xerostomia_score).LT(2)*100,NTCP(xerostomia_score).LT(3)*100');
%         %         fprintf('In the RT Parotid:\n For a constant RBE: %.8f %%\n For McNamara''s RBE: %.8f %%\n For the UCM''s RBE: %.8f %%\n\n',...
%         %             NTCP(xerostomia_score).RT(1)*100,NTCP(xerostomia_score).RT(2)*100,NTCP(xerostomia_score).RT(3)*100');
%         
%         NTCP_final(xerostomia_score).xeroscore = xerostomia_score;
%         NTCP_final(xerostomia_score).NTCP_ConstRBE = NTCP(xerostomia_score).NTCP(1);
%         NTCP_final(xerostomia_score).NTCP_RBEMCN = NTCP(xerostomia_score).NTCP(2);
%         NTCP_final(xerostomia_score).NTCP_RBEUCM = NTCP(dxerostomia_score).NTCP(3);
%     end
% end

end