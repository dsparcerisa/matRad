%% Funcion para el calculo de los valores de NTCP

% Modelos implementados
% 1- Fukahori (2016)
% 2 - Burman (1991)
% 3 - Cheung (2004) y (2007)
% 4 - Liu (2010)
% 5 - Tucker (2007)
% 6 - Peeters (2006)
% 7 - Shaacke (2016)

function [NTCP] = prueba_NTCPcalc (pln, cst, phantomtype, Dose)

% Modelo LKB
LKM = @(x)(exp (-0.5.*x.^2))*(1/sqrt(2*pi));

% Example LKB model
% Dose = ResultConstRBE.Optimized.resultGUI.RBExD;
% numOfVois = size(cst,1);
%
% for j = 1:numOfVois
%     if strcmpi('Rectum',cst{j,2}) > 0
%         indices     = cst{j,4}{1};
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
    NTCP.Fukahori.VOI = 'Rectum';
    
    numOfVois = size(cst,1);
    for Grade = 1:2
        NTCP.Fukahori.NTCP(Grade).Grade = Grade;
        if Grade == 1
            n = 0.035;
            m = 0.10;
            TD_50 = 63.6; % Gy(RBE)
            
        elseif Grade == 2
            n = 0.012;
            m = 0.046;
            TD_50 = 69.1; % Gy(RBE)
        end
        for j = 1:numOfVois
            if strcmpi('Rectum',cst{j,2}) > 0
                indices     = cst{j,4}{1};
                EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
                t =  (EUD - TD_50)/(m * TD_50);
                NTCP.Fukahori.NTCP(Grade).NTCP = integral(LKM, -inf, t);
            end
        end
    end 
    
    % Burman (1991) -> multiple NTCP 
    NTCP.Burman.Model = 'Burman'' Model';
    
    NTCP.Burman.Rectum.VOI = 'Rectum';
    NTCP.Burman.Rectum.Risk = 'Sever proctitis/ necrosis/ stenosis fistula';
    
    NTCP.Burman.Bladder.VOI = 'Bladder';
    NTCP.Burman.Bladder.Risk = 'Synmptomatic bladder contracture & volume loss';
    
    NTCP.Burman.Femoral_Head.VOI = 'Femoral Head';
    NTCP.Burman.Femoral_Head.Risk = 'Necrosis';
    
    numOfVois = size(cst,1);
    for j = 1:numOfVois
        if strcmpi('Rectum',cst{j,2}) > 0
            n = 0.12;
            m = 0.15;
            TD_50 = 80;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Burman.Rectum.NTCP = integral(LKM, -inf, t);
            
        elseif strcmpi('Bladder',cst{j,2}) > 0
            n = 0.5;
            m = 0.11;
            TD_50 = 80;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Burman.Bladder.NTCP = integral(LKM, -inf, t);
            
            
        elseif strcmpi('Rt femoral head',cst{j,2}) > 0 || strcmpi('Lt femoral head',cst{j,2}) > 0          
            n = 0.25;
            m = 0.12;
            TD_50 = 65;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            if strcmpi('Rt femoral head',cst{j,2}) > 0
                NTCP.Burman.Femoral_Head.NTCP_RT = integral(LKM, -inf, t);
            elseif strcmpi('Lt femoral head',cst{j,2}) > 0
                NTCP.Burman.Femoral_Head.NTCP_LT = integral(LKM, -inf, t);
            end
        end
    end
    
    % Cheung (2004) -> Rectum
    % Cheung (2007) -> Bladder
    NTCP.Cheung.Model = 'Cheung'' Model';
    
    NTCP.Cheung.Rectum.VOI = 'Rectum';
    NTCP.Cheung.Rectum.Risk = 'Rectal bleeding';
    
    NTCP.Cheung.Bladder.VOI = 'Bladder';
    NTCP.Cheung.Bladder.Risk = 'Genitourinary toxicity Grade >= 1';
    
    numOfVois = size(cst,1);
    for j = 1:numOfVois
        if strcmpi('Bladder',cst{j,2}) > 0
            n = 0.5;
            m = 0.11;
            TD_50 = 80;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Cheung.Bladder.NTCP = integral(LKM, -inf, t);
            
        elseif strcmpi('Rectum',cst{j,2}) > 0
            for delta = 1:2
                if delta == 1
                    n = 3.91;
                    m = 0.156;
                    TD_50 = 53.6;
                elseif delta == 2
                    n = 0.746;
                    m = 0.092;
                    TD_50 = 56.7;
                end
                
                indices     = cst{j,4}{1};
                EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
                t =  (EUD - TD_50)/(m * TD_50);
                NTCP.Cheung.Rectum.NTCP(delta).Hemorrhoids_Pres = delta-1;
                NTCP.Cheung.Rectum.NTCP(delta).NTCP = integral(LKM, -inf, t);
            end
        end
    end
        
    
    % Liu (2010) -> NTCP Rectal bleeding   
    NTCP.Liu.Model = 'Liu'' Model';
    NTCP.Liu.Risk = 'Rectal Bleeding Grade 2';
    NTCP.Liu.VOI = 'Rectum';
    
    numOfVois = size(cst,1);
    for j = 1:numOfVois
        if strcmpi('Rectum',cst{j,2}) > 0
            n = 0.12;
            m = 0.15;
            TD_50 = 80;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Liu.NTCP = integral(LKM, -inf, t);
        end
    end
    
    % Tucker (2007)
    NTCP.Tucker.Model = 'Tucker'' Model';
    NTCP.Tucker.Risk = 'Rectal Bleeding Grade 2 RTOG';
    NTCP.Tucker.Grade2Note = 'Moderate diarrhea and colic; bowel movement > 5 times daily; excesive rectal mucus or intermittent bleeding';
    NTCP.Tucker.VOI = 'Rectum';
    
    numOfVois = size(cst,1);
    for j = 1:numOfVois
        if strcmpi('Rectum',cst{j,2}) > 0
            n = 0.08; 
            m = 0.14;
            TD_50 = 78;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Tucker.NTCP = integral(LKM, -inf, t);
        end
    end
    
    % Peeters (2006)
    NTCP.Peeters.Model = 'Peeters'' Model';
    NTCP.Peeters.VOI = 'Rectum';
    
    % bleeding
    NTCP.Peeters.Bleeding.Risk = 'Bleeding';
    numOfVois = size(cst,1);
    for j = 1:numOfVois
        if strcmpi('Rectum',cst{j,2}) > 0
            n = 0.13;
            m = 0.14;
            TD_50 = 81;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Peeters.Bleeding.NTCP = integral(LKM, -inf, t);
        end
    end
    
    % Frequecy increase
    NTCP.Peeters.Dep_Freq_Incr.Risk = 'Deposition Frequecy increase';
    numOfVois = size(cst,1);
    for j = 1:numOfVois
        if strcmpi('Rectum',cst{j,2}) > 0
            n = 0.39;
            m = 0.24;
            TD_50 = 84;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Peeters.Dep_Freq_Incr.NTCP = integral(LKM, -inf, t);
        end
    end
   
        % Fecal Incontinence
    NTCP.Peeters.Fecal_Inc.Risk = 'Fecal Incontinence';
    numOfVois = size(cst,1);
    for j = 1:numOfVois
        if strcmpi('Rectum',cst{j,2}) > 0
            n = 7.48;
            m = 0.46;
            TD_50 = 105;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Peeters.Fecal_Inc.NTCP = integral(LKM, -inf, t);
        end
    end
    
    
    
    % Schaake (2016) -> NTCP rectal bleeding
    
    NTCP.Schaake.Model = 'Shaake''s model';
    NTCP.Schaake.Risk = 'Rectal bleeding grade 2';
    NTCP.Schaake.VOI = 'Rectum';
    
    refGy = 70;
    for i = 1:size(cst,1)
        if strcmpi(cst{i,2},'Rectum')          
                indices     = cst{i,4}{1};
                numOfVoxels = numel(indices);
                doseInVoi = pln.numOfFractions .* Dose(indices);
                VX = @(x) numel(doseInVoi(doseInVoi >= x)) / numOfVoxels;
                V70_anorectum = VX(refGy) * 100;
        end
     end
     for delta = 1:2
         S = -8.09 + 0.32 * V70_anorectum + 1.19 .* (delta-1);
         NTCP.Schaake.NTCP(delta).AnticUse = delta-1;
         NTCP.Schaake.NTCP(delta).NTCP = 1/(1+exp(-S));      
     end
elseif strcmp (phantomtype, 'Head and Neck')>0
    
end

end