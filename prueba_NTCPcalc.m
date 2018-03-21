%% Funcion para el calculo de los valores de NTCP

% Modelos implementados

% Prostate
% 1 - Fukahori (2016)
% 2 - Burman (1991)
% 3 - Cheung (2004) y (2007)
% 4 - Liu (2010)
% 5 - Tucker (2007)
% 6 - Peeters (2006)
% 7 - Shaacke (2016)

%Head and Neck
% 1 - Semenenko (2008)
% 2 - Burman (1991)
% 3 - Eisbruch (1999)
% 4 - Luxton (2008)
% 5 - Gay (2007)
% 6 - Roesink (2004)
  
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
    
    % Bleeding
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
    
    % Semenenko (2008) -> NTCP xerostomia
    
    NTCP.Semenenko.Model = 'Semenenko''s Model';
    NTCP.Semenenko.Risk = 'Xerostomia';
    NCTP.Semenenko.VOI = 'Parotid gland';
    numOfVois = size(cst,1);
    
    n = 1;
    m = 0.53;
    TD_50 = 31.4;
    
     for j = 1:numOfVois
        if strcmpi('Parotid_LT',cst{j,2}) > 0
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Semenenko.NTCP_Left = integral(LKM, -inf, t);
        elseif strcmpi('Parotid_RT',cst{j,2}) > 0
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Semenenko.NTCP_Right = integral(LKM, -inf, t);
        end
        
    end
    
    % Burman (1991) -> NTCP xerostomia
    
    NTCP.Burman.Model = 'Burman''s Model';
    NTCP.Burman.Risk = 'Xerostomia';
    NTCP.Burman.VOI = 'Parotid gland';
    numOfVois = size(cst,1);
    
    n = 0.7;
    m = 0.18;
    TD_50 = 46;
    
    for j = 1:numOfVois
        if strcmpi('Parotid_LT',cst{j,2}) > 0
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Burman.NTCP_Left = integral(LKM, -inf, t);
        elseif strcmpi('Parotid_RT',cst{j,2}) > 0
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Burman.NTCP_Right = integral(LKM, -inf, t);
        end
    end
    
    % Eisbruch (1999) -> RTOG grade 4 (NTCP Fibrosis)
    
    NTCP.Eisbruch.Model = 'Eisbruch''s Model';
    NTCP.Eisbruch.Risk = 'RTOG grade 4 (fibrosis)';
    NTCP.Eisbruch.VOI = 'Parotid gland';
    numOfVois = size(cst,1);
    
    n = 1;
    m = 0.18;
    TD_50 = 28.4;
    
    for j = 1:numOfVois
        if strcmpi('Parotid_LT',cst{j,2}) > 0
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Eisbruch.NTCP_Left = integral(LKM, -inf, t);
        elseif strcmpi('Parotid_RT',cst{j,2}) > 0
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Eisbruch.NTCP_Right = integral(LKM, -inf, t);
        end
    end

    % Luxton (2008) -> Multiple NTCP
    
    NTCP.Eisbruch.Model = 'Luxton''s Model';
    numOfVois = size(cst,1);
    
    for j = 1:numOfVois
        if strcmpi('Parotid_LT',cst{j,2}) > 0
            NTCP.Luxton.Parotid.Risk = 'Xerostomia';
            NTCP.Luxton.Parotid.VOI = 'Parotid gland';
            n = 0.7;
            m = 0.18;
            TD_50 = 46;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Luxton.Parotid.NTCP_Left = integral(LKM, -inf, t);
        elseif strcmpi('Parotid_RT',cst{j,2}) > 0
            NTCP.Luxton.Parotid.Risk = 'Xerostomia';
            NTCP.Luxton.Parotid.VOI = 'Parotid gland';
            n = 0.7;
            m = 0.18;
            TD_50 = 46;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Luxton.Parotid.NTCP_Right = integral(LKM, -inf, t);
            
        elseif strcmpi('Brain_stem',cst{j,2}) > 0
            NTCP.Luxton.Brain_Stem.Risk = 'Necrosis/Infraction';
            NTCP.Luxton.Brain_Stem.VOI = 'Brain Stem';
            n = 0.16;
            m = 0.15;
            TD_50 = 65;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Luxton.Brain_Stem.NTCP = integral(LKM, -inf, t);
            
        elseif strcmpi('Cerebellum',cst{j,2}) > 0
            NTCP.Luxton.Cerebellum.Risk = 'Necrosis/Infraction';
            NTCP.Luxton.Cerebellum.VOI = 'Brain/Cerebellum';
            n = 0.25;
            m = 0.15;
            TD_50 = 60;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Luxton.Cerebellum.NTCP = integral(LKM, -inf, t);
            
        elseif strcmpi('Larynx',cst{j,2}) > 0
            NTCP.Luxton.Larynx.Necro.Risk = 'Cartilague Necrosis';
            NTCP.Luxton.Larynx.VOI = 'Larynx';
            n_1 = 0.08;
            m_1 = 0.17;
            TD_50_1 = 70;
            indices     = cst{j,4}{1};
            EUD_1 = (sum(Dose(indices).^(1/n_1))/numel(indices))^n_1 *(pln.numOfFractions);
            t_1 =  (EUD_1 - TD_50_1)/(m_1 * TD_50_1);
            NTCP.Luxton.Larynx.Necro.NTCP = integral(LKM, -inf, t_1);
            
            NTCP.Luxton.Larynx.Edema.Risk = 'Laryngeal Edema';
            n_2 = 0.11;
            m_2 = 0.075;
            TD_50_2 = 80;
            indices     = cst{j,4}{1};
            EUD_2 = (sum(Dose(indices).^(1/n_2))/numel(indices))^n_2 *(pln.numOfFractions);
            t_2 =  (EUD_2 - TD_50_2)/(m_2 * TD_50_2);
            NTCP.Luxton.Larynx.Edema.NTCP = integral(LKM, -inf, t_2);
            
        elseif strcmpi('Lens_LT',cst{j,2}) > 0
            NTCP.Luxton.OcularLens.Risk = 'Cataract Requiring Intervention';
            NTCP.Luxton.OcularLens.VOI = 'Ocular Lens';
            n = 0.3;
            m = 0.27;
            TD_50 = 18;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Luxton.OcularLens.NTCP_Left = integral(LKM, -inf, t);
        elseif strcmpi('Lens_RT',cst{j,2}) > 0
            NTCP.Luxton.Parotid.Risk = 'Cataract Requiring Intervention';
            NTCP.Luxton.Parotid.VOI = 'Ocular Lens';
            n = 0.3;
            m = 0.27;
            TD_50 = 18;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Luxton.OcularLens.NTCP_Right = integral(LKM, -inf, t);
            
        elseif strcmpi('Optic_NRV_LT',cst{j,2}) > 0
            NTCP.Luxton.OpticNerve.Risk = 'Blindness';
            NTCP.Luxton.OpticNerve.VOI = 'Optic Nerve';
            n = 0.25;
            m = 0.14;
            TD_50 = 65;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Luxton.OpticNerve.NTCP_Left = integral(LKM, -inf, t);
        elseif strcmpi('Optic_NRV_RT',cst{j,2}) > 0
            NTCP.Luxton.OpticNerve.Risk = 'Blindness';
            NTCP.Luxton.OpticNerve.VOI = 'Optic Nerve';
            n = 0.25;
            m = 0.14;
            TD_50 = 65;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Luxton.OpticNerve.NTCP_Right = integral(LKM, -inf, t);
            
        elseif strcmpi('Chiasma',cst{j,2}) > 0
            NTCP.Luxton.Chiasma.Risk = 'Blindness';
            NTCP.Luxton.Chiasma.VOI = 'Optic Chiasma';
            n = 0.25;
            m = 0.14;
            TD_50 = 65;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Luxton.Chiasma.NTCP = integral(LKM, -inf, t);
            
        elseif strcmpi('Skin',cst{j,2}) > 0
            NTCP.Luxton.Skin.Risk = 'Necrosis/Ulceration';
            NTCP.Luxton.Skin.VOI = 'Skin';
            n = 0.1;
            m = 0.12;
            TD_50 = 70;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Luxton.Skin.NTCP = integral(LKM, -inf, t);
            
        elseif strcmpi('Spinal_Cord',cst{j,2}) > 0
            NTCP.Luxton.Spinal_Cord.Risk = 'Myelitis/Necrosis';
            NTCP.Luxton.Spinal_Cord.VOI = 'Spinal Cord';
            n = 0.05;
            m = 0.175;
            TD_50 = 66.5;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Luxton.Spinal_Cord.NTCP = integral(LKM, -inf, t);
            
        elseif strcmpi('TM_JOINT_LT',cst{j,2}) > 0
            NTCP.Luxton.Tm_Joint.Risk = 'Marked Limit. Joint Funct. ';
            NTCP.Luxton.Tm_Joint.VOI = 'Temporomandibular Joint';
            n = 0.07;
            m = 0.1;
            TD_50 = 72;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Luxton.Tm_Joint.NTCP_Left = integral(LKM, -inf, t);
        elseif strcmpi('TM_JOINT_RT',cst{j,2}) > 0
            NTCP.Luxton.Tm_Joint.Risk = 'Marked Limit. Joint Funct. ';
            NTCP.Luxton.Tm_Joint.VOI = 'Temporomandibular Joint';
            n = 0.07;
            m = 0.1;
            TD_50 = 72;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Luxton.Tm_Joint.NTCP_Right = integral(LKM, -inf, t);
            
        end
    end
    
    % Gay (2007) -> Multiple NTCP
    
    NTCP.Gay.Model = 'Gay''s Model';
    numOfVois = size(cst,1);
    
    for j = 1:numOfVois
        if strcmpi('Brain_stem',cst{j,2}) > 0
            NTCP.Gay.Brain_Stem.Risk = 'Necrosis';
            NTCP.Gay.Brain_Stem.VOI = 'Brain Stem';
            a = 7;
            gamma_50 = 3;
            TD_50 = 65;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^a)/numel(indices))^(1/a) *(pln.numOfFractions);
            NTCP.Gay.Brain_Stem.NTCP = 1/(1 + (TD_50/EUD)^(4 * gamma_50));
            
        elseif strcmpi('Cerebellum',cst{j,2}) > 0
            NTCP.Gay.Cerebellum.Risk = 'Necrosis/';
            NTCP.Gay.Cerebellum.VOI = 'Brain/Cerebellum';
            a = 5;
            gamma_50 = 3;
            TD_50 = 60;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^a)/numel(indices))^(1/a) *(pln.numOfFractions);
            NTCP.Gay.Cerebellum.NTCP = 1/(1 + (TD_50/EUD)^(4 * gamma_50));
            
        elseif strcmpi('Chiasma',cst{j,2}) > 0
            NTCP.Gay.Chiasma.Risk = 'Blindness';
            NTCP.Gay.Chiasma.VOI = 'Optic Chiasma';
            a = 25;
            gamma_50 = 3;
            TD_50 = 65;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^a)/numel(indices))^(1/a) *(pln.numOfFractions);
            NTCP.Gay.Chiasma.NTCP = 1/(1 + (TD_50/EUD)^(4 * gamma_50));
            
        elseif strcmpi('Optic_NRV_LT',cst{j,2}) > 0
            NTCP.Gay.OpticNerve.Risk = 'Blindness';
            NTCP.Gay.OpticNerve.VOI = 'Optic Nerve';
            a = 25;
            gamma_50 = 3;
            TD_50 = 65;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^a)/numel(indices))^(1/a) *(pln.numOfFractions);
            NTCP.Gay.OpticNerve.NTCP_Left = 1/(1 + (TD_50/EUD)^(4 * gamma_50));
        elseif strcmpi('Optic_NRV_RT',cst{j,2}) > 0
            NTCP.Gay.OpticNerve.Risk = 'Blindness';
            NTCP.Gay.OpticNerve.VOI = 'Optic Nerve';
            a = 25;
            gamma_50 = 3;
            TD_50 = 65;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^a)/numel(indices))^(1/a) *(pln.numOfFractions);
            NTCP.Gay.OpticNerve.NTCP_Right = 1/(1 + (TD_50/EUD)^(4 * gamma_50));
            
        elseif strcmpi('Lens_LT',cst{j,2}) > 0
            NTCP.Gay.OcularLens.Risk = 'Cataract Requiring Intervention';
            NTCP.Gay.OcularLens.VOI = 'Ocular Lens';
            a = 3;
            gamma_50 = 1;
            TD_50 = 18;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^a)/numel(indices))^(1/a) *(pln.numOfFractions);
            NTCP.Gay.OcularLens.NTCP_Left = 1/(1 + (TD_50/EUD)^(4 * gamma_50));
        elseif strcmpi('Lens_RT',cst{j,2}) > 0
            NTCP.Gay.OcularLens.Risk = 'Cataract Requiring Intervention';
            NTCP.Gay.OcularLens.VOI = 'Ocular Lens';
            a = 3;
            gamma_50 = 1;
            TD_50 = 18;
            indices     = cst{j,4}{1};
            EUD = (sum(Dose(indices).^a)/numel(indices))^(1/a) *(pln.numOfFractions);
            NTCP.Gay.OcularLens.NTCP_Right = 1/(1 + (TD_50/EUD)^(4 * gamma_50));
            
        end
    end
    
        % Roesink (2004) -> Salivary Excretion
    
        NTCP.Roesink.Model = 'Roesink''s Model';
        NTCP.Roesink.Risk = 'Xerostomia';
        NCTP.Roesink.VOI = 'Parotid gland';
        NTCP.Roesink.Note = 'SEFX = Salivary Excretion Factor < X';
        numOfVois = size(cst,1);
        
        for j = 1:numOfVois
            if strcmpi('Parotid_LT',cst{j,2}) > 0
                n = 1;
                m = 0.42;
                TD_50 = 52;
                indices     = cst{j,4}{1};
                EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
                t =  (EUD - TD_50)/(m * TD_50);
                NTCP.Roesink.Left.NTCP_SEF25 = integral(LKM, -inf, t);   
                clear n m TD_50 indices EUD t
                n = 1;
                m = 0.43;
                TD_50 = 47;
                indices     = cst{j,4}{1};
                EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
                t =  (EUD - TD_50)/(m * TD_50);
                NTCP.Roesink.Left.NTCP_SEF35 = integral(LKM, -inf, t);
                clear n m TD_50 indices EUD t
                n = 1;
                m = 0.53;
                TD_50 = 29;
                indices     = cst{j,4}{1};
                EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
                t =  (EUD - TD_50)/(m * TD_50);
                NTCP.Roesink.Left.NTCP_SEF45 = integral(LKM, -inf, t);
                clear n m TD_50 indices EUD t
                n = 1;
                m = 0.59;
                TD_50 = 25;
                indices     = cst{j,4}{1};
                EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
                t =  (EUD - TD_50)/(m * TD_50);
                NTCP.Roesink.Left.NTCP_SEF55 = integral(LKM, -inf, t);
                clear n m TD_50 indices EUD t
                
                
            elseif strcmpi('Parotid_RT',cst{j,2}) > 0
                n = 1;
                m = 0.42;
                TD_50 = 52;
                indices     = cst{j,4}{1};
                EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
                t =  (EUD - TD_50)/(m * TD_50);
                NTCP.Roesink.Right.NTCP_SEF25 = integral(LKM, -inf, t);   
                clear n m TD_50 indices EUD t
                n = 1;
                m = 0.43;
                TD_50 = 47;
                indices     = cst{j,4}{1};
                EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
                t =  (EUD - TD_50)/(m * TD_50);
                NTCP.Roesink.Right.NTCP_SEF35 = integral(LKM, -inf, t);
                clear n m TD_50 indices EUD t
                n = 1;
                m = 0.53;
                TD_50 = 29;
                indices     = cst{j,4}{1};
                EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
                t =  (EUD - TD_50)/(m * TD_50);
                NTCP.Roesink.Right.NTCP_SEF45 = integral(LKM, -inf, t);
                clear n m TD_50 indices EUD t
                n = 1;
                m = 0.59;
                TD_50 = 25;
                indices     = cst{j,4}{1};
                EUD = (sum(Dose(indices).^(1/n))/numel(indices))^n *(pln.numOfFractions);
                t =  (EUD - TD_50)/(m * TD_50);
                NTCP.Roesink.Right.NTCP_SEF55 = integral(LKM, -inf, t);
                clear n m TD_50 indices EUD t
            end

        end
end

end