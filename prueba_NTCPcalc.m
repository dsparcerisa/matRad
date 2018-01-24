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
    NTCP.Fukahori.VOI = 'Rectum';
    
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
            numOfVoxels = numel(indices);
            doseInVoi   = pln.numOfFractions.*Dose(indices);
            EUD = (1/numel(indices).*sum(Dose(indices).^(1/n))^n);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Burman.Rectum.NTCP = integral(LKM, -inf, t);
            
        elseif strcmpi('Bladder',cst{j,2}) > 0
            n = 0.5;
            m = 0.11;
            TD_50 = 80;
            indices     = cst{j,4}{1};
            numOfVoxels = numel(indices);
            doseInVoi   = pln.numOfFractions.*Dose(indices);
            EUD = (1/numel(indices).*sum(Dose(indices).^(1/n))^n);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Burman.Bladder.NTCP = integral(LKM, -inf, t);
            
            
        elseif strcmpi('Rt femoral head',cst{j,2}) > 0 || strcmpi('Lt femoral head',cst{j,2}) > 0;          
            n = 0.25;
            m = 0.12;
            TD_50 = 65;
            indices     = cst{j,4}{1};
            numOfVoxels = numel(indices);
            doseInVoi   = pln.numOfFractions.*Dose(indices);
            EUD = (1/numel(indices).*sum(Dose(indices).^(1/n))^n);
            t =  (EUD - TD_50)/(m * TD_50);
            if strcmpi('Rt femoral head',cst{j,2}) > 0
                NTCP.Burman.Bladder.NTCP_RT = integral(LKM, -inf, t);
            elseif strcmpi('Lt femoral head',cst{j,2}) > 0
                NTCP.Burman.Bladder.NTCP_LT = integral(LKM, -inf, t);
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
    if strcmpi('Bladder',cst{j,2}) > 0
        n = 0.5;
        m = 0.11;
        TD_50 = 80;
        indices     = cst{j,4}{1};
        numOfVoxels = numel(indices);
        doseInVoi   = pln.numOfFractions.*Dose(indices);
        EUD = (1/numel(indices).*sum(Dose(indices).^(1/n))^n);
        t =  (EUD - TD_50)/(m * TD_50);
        NTCP.Cheung.Bladder.NTCP = integral(LKM, -inf, t);
        
    elseif strcmpi('Rectum',cst{j,2}) > 0
        for delta = 1:2
            if delta == 1
                n = 0.746;
                m = 0.092;
                TD_50 = 56.7;
            elseif delta == 2
                n = 3.91;
                m = 0.156;
                TD_50 = 53.6;
            end
                
            indices     = cst{j,4}{1};
            numOfVoxels = numel(indices);
            doseInVoi   = pln.numOfFractions.*Dose(indices);
            EUD = (1/numel(indices).*sum(Dose(indices).^(1/n))^n);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Cheung.Rectum.NTCP(delta).With_Hemorrhoids = delta-1;
            NTCP.Cheung.Rectum.NTCP(delta).NTCP = integral(LKM, -inf, t);
        end
        
            n = 0.12;
            m = 0.15;
            TD_50 = 80;
            indices     = cst{j,4}{1};
            numOfVoxels = numel(indices);
            doseInVoi   = pln.numOfFractions.*Dose(indices);
            EUD = (1/numel(indices).*sum(Dose(indices).^(1/n))^n);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Burman.Rectum.NTCP = integral(LKM, -inf, t);
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
            numOfVoxels = numel(indices);
            doseInVoi   = pln.numOfFractions.*Dose(indices);
            EUD = (1/numel(indices).*sum(Dose(indices).^(1/n))^n);
            t =  (EUD - TD_50)/(m * TD_50);
            NTCP.Liu.NTCP = integral(LKM, -inf, t);
        end
    end
         
   
    
    
    
    
    
    % Schaake (2016) -> NTCP rectal bleeding
    
    NTCP.Schaake.Model = 'Shaake''s model';
    NTCP.Schaake.Risk = 'Rectal bleeding';
    NTCP.Schaake.VOI = 'Rectum';
    
    refGy = [70];
    for i = 1:size(cst,1)
        if strcmpi(cst{i,2},'Rectum')          
                indices     = cst{i,4}{1};
                numOfVoxels = numel(indices);
                doseInVoi    = pln.numOfFractions.*Dose(indices);
                VX = @(x) numel(doseInVoi(doseInVoi >= x)) / numOfVoxels;
                for runVX = 1:numel(refGy)
                    sRefGy = num2str(refGy(runVX),3);
                    DoseStatistics.V70Gy = VX(refGy(runVX));
                end      
                
            V70_anorectum = DoseStatistics.V70Gy;%
        end
     end
     for delta = 1:2
         S = -8.09 + 0.32 .* V70_anorectum .*100 + 1.19 .* (delta-1);
         NTCP.Schaake.NTCP(delta).AnticUse = delta-1;
         NTCP.Schaake.NTCP(delta).NTCP = 1./(1+exp(-S));      
     end
    
end

end