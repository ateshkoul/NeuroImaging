function batch_aal(job)
% ______________________________________________________________________
% 
% Batch mode file for AAL 
% The file takes input from Batch toolbox in spm8 
% It is called by config file tbx_cfg_aal
%
% Written by Atesh Koul, National Brain Research Center, India
% ___________________________________________________________________________
% 

% Contrast from config file
cspec = job.conspec;
% Atlas file from config file
atlas = job.atl;
% Radius for extended local maxima labelling (used only for extended local maxima only)
ext = job.extn;
% Whether to print or not
pr = job.pr;
% Name of the print file
name = job.name;

for k = 1:numel(cspec)
    job.conspec=cspec(k);
    if (numel(cspec(k).contrasts) == 1) && isinf(cspec(k).contrasts)
        tmp    = load(job.spmmat{1});
        cspec1 = repmat(cspec(k),size(tmp.SPM.xCon));
        for l=1:numel(tmp.SPM.xCon)
            cspec1(l).contrasts = l;
        end
        job1           = job;
        job1.conspec   = cspec1;
        batch_aal(job1);  
        
    else
        xSPM.swd       = spm_str_manip(job.spmmat{1},'H');
        xSPM.Ic        = job.conspec.contrasts;
        xSPM.u         = job.conspec.thresh;
        xSPM.Im        = [];
        if ~isempty(job.conspec.mask)
            xSPM.Im    = job.conspec.mask.contrasts;
            xSPM.pm    = job.conspec.mask.thresh;
            xSPM.Ex    = job.conspec.mask.mtype;
        end
        xSPM.thresDesc = job.conspec.threshdesc;
        xSPM.title     = job.conspec.titlestr;
        xSPM.k         = job.conspec.extent;   
       
        switch job.units
            case 1
                xSPM.units = {'mm' 'mm' 'mm'};
            case 2
                xSPM.units = {'mm' 'mm' 'ms'};
            case 3
                xSPM.units = {'mm' 'mm' 'Hz'};
            case 4
                xSPM.units = {'Hz' 'ms' ''};
            case 5
                xSPM.units = {'Hz' 'Hz' ''};
            otherwise
                error('Unknown data type.');
        end
        [hReg xSPM SPM] = spm_results_ui('Setup',xSPM); 
        TabDat = spm_list('List',xSPM,hReg);
        if pr==2 || pr ==4
            spm_figure('Print','Graphics',['Co-or -' name]);
        end        
        assignin('base','TabDat',TabDat);
        assignin('base', 'hReg', hReg);
        assignin('base', 'xSPM', xSPM);
        assignin('base', 'SPM',  SPM);
        switch job.option
            case 1
                %Local Maxima Labeling
                gin_list_dlabels('List',xSPM,atlas,pr); 
            case 2
                %Extended Local Maxima Labeling
                gin_list_plabels('List',xSPM,atlas,ext,pr);
            case 3
                %Cluster Labeling
                gin_clusters_plabels('List',xSPM,atlas,pr);
        end
        
    end
end

