function MSrate = get_MSrate_EXPcondition(s_data, EXPmask, cfg)

% mandatory:
% cfg.fsample = (Hz)
% cfg.movingwin = (s)

dat = s_data.lgcl_mask_MS;

MSrate.allMS = local_computeMSrate(dat, EXPmask, cfg);

if isfield(s_data, 'cell_labeled_MS')
    
    nDirSacc = length(s_data.cell_labeled_MS);
    MSrate.labeled_MS = cell(nDirSacc, 1);
    
    for iDir = 1:nDirSacc
    
        thisdat = s_data.cell_labeled_MS{iDir};
        MSrate.labeled_MS{iDir} = local_computeMSrate(thisdat, EXPmask, cfg);
        
    end
        
end

end


%%

function out = local_computeMSrate(dat, EXPmask, cfg)

this_cond = dat(:, EXPmask);
lwin = round(cfg.fsample*cfg.movingwin);
y = movavg(double(this_cond), 'linear', lwin);
out = mean(y, 2)./cfg.movingwin;

end