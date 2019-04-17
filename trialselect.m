function [out_data] = trialselect(s_data, cfg)
% simplify trial selection

if isfield(cfg, 'trial')
    
    out_data = s_data;
    out_data.trialinfo = s_data.trialinfo(cfg.trial, :);
    out_data.trial = s_data.trial(:,:,cfg.trial);
    out_data.velocity = s_data.velocity(:,:,cfg.trial);
    out_data.MSthresh = s_data.MSthresh(:,:,cfg.trial);
    out_data.lgcl_mask_MS = s_data.lgcl_mask_MS(:,cfg.trial);
    out_data.angle = s_data.angle(:,cfg.trial);
    out_data.lgcl_MS_onset = s_data.lgcl_MS_onset(:,cfg.trial);
%    out_data.old_lgcl_mask = s_data.old_lgcl_mask(:,cfg.trial);

else
    
    error('cfg requires a "trial" subfield')
    
end
    
end

