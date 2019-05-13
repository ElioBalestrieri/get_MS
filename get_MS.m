function [s_data] = get_MS(s_data, cfg)
% implement MS detection algorithm described in Engbert & Kliegl (2003)
% currently the only difference in th ealgorithm per se is that, since we
% are collecting data from just one eye, we'll skip the selection of
% microsaccadic event happening for both eyes.
% started by EB on 15-Apr-2019

% subfields necessary
% cfg.screensize = [1920 1080]; 
% cfg.pix_per_va = 51;
% cfg.kernellength = 3;

% step0 convert everything from pixels to va
xchan_idx = find(ismember(s_data.label, 'X_pos'));
ychan_idx = find(ismember(s_data.label, 'Y_pos'));
s_data.trial(:,xchan_idx,:) = (s_data.trial(:,xchan_idx,:)-...
    cfg.screensize(1)/2)/cfg.pix_per_va;
s_data.trial(:,ychan_idx,:) = (s_data.trial(:,ychan_idx,:)-...
    cfg.screensize(2)/2)/cfg.pix_per_va;

% step1 compute velocity
s_data = local_compute_velocity(s_data, xchan_idx, ychan_idx);

% step2 compute threshold
s_data = local_compute_threshold(s_data);

% step 3 select MS events
s_data = local_select_MS(s_data, cfg);

% OPTIONAL step 5 consider only time of interest
if isfield(cfg, 'toi')
    
    s_data = local_select_toi(s_data, cfg);
    
end

% step 6 -added by me-
% a MS is an event extended in time. But having logical TRUE for each
% timepoint showing high velocity might be a confounder, especially when we
% consider the angles: for this reason, we add another subfield containing
% only the first sample of the saccadic event

swap_single_saccade_ON = diff(s_data.lgcl_mask_MS)==1;
s_data.lgcl_MS_onset = logical(cat(1, s_data.lgcl_mask_MS(1,:),...
    swap_single_saccade_ON));

% it has to be noted that the single MS onset might be werid and not than
% informative. Solution: take even the MS offset, and compute the angle
% between the onset and offset of MS, on the data and not on the speed.
swap_single_saccade_OFF = diff(s_data.lgcl_mask_MS) == -1;
s_data.lgcl_MS_offset = logical(cat(1, swap_single_saccade_OFF, ...
    s_data.lgcl_mask_MS(end,:))); % instead of assigning 0s, let's assign the
                                  % last row of the MS, in order to correct
                                  % for MSs started but not conluded on
                                  % time

% ... and only now we compute angles
s_data = local_compute_angles(s_data, xchan_idx, ychan_idx);


end

%% ######################### LOCAL FUNCTIONS ##############################

function s_data = local_compute_velocity(s_data, xchan_idx, ychan_idx)

% create matrices of values for velocity computation
mat_nminus1 = cat(1, nan(1,2,size(s_data.trial, 3)), ...
    s_data.trial(1:end-1, [xchan_idx, ychan_idx], :));
mat_nminus2 = cat(1, nan(2,2,size(s_data.trial, 3)), ...
    s_data.trial(1:end-2, [xchan_idx, ychan_idx], :));
mat_nplus1 = cat(1, s_data.trial(2:end, [xchan_idx, ychan_idx], :),...
    nan(1,2,size(s_data.trial, 3)));
mat_nplus2 = cat(1, s_data.trial(3:end, [xchan_idx, ychan_idx], :),...
    nan(2,2,size(s_data.trial, 3)));

deltaT = mean(diff(s_data.x_time));
s_data.velocity = (mat_nplus2 + mat_nplus1 - mat_nminus1 - mat_nminus2)/...
    (6*deltaT);

end

function s_data = local_compute_threshold(s_data)

% substep 1 --> compute std by applying a median estimator to the TS
% big question: mismatch between 2003 (no square root) and 2006 (sqrt)
% articles. makes more sense to me to have root squared vals though
stdxy = sqrt(nanmedian(s_data.velocity.^2, 1) - ...
    nanmedian(s_data.velocity,1).^2);

% substep 2 --> define threshold
s_data.MSthresh = repmat(stdxy*6, size(s_data.velocity,1), 1, 1);

end

function s_data = local_select_MS(s_data, cfg)

% apply formula as in engbert 2006
numX = squeeze(s_data.velocity(:,1,:));
numY = squeeze(s_data.velocity(:,2,:));
denX = squeeze(s_data.MSthresh(:,1,:));
denY = squeeze(s_data.MSthresh(:,2,:));

s_data.lgcl_mask_MS = ((numX./denX).^2 + (numY./denY).^2) > 1;

% noise reduction: 
% the easiest way to implement a selection criterion is by convolving the
% logical with a kernel of ones. This allows you to have a fast way of
% defining MS vs noise
% open question: which kernel length to use? 3 datapoints -hence variable
% in time- or a fixed tw? 

lgcl_kernel = ones(1,cfg.kernellength);

% preallocate for speed
swap_lgcl_conv = nan(length(s_data.lgcl_mask_MS),...
    size(s_data.lgcl_mask_MS,2));

% start for loop in each trial. Ideally, even this loop could be avoided,
% maybe, by using 2dconv. Nevertheless I did not understand how it works,
% and I prefer to stick to the original
for iTrl = 1:size(s_data.lgcl_mask_MS,2)

    vectConv = conv(s_data.lgcl_mask_MS(:,iTrl), lgcl_kernel)...
        >=cfg.kernellength;
    
    % we have to account for the fact that any value reaching the threshold
    % shows the "spike" at a delay due to the convolution process. There is
    % the need to put ones even at the beginning of the MS itself, done by
    % find
    spikesP = find(vectConv); deltaTcorrector = -(0:cfg.kernellength-1);
    fooidx = repmat(spikesP, 1, numel(deltaTcorrector));
    footmd = repmat(deltaTcorrector, numel(spikesP), 1);
    rightIdx = unique(fooidx+footmd);
    vectConv(rightIdx) = true;
    % cut the last N-1 elements, were N is the kernel length
    vectConv((end-cfg.kernellength+2):end) = [];
    swap_lgcl_conv(:,iTrl) = vectConv;
    
end

% attach the matrix obtained in this way to the data
% s_data.old_lgcl_mask = s_data.lgcl_mask_MS;
s_data.lgcl_mask_MS = logical(swap_lgcl_conv);

end

function s_data = local_select_toi(s_data, cfg)

lgcl_mask_time = s_data.x_time>=min(cfg.toi) & s_data.x_time<=max(cfg.toi);

% now apply this to the subfields of interest
s_data.trial = s_data.trial(lgcl_mask_time,:,:);
s_data.velocity = s_data.velocity(lgcl_mask_time,:,:);
s_data.MSthresh = s_data.MSthresh(lgcl_mask_time,:,:);
s_data.lgcl_mask_MS = s_data.lgcl_mask_MS(lgcl_mask_time,:,:);
s_data.x_time = s_data.x_time(lgcl_mask_time);

end

function s_data = local_compute_angles(s_data, xchan_idx, ychan_idx)

ntrl = length(s_data.trialinfo);

s_data.MS_features = cell(1, ntrl);
s_data.lastMS = nan(ntrl, 6);
 
for iTrl = 1:ntrl
    
    % onset & offset on X axis
    x_onsets = s_data.trial(s_data.lgcl_MS_onset(:,iTrl), xchan_idx, iTrl);
    x_offsets = s_data.trial(s_data.lgcl_MS_offset(:,iTrl), xchan_idx, iTrl);
    
    % onset & offset on Y axis
    y_onsets = s_data.trial(s_data.lgcl_MS_onset(:,iTrl), ychan_idx, iTrl);
    y_offsets = s_data.trial(s_data.lgcl_MS_offset(:,iTrl), ychan_idx, iTrl);
    
    % difference vectors
    diff_vects = [x_offsets, y_offsets] - [x_onsets, y_onsets];
    
    % compute angles
    these_angles = atan2(diff_vects(:,2), diff_vects(:,1));
    
    % store time -saccade offset?
    these_Ts = s_data.x_time(s_data.lgcl_MS_offset(:,iTrl));
    
    % small summary matrix of MS direction and T. only last saccade will be
    % taken into account 
    % colord = 1) Time, 2) x onset, 3) y onset, 4) x offset, 5) y offset, 6) angles
    summat = [these_Ts, x_onsets, y_onsets, x_offsets, y_offsets,...
        these_angles]; 
    
    s_data.MS_features{iTrl} = summat;
    
    if ~isempty(summat)

        s_data.lastMS(iTrl,:) = summat(end,:);
    
    else
        
        s_data.lastMS(iTrl, :) = nan;
        
    end
    
end


end

