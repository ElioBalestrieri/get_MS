function [MSstruct, s_data] = get_MS(EEG, cfg)
% implement MS detection algorithm described in Engbert & Kliegl (2003)
% with further modifications in Engbert and Mergenthaler *PNAS* (2006).
% function optimized to work on data coming from Buschlab
% even if the function is created for handling data from 2 eyes, it does
% not contain an algorithm of selection of MS computed from both eyes
% simultaneously. It is advisable to average the signal from both eyes
% before submitting it to the function (as in Lin, Nobre and Van Ede, 
% Nature Comms 2022).
%
% started by EB on 15-Apr-2019
% modified on 25-08-2019 to allow double channel recording
% modified on 02-05-2022 to adapt to standard buschlab output
%                        to allow easy MS labeling
%
% MANDATORY
% cfg.kernellength = number of contiguous datapoints necessary for defining
%                    a MS. 2-3;
% cfg.lambda = denominator in velocity formula. Suggested value:5 
%              after Engbert and Mergenthaler 2006 PNAS;
%
% OPTIONAL:
% cfg.labelMS = {[-cos, cos], [-cos, cos]}; automatically label MS
%               direction based on the cosine of the MS angle. Cosine was
%               chosen because it makes easy to label "left" vs "right" MS.
%               Please note that this measure is not suitable for "high" vs
%               "low" MS!!! When this field is present, the resulting
%               structure will contain a matrix Ntimepoint X Mtrials where
%               0 indicates absence of MS, 1 MS to 1st direction defined, 2
%               MS to 2nd direction defined and so on...
%
% cfg.toi = select time of interest for the output based on times
%
% cfg.noisesuppress = apply signal smoothing (default = true);
%
% cfg.smoothkernellength = length of the gaussian window for signal
%                          smoothing (default = 5);


%% convert input

% find eye(s) position channels
chanlabels = {'Eyegaze-X', 'Eyegaze-Y'};
tbl_chan = struct2table(EEG.chanlocs);
xy_idxs = [find(ismember(tbl_chan.labels, chanlabels{1})), ...
           find(ismember(tbl_chan.labels, chanlabels{2}))];
       
if length(xy_idxs) < 2
    error('not enough eye channels found')
end

s_data = [];
s_data.trial = permute(EEG.data(xy_idxs, :, :), [2, 1, 3]);
s_data.label = chanlabels;
s_data.x_time = EEG.times';

if isfield(EEG, 'trialinfo')
    s_data.trialinfo = EEG.trialinfo;
end
%% extract MS
lambda = cfg.lambda;

xchan_idx = 1;
ychan_idx = 2;

% define defaults for smoothing
if ~isfield(cfg, 'noisesuppress')
    cfg.noisesuppress = true;
end

if ~isfield(cfg, 'smoothkernellength')
    cfg.smoothkernellength = 5;
end

if ~isfield(cfg, 'detrend')
    cfg.detrend = false;
end


% step 0: smooth (if asked)
if cfg.noisesuppress
    s_data = local_smooth_data(s_data, cfg);    
end

if cfg.detrend
    s_data = local_detrend_data(s_data);    
end


% step1 compute velocity
s_data = local_compute_velocity(s_data, xchan_idx, ychan_idx);

% step2 compute threshold
s_data = local_compute_threshold(s_data, lambda);

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
s_data.lgcl_MS_onset = logical(cat(1, s_data.lgcl_mask_MS(1,:,:),...
    swap_single_saccade_ON));

% it has to be noted that the single MS onset might be werid and not than
% informative. Solution: take even the MS offset, and compute the angle
% between the onset and offset of MS, on the data and not on the speed.
swap_single_saccade_OFF = diff(s_data.lgcl_mask_MS) == -1;
s_data.lgcl_MS_offset = logical(cat(1, swap_single_saccade_OFF, ...
    s_data.lgcl_mask_MS(end,:,:))); % instead of assigning 0s, let's assign the
                                    % last row of the MS, in order to correct
                                    % for MSs started but not conluded on
                                    % time

% ... and only now we compute angles
s_data = local_compute_angles(s_data, xchan_idx, ychan_idx);

% finally, assign directional labels to MS (if required)
if isfield(cfg, 'labelMS')
    
    s_data = local_labelMS(s_data, cfg);
    
end

% provide a reduced output to impreve readability
% (but maintain the chance to full s_data output)
MSstruct = [];
MSstruct.lgcl_mask_MS = s_data.lgcl_mask_MS;
MSstruct.labeled_MS = s_data.labeled_MS;
MSstruct.angle_maskMS = s_data.angle_maskMS;
MSstruct.timerange = minmax(s_data.x_time');
MSstruct.trialinfo = s_data.trialinfo;

end

%% ######################### LOCAL FUNCTIONS ##############################
function s_data = local_smooth_data(s_data, cfg)

dat = s_data.trial;

nchans = size(dat, 2); ntp = size(dat, 1); ntrls = size(dat, 3);
nzeros = cfg.smoothkernellength;

zeropad = zeros(nzeros, nchans, ntrls);
dat = cat(1, zeropad, dat, zeropad);

gkrnl = gausswin(cfg.smoothkernellength); den = sum(gkrnl);

for itrl = 1:ntrls
    
    for ichan = 1:nchans
        
        vect = dat(:, ichan, itrl);
        
        dat(:,ichan,itrl) = conv(vect, gkrnl, 'same')./den;

    end
    
end

s_data.trial = dat(nzeros+1:end-nzeros, :, :);

end

function s_data = local_detrend_data(s_data)

nchans = size(s_data.trial, 2);

for ichan = 1:nchans
    
    mat_ = squeeze(s_data.trial(:, ichan, :));
    detmat_ = detrend(mat_);
    s_data.trial(:, ichan, :) = detmat_;

end

end

function s_data = local_compute_velocity(s_data, xchan_idx, ychan_idx)

totchan = [xchan_idx, ychan_idx]; nchan = numel(totchan);

% create matrices of values for velocity computation
mat_nminus1 = cat(1, nan(1,nchan,size(s_data.trial, 3)), ...
    s_data.trial(1:end-1, [xchan_idx, ychan_idx], :));
mat_nminus2 = cat(1, nan(2,nchan,size(s_data.trial, 3)), ...
    s_data.trial(1:end-2, [xchan_idx, ychan_idx], :));
mat_nplus1 = cat(1, s_data.trial(2:end, [xchan_idx, ychan_idx], :),...
    nan(1,nchan,size(s_data.trial, 3)));
mat_nplus2 = cat(1, s_data.trial(3:end, [xchan_idx, ychan_idx], :),...
    nan(2,nchan,size(s_data.trial, 3)));

deltaT = mean(diff(s_data.x_time));
s_data.velocity = (mat_nplus2 + mat_nplus1 - mat_nminus1 - mat_nminus2)/...
    (6*deltaT);

end

function s_data = local_compute_threshold(s_data, lambda)

% substep 1 --> compute std by applying a median estimator to the TS
% big question: mismatch between 2003 (no square root) and 2006 (sqrt)
% articles. makes more sense to me to have root squared vals though
stdxy = sqrt(nanmedian(s_data.velocity.^2, 1) - ...
    nanmedian(s_data.velocity,1).^2);

% substep 2 --> define threshold
s_data.MSthresh = repmat(stdxy*lambda, size(s_data.velocity,1), 1, 1);

end

function s_data = local_select_MS(s_data, cfg)

% apply formula as in engbert 2006
numX = squeeze(s_data.velocity(:,1,:));
numY = squeeze(s_data.velocity(:,2,:));
denX = squeeze(s_data.MSthresh(:,1,:));
denY = squeeze(s_data.MSthresh(:,2,:));

% now compare info from 2 eyes (if present)
if size(s_data.MSthresh, 2)>2
    
    % still compute the previous mask
    mask_eye_1 = ((numX./denX).^2 + (numY./denY).^2) > 1;
    
    numX2 = squeeze(s_data.velocity(:,3,:));
    numY2 = squeeze(s_data.velocity(:,4,:));
    denX2 = squeeze(s_data.MSthresh(:,3,:));
    denY2 = squeeze(s_data.MSthresh(:,4,:));

    mask_eye_2 = ((numX2./denX2).^2 + (numY2./denY2).^2) > 1;  
    
    s_data.lgcl_mask_MS = cat(4,mask_eye_1, mask_eye_2);

else
    
    s_data.lgcl_mask_MS = ((numX./denX).^2 + (numY./denY).^2) > 1;

end

% the easiest way to implement a selection criterion is by convolving the
% logical with a kernel of ones. This allows you to have a fast way of
% defining MS vs noise

lgcl_kernel = ones(1,cfg.kernellength);

% preallocate for speed
swap_lgcl_conv = nan(size(s_data.lgcl_mask_MS,1),...
    size(s_data.lgcl_mask_MS,2),...
    size(s_data.lgcl_mask_MS,4));

% start for loop in each trial. Ideally, even this loop could be avoided,
% maybe, by using 2dconv. Nevertheless I did not understand how it works,
% and I prefer to stick to the original

% update: keep separated matrices for the two eyes
for iEye = 1:size(s_data.lgcl_mask_MS,4)
    for iTrl = 1:size(s_data.lgcl_mask_MS,2)

        vectConv = conv(squeeze(s_data.lgcl_mask_MS(:,iTrl, iEye)), lgcl_kernel)...
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
        swap_lgcl_conv(:,iTrl,iEye) = vectConv;

    end
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

% number of eyetracked is given by the number of x channels
neye = numel(xchan_idx);

s_data.MS_features = cell(neye, ntrl);
s_data.lastMS = nan(ntrl, 6, neye);
s_data.avgMS = nan(ntrl, 3, neye);
[s_data.T_onoff, s_data.MS_angles] = deal(cell(ntrl, neye));
s_data.angle_maskMS = nan(length(s_data.x_time), ntrl);

for iEye = 1:neye

    for iTrl = 1:ntrl

        % onset & offset on X axis
        x_onsets = s_data.trial(s_data.lgcl_MS_onset(:,iTrl, iEye), xchan_idx(iEye), iTrl);
        x_offsets = s_data.trial(s_data.lgcl_MS_offset(:,iTrl, iEye), xchan_idx(iEye), iTrl);

        % onset & offset on Y axis
        y_onsets = s_data.trial(s_data.lgcl_MS_onset(:,iTrl, iEye), ychan_idx(iEye), iTrl);
        y_offsets = s_data.trial(s_data.lgcl_MS_offset(:,iTrl, iEye), ychan_idx(iEye), iTrl);

        % difference vectors
        diff_vects = [x_offsets, y_offsets] - [x_onsets, y_onsets];

        % compute "Average MS" as the sum of MS vectors, and compute angles.
        % Then store it in a corresponding subfield in s_data.
        res_vect = sum(diff_vects,1);
        if ~isempty(res_vect)     
            angle_avg = atan2(res_vect(2), res_vect(1)); % Y before, x after
            % colord: 1) angle, 2) X, 3) Y
            s_data.avgMS(iTrl, :, iEye) = [angle_avg, res_vect];
        end

        % compute angles
        these_angles = atan2(diff_vects(:,2), diff_vects(:,1));
        
        % create another matrix where MS are defined by their angle and not
        % by simple ones
        idxs_onset = find(s_data.lgcl_MS_onset(:, iTrl, iEye));
        idxs_offset = find(s_data.lgcl_MS_offset(:, iTrl, iEye));
        
        acc = 0;
        for thisangle = these_angles'        
            acc = acc+1;
            s_data.angle_maskMS(idxs_onset(acc):idxs_offset(acc), iTrl) = thisangle;
        end

        % store time -saccade offset?
        these_Ts = s_data.x_time(s_data.lgcl_MS_offset(:,iTrl, iEye));
        
        % create another subfield specifically for time onsets and offsets
        % for each mMS
        onset_offset = [s_data.x_time(s_data.lgcl_MS_onset(:,iTrl, iEye)), ...
            s_data.x_time(s_data.lgcl_MS_offset(:,iTrl, iEye))];

        % small summary matrix of MS direction and T. only last saccade will be
        % taken into account 
        % colord = 1) Time, 2) x onset, 3) y onset, 4) x offset, 5) y offset, 6) angles
        
        summat = [these_Ts, x_onsets, y_onsets, x_offsets, y_offsets,...
            these_angles]; 

        % attech to the structure
        s_data.MS_features{iEye, iTrl} = summat;
        s_data.T_onoff{iTrl, iEye} = onset_offset;
        s_data.MS_angles{iTrl, iEye} = these_angles;

        
        
        if ~isempty(summat)

            s_data.lastMS(iTrl,:, iEye) = summat(end,:);

        else

            s_data.lastMS(iTrl, :, iEye) = nan;

        end

    end

end
end

function s_data = local_labelMS(s_data, cfg)

acc = 0;
labeled_ms = zeros(size(s_data.lgcl_mask_MS));
xax_angle = cos(s_data.angle_maskMS);

for iboundary = cfg.labelMS
    
    acc = acc+1;
    lowbd = min(iboundary{1}); upbd = max(iboundary{1}); 
    thismask = ((xax_angle>=lowbd) & (xax_angle<=upbd))*acc;
    
    labeled_ms = labeled_ms + thismask;
    
end

s_data.labeled_MS = labeled_ms;

% store separately the directional MSs as logicals
acc = 0; cell_labeled_MS = cell(length(cfg.labelMS), 1);
for iboundary = cfg.labelMS
    
    acc = acc+1;
    cell_labeled_MS{acc} = s_data.labeled_MS == acc;
    
end

s_data.cell_labeled_MS = cell_labeled_MS;

end