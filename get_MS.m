function [s_data] = get_MS(s_data, cfg)
% implement MS detection algorithm described in Engbert & Kliegl (2003)
% currently the only difference in th ealgorithm per se is that, since we
% are collecting data from just one eye, we'll skip the selection of
% microsaccadic event happening for both eyes.
% started by EB on 15-Apr-2019

% subfields necessary
% cfg.screensize = [1920 1080]; 
% cfg.pix_per_va = 51;

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
s_data = local_select_MS(s_data);

% step 4 compute angles
s_data.angle = atan(squeeze(s_data.velocity(:,2,:))./...
    squeeze(s_data.velocity(:,1,:)));


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

deltaT = s_data.x_time(2)-s_data.x_time(1);
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

function s_data = local_select_MS(s_data)

% apply formula as in engbert 2006
num1 = squeeze(s_data.velocity(:,1,:));
num2 = squeeze(s_data.velocity(:,2,:));
den1 = squeeze(s_data.MSthresh(:,1,:));
den2 = squeeze(s_data.MSthresh(:,2,:));

s_data.lgcl_mask_MS = ((num1./den1).^2 + (num2./den2).^2) > 1;

end

function s_data = local_visualize_selection(s_data)

figure
for iPlot=1:size(s_data.velocity,3)

    this_mat_plot = [squeeze(s_data.trial(:,[2 3],iPlot)),...
        s_data.lgcl_mask_MS(:,iPlot)];
    
    plot(this_mat_plot, 'LineWidth', 2)
    waitforbuttonpress
         
end


%% scatterplot single MSs

figure
for iPlot=1:size(s_data.velocity,3)

    velX = squeeze(s_data.velocity(:,1,iPlot));
    velY = squeeze(s_data.velocity(:,2,iPlot));

    MSs_ev_X = velX;
    MSs_ev_X(~s_data.lgcl_mask_MS(:,iPlot))= nan;
    MSs_ev_Y = velY;
    MSs_ev_Y(~s_data.lgcl_mask_MS(:,iPlot))= nan;

    plot(velX, velY, 'k'); hold on
    plot(MSs_ev_X, MSs_ev_Y, 'r', 'LineWidth', 4)
    hold off
    
    waitforbuttonpress
    
end

polarplot(s_data.angle(:,1), sqrt(squeeze(s_data.velocity(:,1,1)).^2 +...
    squeeze(s_data.velocity(:,2,1)).^2))

polarplot(s_data.angle(:,1), ones(3001,2))

end