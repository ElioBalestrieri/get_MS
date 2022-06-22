function plot_summary(s_data)

figure; 
subplot(2, 3, 1:2); hold on

sum_mask = sum(s_data.lgcl_mask_MS, 2);
plot(s_data.x_time, sum_mask)

if isfield(s_data, 'cell_labeled_MS')
    
    for imat = s_data.cell_labeled_MS'

        sum_mask_label = sum(imat{1}, 2);
        plot(s_data.x_time, sum_mask_label)    

    end

end    
    
xlabel(minmax(s_data.x_time));
ylabel('N microsaccades detected')
xlabel('time')

subplot(2, 3, 3)
vectangle = cat(1, s_data.MS_angles{:});
polarhistogram(vectangle)

subplot(2, 3, 4:5)

% MS rate
cfg_MSrate = [];
cfg_MSrate.fsample = 256;
cfg_MSrate.movingwin = .1;

MSrate = get_MSrate_EXPcondition(s_data, true(length(s_data.trialinfo), 1), cfg_MSrate);

plot(s_data.x_time, MSrate.allMS)        
xlabel(minmax(s_data.x_time));
ylabel('MS rate /s')
xlabel('time')

subplot(2, 3, 6)

% single trial
for itrl = 1:size(s_data.trial, 3)
    
    thistrl = s_data.trial(:, :, itrl);
    plot(thistrl(:, 1),thistrl(:, 2), 'k'); hold on
    
    theseMS = s_data.lgcl_mask_MS(:, itrl);
    only_ms = thistrl;
    only_ms(~theseMS, :) = nan;
    
    plot(only_ms(:, 1), only_ms(:, 2), 'r', 'LineWidth', 2)
    
    theseMSoffsets = find(s_data.lgcl_MS_offset(:, itrl))';
    theseAngles = s_data.MS_angles{itrl};
    
    acc = 0;
    
    for idx = theseMSoffsets
        
        acc = acc+1;
        
        xoffs = double(thistrl(idx, 1));
        yoffs = double(thistrl(idx, 2));

        text(xoffs, yoffs, num2str(theseAngles(acc)), 'Color', [1, 0, 0])
        
    end
    
    hold off
    
    waitforbuttonpress
    
    
end




end