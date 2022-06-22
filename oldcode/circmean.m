function m_angle = circmean(theta)
% compute mean of circular quantities

m_angle = atan2(nanmean(sin(theta)), nanmean(cos(theta)));


end