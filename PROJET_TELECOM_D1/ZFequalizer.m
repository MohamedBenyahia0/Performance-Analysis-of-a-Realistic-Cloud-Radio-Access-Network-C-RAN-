function [s_est,loc]=ZFequalizer(y,H,mod,M)
z=pinv(H)*y;
[s_est,loc] = threshold_detector(z,mod,M);

end
