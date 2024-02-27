function [s_est,loc]=ZFequalizer(y,H,mod,M)
z=pinv(H)*y;
[n,m ]=size(z);
z=reshape(z,[n*m,1]);
[s_est,loc] = threshold_detector(z,mod,M);

end
