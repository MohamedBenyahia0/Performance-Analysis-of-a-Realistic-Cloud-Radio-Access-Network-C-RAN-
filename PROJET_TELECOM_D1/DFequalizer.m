function [s_est,loc]=DFequalizer(y,H,mod,M)
[Q,R ]=qr(H);
z=Q'*y; % Q' hermitian of Q
N=length(z);

s_est=zeros(N,1); %s_est=[ s_est(N)...s_est(1)]^T
s_est(N)=threshold_detector(z(N)/R(N,N),mod,M);
for i=1:N-1
    z_prim_i=z(i);
    for k=i+1:N
        z_prim_i=z_prim_i-R(i,k)*s_est(k);
    end
    s_est(i)=threshold_detector(z_prim_i/R(i,i),mod,M);


end
