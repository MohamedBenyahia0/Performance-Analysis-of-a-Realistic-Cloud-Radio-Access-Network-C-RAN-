function [H1,H2,H3] =plot_channel()
tau = [0,0.5,1,1.5,2];
A_1 = [1,0.1,0.1,0.1,0.1];
A_2 = [1,0.8,0.6,0.4,0.2];
A_3 = [1,0.8,0.8,0.8,0.8];
Ts = 0.05e-6;
L = 6;
m = -L:1:L;
N=100;
h1 = filtre_canal(m,A_1,tau,Ts,L);
h1_norm = h1/norm(h1);
H1=toeplitzMatrix(h1_norm,N,L);
h2 = filtre_canal(m,A_2,tau,Ts,L);
h2_norm = h2/norm(h2);
H2=toeplitzMatrix(h2_norm,N,L);
h3 = filtre_canal(m,A_3,tau,Ts,L);
h3_norm = h3/norm(h3);
H3=toeplitzMatrix(h3_norm,N,L);

plot(m,h1_norm,'-o',"LineWidth",3.0);
title('Channels L=6 N=100');
xlabel("k" );
ylabel("h(k)");

hold on;
plot(m,h2_norm,'-*',"LineWidth",3.0);

plot(m,h3_norm,'-x',"LineWidth",3.0);
legend('Channel 1','Channel 2','Channel 3');
hold off;
