tau = [0,0.5,1,1.5,2];
A_1 = [1,0.1,0.1,0.1,0.1];
A_2 = [1,0.8,0.6,0.4,0.2];
A_3 = [1,0.8,0.8,0.8,0.8];
Ts = 0.05e-6;
L = 6;
m = -L:1:L;

h1 = filtre_canal(m,A_1,tau,Ts,L);
h1_norm = h1/norm(h1);

h2 = filtre_canal(m,A_2,tau,Ts,L);
h2_norm = h2/norm(h2);

h3 = filtre_canal(m,A_3,tau,Ts,L);
h3_norm = h3/norm(h3);

%Example of an ugly plot

plot(m,h1_norm);