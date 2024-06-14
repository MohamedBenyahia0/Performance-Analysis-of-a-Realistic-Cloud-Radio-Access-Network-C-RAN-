[K,rej_zf] = d4_perfs_students('40k','zf');
[K,rej_dfe] = d4_perfs_students('40k','dfe');
figure();
plot(K,rej_zf,"-o",'LineWidth',3.0);
hold on;
plot(K,rej_dfe,"-*",'LineWidth',3.0);
xlabel('Number of users')
ylabel('Rejection Rate')
title('Rejection rate versus number of users 40k'); 
legend('zf','dfe')
grid on;
hold off;
%%
[K,rej_zf] = d4_perfs_students('400k','zf');
[K,rej_dfe] = d4_perfs_students('400k','dfe');
figure();
plot(K,rej_zf,"-o",'LineWidth',3.0);
hold on;
plot(K,rej_dfe,"-*",'LineWidth',3.0);
xlabel('Number of users')
ylabel('Rejection Rate')
title('Rejection rate versus number of users 400k'); 
legend('zf','dfe')
grid on;
hold off;
%%
[K,rej_zf] = d4_perfs_students('4M','zf');
[K,rej_dfe] = d4_perfs_students('4M','dfe');
figure();
plot(K,rej_zf,"-o",'LineWidth',3.0);
hold on;
plot(K,rej_dfe,"-*",'LineWidth',3.0);
xlabel('Number of users')
ylabel('Rejection Rate')
title('Rejection rate versus number of users 4M'); 
legend('zf','dfe')
grid on;
hold off;
%%
[K,rej_zf] = d4_perfs_students('40M','zf');
[K,rej_dfe] = d4_perfs_students('40M','dfe');
figure();
plot(K,rej_zf,"-o",'LineWidth',3.0);
hold on;
plot(K,rej_dfe,"-*",'LineWidth',3.0);
xlabel('Number of users')
ylabel('Rejection Rate')
title('Rejection rate versus number of users 40M'); 
legend('zf','dfe')
grid on;
hold off;
%%
% decode and forward the basestation 80Mbps
%quantize and forward frequence echantillonnage * nombre de bits * 2 voies
%pour 16 QAm = 30*10*2=600Mbps

