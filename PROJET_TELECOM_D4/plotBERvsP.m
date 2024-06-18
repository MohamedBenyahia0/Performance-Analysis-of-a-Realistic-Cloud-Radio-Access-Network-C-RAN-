
MC=50;%number of Monte-Carlo simulations
M=2;%% taille de la constellation
mod = 'PSK';
Ts=1/(20e6);
lamda = 3e8/2.4e9;
No = -174;

N_linear =1e-3*(10^(No/10))*B;
N_dB = 10*log10(N_linear);
c=3e8;
f0=2.4e9;
num_user = 3;

P=0:2:20; %% Eb/No en dBm
nb_err = zeros(1,length(P));
empiricalBER_without_interf=zeros(1,length(P));
Energy_symbol = 1;

N_bits= 1e3; % nombre de bits d'information



parfor ii = 1:length(P)
    
    nb_err_sim=0;
    for mm=1:MC
        bits = randi([0,1],num_user,N_bits); % generation des bits d'information
        symboles_BPSK = zeros(num_user,N_bits);
    
        
        P_lin = 1e-3*10^(P(ii)/10);
        for i =1:num_user
            symboles_BPSK(i,:) = bits2symbols(bits(i,:),mod,M);

        end
        P_max=1e-3*10^(20/10);
        symboles_BPSK(1,:)=symboles_BPSK(1,:).*sqrt(P_max);
        symboles_BPSK(3,:)=symboles_BPSK(3,:).*sqrt(P_max);
        symboles_BPSK(2,:)=symboles_BPSK(2,:).*sqrt(P_lin);
        
        
        
        
        % Centers of the three cells
        centers =[-2 0
                   0 0
                   2 0];
    
        % Preallocate arrays for user coordinates
        x_coords = zeros(1, 3);
        y_coords = zeros(1, 3);
        h = zeros(1,3);
    
        % Generate uniformly distributed points within each cell
        for i = 1:num_user
            center_x = centers(i, 1);
            center_y = centers(i, 2);
    
            % Generate uniformly distributed points in the range [-semi_length, semi_length]
            x_coords(i) = center_x + (-1 + 2 * 1 * rand(1, 1));
            y_coords(i) = center_y + (-1 + 2 * 1 * rand(1, 1));
    
            d2=x_coords(i).^2+y_coords(i).^2;%% square-distance from the origin
            h(i)=sqrt(min(1, (c/(4*pi*f0))^2./(d2*1e6)));% square magnitude attenuation in Friis equation 
        end
    
        
        W2 =sqrt(N_linear/2)* (randn(1,N_bits) + 1i*randn(1,N_bits));  % noise
        H=[0 h(2) 0];
      
        
    
        Z= H*symboles_BPSK+ W2;  % signal recu après FFT
        z_vec= reshape(Z,[N_bits,1]);
        z_vec=z_vec./(h(2)*sqrt(P_lin));
        me=mean(z_vec);
        z_vec=(z_vec-mean(z_vec))/std(z_vec);
    
        
        s_est = threshold_detector(z_vec, mod, M) ; % prise de décision
        bits_est = symbols2bits(s_est,mod,M); % convertir les symboles en bits
    
        
        nb_err_sim=nb_err_sim+sum(abs(bits_est-bits(2,:))); % compter le nombre d'erreur
        
    end
    nb_err(ii)=nb_err_sim/MC;
    empiricalBER_without_interf(ii) = nb_err(ii)/(N_bits);  % calcul du taux d'erreur binaire

end
%% with interference
P=0:2:20; %% Eb/No en dBm
nb_err = zeros(1,length(P));
empiricalBER_with_interf=zeros(1,length(P));
Energy_symbol = 1;

N_bits= 1e3; % nombre de bits d'information



parfor ii = 1:length(P)
    
    nb_err_sim=0;
    for mm=1:MC
        bits = randi([0,1],num_user,N_bits); % generation des bits d'information
        symboles_BPSK = zeros(num_user,N_bits);
    
        
        P_lin = 1e-3*10^(P(ii)/10);
        for i =1:num_user
            symboles_BPSK(i,:) = bits2symbols(bits(i,:),mod,M);

        end
        P_max=1e-3*10^(20/10);
        symboles_BPSK(1,:)=symboles_BPSK(1,:).*sqrt(P_max);
        symboles_BPSK(3,:)=symboles_BPSK(3,:).*sqrt(P_max);
        symboles_BPSK(2,:)=symboles_BPSK(2,:).*sqrt(P_lin);
        
        
        
        
        % Centers of the three cells
        centers =[-2 0
                   0 0
                   2 0];
    
        % Preallocate arrays for user coordinates
        x_coords = zeros(1, 3);
        y_coords = zeros(1, 3);
        h = zeros(1,3);
    
        % Generate uniformly distributed points within each cell
        for i = 1:num_user
            center_x = centers(i, 1);
            center_y = centers(i, 2);
    
            % Generate uniformly distributed points in the range [-semi_length, semi_length]
            x_coords(i) = center_x + (-1 + 2 * 1 * rand(1, 1));
            y_coords(i) = center_y + (-1 + 2 * 1 * rand(1, 1));
    
            d2=x_coords(i).^2+y_coords(i).^2;%% square-distance from the origin
            h(i)=sqrt(min(1, (c/(4*pi*f0))^2./(d2*1e6)));% square magnitude attenuation in Friis equation 
        end
    
        
        W2 =sqrt(N_linear/2)* (randn(1,N_bits) + 1i*randn(1,N_bits));  % noise
        H=[h(1) h(2) h(3)];
      
        
    
        Z= H*symboles_BPSK+ W2;  % signal recu après FFT
        z_vec= reshape(Z,[N_bits,1]);
        z_vec=z_vec./(h(2)*sqrt(P_lin));
        me=mean(z_vec);
        z_vec=(z_vec-mean(z_vec))/std(z_vec);
    
        
        s_est = threshold_detector(z_vec, mod, M) ; % prise de décision
        bits_est = symbols2bits(s_est,mod,M); % convertir les symboles en bits
    
        
        nb_err_sim=nb_err_sim+sum(abs(bits_est-bits(2,:))); % compter le nombre d'erreur
        
    end
    nb_err(ii)=nb_err_sim/MC;
    empiricalBER_with_interf(ii) = nb_err(ii)/(N_bits);  % calcul du taux d'erreur binaire

end

%%
figure;
semilogy(P,empiricalBER_without_interf,'-x','LineWidth',2.0);
hold on;
semilogy(P,empiricalBER_with_interf,'-x','LineWidth',2.0)
grid on
xlabel('P (in dBm)');ylabel('BER'); title('BER vs P');
legend('without interference','with interference')
grid on;
hold off;