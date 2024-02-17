function m_hat = decoder_2errors(c,g)
n=31;
d=length(g)-1;
nb_1err=n;
nb_2err=nchoosek(n,2);

S2=zeros(nb_1err+nb_2err,d);
for i=1:n 
    input_bits=zeros(1,n);
    input_bits(i)=1;
    t=encoder(input_bits,g);
    S2(i,:)=t(1:d);
end
k=nb_1err+1;
for i=1:n-1
    for j=i+1:n
        input_bits=zeros(1,n);
        input_bits(i)=1;
        input_bits(j)=1;
        t=encoder(input_bits,g);
        S2(k,:)=t(1:d);
        k=k+1;
    end
end
e=encoder(c,g);
if 1 && all(e==0)
    m_hat=c(d+1:n);
else 
    s=find(ismember(e,S2),1); % retourne indice de la ligne correspondante a l'erreur e
    if s<=n  % 1 erreur
        c(s)=c(s)+1;
        m_hat=c(d+1:n);
    else % 2 erreur
        for i=1:n-1
            for j=i+1:n
                if S2(i,:)+S2(j,:)==S2(s,:)
                    c(i)=c(i)+1;
                    c(j)=c(j)+1;
                    m_hat=c(d+1,n);
                end
            end
        end
    end


end