function m_hat= decoder_1error(c,g)
n=31;
d=length(g)-1;
S1=zeros(31,d );

for i=1:n 
    input_bits=zeros(1,n);
    input_bits(i)=1;
    t=encoder(input_bits,g);
    S1(i,:)=t(1:d);
end
e=encoder(c,g);
e=e(1:d);
if 1 && all(e==0) 
    
    m_hat=c(d+1:n);
else
    
    err_position=find(ismember(e,S1),1); % retourne indice de la ligne correspondante a l'erreur e
    
    c(err_position)=c(err_position)+1;
    
    m_hat=c(d+1:n);

end