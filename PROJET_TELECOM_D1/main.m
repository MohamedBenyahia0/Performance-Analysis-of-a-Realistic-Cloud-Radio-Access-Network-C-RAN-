g1=[1,0,1,0,0,1];%g=x^5+x^2+1 de degre n-k=5
k1=26;%k=31-5
m1=ones(1,k1);
c1=encoder(m1,g1);

c1_prim=c1;
c1_prim(4)=c1_prim(4)+1;
m1_hat=decoder_1error(c1_prim,g1);

%%%%%%%%%%%%%
g2=[1,0,0,1,0,1,1,0,1,1,1];%g=x^10+x^9+x^8+x^6+x^5+x^3+1
k2=31-10;
m2=ones(1,k2);
c2=encoder(m2,g2);

c2_prim=c2;
c2_prim(3)=c2_prim(3)+1;
c2_prim(6)=c2_prim(6)+1;
m2_hat=decoder_2errors(c2_prim,g2);


