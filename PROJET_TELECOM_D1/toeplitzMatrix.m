function H = toeplitzMatrix(h,N,L)
H=zeros(N,N);
for i=1:N
    for k=1:N
        if L+1+k-i>=1 && L+1+k-i<=2*L+1
            H(i,k)=h(L+1+k-i);
        end
    end
  

end
end