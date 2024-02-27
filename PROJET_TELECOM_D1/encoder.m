function c = encoder(m, g)
    new_val = zeros(1, length(g)-1);
    for i = length(m):-1:1
        old_val = new_val;
        new_val(1) = bitxor(m(i), old_val(end));
        for j = 2:length(new_val)
            if (g(j) == 1)
                new_val(j) = bitxor(old_val(j-1), old_val(end));
            else
                new_val(j) = old_val(j-1);
            end
        end
        
    end
    m_c=new_val;

    c=cat(2,m_c,m);% On concatene m_c et m 
    

    
