function N=taylormatrixR(S,pKa,pH,T,I_m,constant)

[m,n]=size(S);
R = 8.31e-3; % kJ/mol/K
no_pseudo=[];
for i=1:m
   no_pseudo=[no_pseudo;length(diag(pKa(i).pKas, 1)')];
end
pKa_M=[];
for i=1:m
   pKa_M=[pKa_M;(diag(pKa(i).pKas, 1))];
end
t=length(pKa_M);
N=zeros(n,1+t+n+m);
pka_m=zeros(m ,max(no_pseudo));
for i=1:size(pka_m,1)
        kegg=diag(pKa(i).pKas,1);
    for j=1:no_pseudo(i)
        pka_m(i,j)=kegg(j);
    end
end
pka_m(isnan(pka_m)) = 0 ; 
pka=sym('pka',[m max(no_pseudo)]);
I=sym('I',[n 1]);
I_m(isnan(I_m)) = 0.25 ; 
for i=1:n
    inds=find(S(:,i));
    B=sym('B',[length(inds) 1]);
    for j=1:length(inds)
        diss = pKa(inds(j));
        T=T(i);
        pH=pH(i);
        dG0s = cumsum(-[0, pka(inds(j),1:no_pseudo(inds(j)))] * R * T * log(10));
        dG0s = dG0s - dG0s(diss.majorMSpH7);
        pseudoisomers = [dG0s(:), diss.nHs(:), diss.zs(:)];
        alpha = (9.20483*T)/10^3 - (1.284668*T^2)/10^5 + (4.95199*T^3)/10^8; % Approximation of the temperature dependency of ionic strength effects
        if I_m(i)>0
            DH = (alpha * sqrt(I(i))) / (1 + 1.6 * sqrt(I(i))); % Debye Huckel
        else
            DH=0;
        end
        dG0_prime_vector = pseudoisomers(:, 1) + ...
                           pseudoisomers(:, 2) * (R*T*log(10)*pH + DH) - ...
                           pseudoisomers(:, 3).^2 * DH;
        B(j)=-R*T*log(sum(exp(dG0_prime_vector / (-R * T))));
    end
        f=S(inds,i)'*B;
        if ~isequal(f,sym('0'))
            g=subs(f,pka,pka_m);
            N(i,1)=vpa(subs(g,I,I_m));
            if ~constant %For linear data reconciliation, only constant terms are required
                    F=symvar(f);
                    h=jacobian(f,F);
                      if ~isempty(F)
                          h_1=subs(h,pka,pka_m);
                          diffs=vpa(subs(h_1,I,I_m));
                          N(i,1+t+i)=diffs(1);
                             for k=2:length(F)
                                 str=char(F(k));
                                 index_c=str2double(str(strfind(str,'a')+1:strfind(str,'_')-1));
                                 index_pseudo=str2double(str((strfind(str,'_')+1):end));
                                 L=sum(no_pseudo(1:index_c-1));
                                 N(i,L+index_pseudo)=diffs(k);
                             end
                      end
            end          

        end 
        
end

end


            

            
            
        
        
        
        
        
        
        
        
                       
               
               
        
        
        
        
    
    