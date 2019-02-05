function J=costFunction(x,b,Gf_m,G_m,sigma)
sigma=diag(sigma);
sigma=[sigma;ones(length(Gf_m),1)];
sigma=diag(sigma);
A=(1:length(b))';
B=((length(b)+1):(length(b)+length(Gf_m)))';
C=((length(b)+length(Gf_m)+1):length(x))';
temp=[(b-x(A));(Gf_m-x(B));(x(B)-((G_m)*x(C)))];

J=(1/(2*size(temp,1)))*sum(temp'*sigma*temp);


