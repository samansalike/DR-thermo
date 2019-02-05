function [Recon_var,index_unobservable]=recon_l(S,y,index_rm,index_fm,sigma,use_gc)
                % linear data recociliation  
                % reconciles all available measurements and estimates observable missing data
                % inputs : G_r is a vector consisting of all measured reaction Gibbs energies
                %          G_f is a vector consisting of all the measured Gibbs energies of formation
                %          sigma is the inverse of the variance matrix
                %          index_m is a vector of indices of the rows of S corresponding to the available measurements
                %          ( here only missing formation energies are considered )
                % output : x= [x_m x_e]', where  
                %          x_m is the vector representing the reconciled estimates of all measured variables 
                %          x_e is the vector consisting of the estimates of missing formation energies that are observable
                %          index_observable are the indices corresponding to rows of S that represent the elements of x_e

if nargin<6
    use_gc=false;
     % use_gc denotes the use of group contribution estimates for unobservable variables
     % if use_gc is false, only the estimates of the observable variables are calculated
end

[m,n]=size(S);
index_m=[index_rm;n+index_fm];
index_fest=1:m;
index_fest(index_fm)=[]; 
index_rest=1:n;
index_rest(index_rm)=[]; 
index_est=[index_rest';n+index_fest'];  % indices of variables to be estimated

%constraint matrix
A=[eye(n) -S']; %model: [I -S']=[dG0r dG0f]';
%Splitting the constraint matrix according to measured and unmeasured Gibbs
%energies
Ax=A(:,index_m);
Au=A(:,index_est);
[Q,R,E]=qr(Au);  %QR decomposition of Au
s=rank(full(Au));%number of independant coloumns of Au , these variables are observable 
Q1=Q(:,(1:s));       
Q2=Q(:,(s+1:end));
R1=R((1:s),(1:s));
R2=R((1:s),(s+1:end));
P=Q2'; 
PIu=inv(E); % permutation matrix
index_nums=PIu*(1:length(index_est))'; % rearranged indices
index_obs=index_nums(1:s); % these variables can be uniquely estimated 
index_unobs=index_nums(s+1:end); % estimates for these variables are required to be specified
Ru=inv(R1)*R2; 
%% Classification of variables using Ru
        I=[];
        J=[];
        Ru=round(Ru,5);
        for i=1:size(Ru,1)
           if isequal(Ru(i,:),zeros(1,size(Ru,2))) % Only the corresponding variables are observable
              I=[I;i];
           end
           if ~isequal(Ru(i,:),zeros(1,size(Ru,2))) % The corresponding variables are also unobservable
              J=[J;i];
           end
        end
index_observable=index_est(sort(index_obs(I)));
index_unobservable=index_est(sort(index_obs(J)));
%%
x_m=y-((sigma*transpose(P*Ax)*inv((P*Ax)*sigma*transpose(P*Ax))*(P*Ax))*y); % reconciled measurements

if isempty(index_est)
    Recon_var=x_m;
    index_unobservable=[];
else
  if s==min(size(Au)) % if all unknown variables are observable
   x_e=-E*inv(R1)*(transpose(Q1)*Ax*x_m); %reconciled estimates of unmeasured variables
   Recon_var(index_m)=x_m;
   Recon_var(index_est)=x_e;
  else
    if use_gc
        %
         group_cache_file = '../cache/G.mat';
         if ~exist(group_cache_file)
            fprintf('Group incidence matrix not found in cache \n')
         end
         load(group_cache_file);
        %
       [G_grp,I,index_fout_range]=Groupcontributions(P*Ax,y(1:length(index_rm)),y(length(index_rm)+1:end),sigma,G,index_fest);
       G_f_gc=G(:,I)*G_grp; %group contribution estimates of formation Gibbs energies
       %Here G(:,I) is used because not all group contribution estimates
       %may be available. When the formation energy of a compound is
       %dependant on an unavialable group contibution, its corresponding
       %values in G_f_gc is incorrect and is thus made NaN subsequently
       G_r_gc=S'*G_f_gc; %group contribution estimates of reaction Gibbs energies
       G_f_gc(index_fout_range)=0; %these formation energies cannot be estimated using available group contributions
       index_rout_range=[];
       for i=1:n
           if intersect(find(S(:,i)),index_fout_range) 
               %i.e if any the reaction Gibbs energies 
               %depend on uncovered formation energies
              index_rout_range=[index_rout_range;i];
           end
       end
       G_r_gc(index_rout_range)=0;
       G_gc=[G_r_gc;G_f_gc];
       u_ps=zeros(length(index_unobs),1);
       for i=1:size(index_unobs,1)
           u_ps(i)=G_gc(index_est(index_unobs(i))); % u_ps is a vector consisting of the p-s variables that are unobservable
       end
       index_estimate=index_est;
       index_out_range=[index_rout_range;n+index_fout_range];
       index_GC_out=intersect(index_est(index_unobs),index_out_range);
       [~,IA]=intersect(index_estimate,index_GC_out);
       index_estimate(IA)=[];
       
    else 
        index_estimate=index_observable;
        u_ps=zeros(length(index_unobs),1);
    end
    u_s=-inv(R1)*((transpose(Q1)*Ax*x_m)+R2*u_ps); 
    x_e=zeros(length(index_est),1); % vector to store all estimated values
    u_jumbled=[u_s;u_ps];
    for i=1:length(index_nums)
        x_e(index_nums(i))=u_jumbled(i);
    end
      Recon_var(index_m)=x_m;
      Recon_var(index_est)=x_e;
   if use_gc
      Recon_var(index_GC_out)=nan;
   else
      Recon_var(setdiff(index_est,index_observable))=nan;
   end
   index_unobservable=sort([index_est(index_unobs);index_unobservable]);
  end
end 

