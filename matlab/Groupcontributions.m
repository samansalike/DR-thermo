function  [G_grp,I,index_out_range]=Groupcontributions(A,G_r,G_f,sigma,G,index_fest) 
% here A is the constraint matrix such that A.x=0 where x is the vector [G_r G_f]' of 
% measured reaction and formation Gibbs energy values.
% index_fest are the indices of the unmeasured formation energies present
% index_out_range are the indices of the unmeasured formation energies that
% cannot be estimated using group contributions due to lack of coverage

if nargin<6
    index_fest=[];
end

%group_cache_file = '../cache/G.mat';
%if ~exist(group_cache_file)
%    fprintf('Group incidence matrix not found in cache \n')
%end
%load(group_cache_file);
G_temp=G;
G_temp(index_fest,:)=[]; %Removing the formation energies that are unmeasured
I=[];
for i=1:size(G_temp,2)
       if ~isempty(find(G_temp(:,i))) %The groups that can be estimated after removing the unmeasured formation energies
       I=[I;i]; %these groups can be estimated from available measurements as they are not out of range
       end
end
G_temp=G_temp(:,I); %removing the columns that are zero vectors
index_out_range=[];
for i=1:length(index_fest)
    grp_index_fest=find(G(index_fest(i),:));
    if isempty(intersect(grp_index_fest,I))
        index_out_range=[index_out_range;index_fest(i)]; 
        %these are the unmeasured formation energies that are not covered 
        %by estimable group contributions
    end
end

n=size(A,1);
[~,gn]=rref(G_temp);
tildeG=G_temp(:,gn); 
Pg=pinv(tildeG)*G_temp;
g_full=size(G_temp,2); %number of groups that can be independently measured
inds=1:g_full;
inds(gn)=[];
g=size(tildeG,2);
x0=zeros((length(G_r)+length(G_f)+g),1);
Z1=zeros(n,g);
Aeq=[A Z1];
Beq=zeros(size(Aeq,1),1);
%%
%performing fmincon
options=optimoptions('fmincon','MaxFunctionEvaluations',5000000,'MaxIterations',1000,'OptimalityTolerance',1e-10);
[x,~]=fmincon(@(x)costFunction(x,G_r,G_f,tildeG,sigma),x0,[],[],Aeq,Beq,[],[],[],options);
G_grp=pinv(Pg)*x(length(G_r)+length(G_f)+1:end);


