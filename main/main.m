cd ../matlab
%% loading the reactions and measured data
reaction_file = '../reactions/reactions.txt'; %reaction set in CID notation saved in reactions.txt
dG0r_measured_file = '../data/reaction_dG.txt'; %
dG0f_measured_file = '../data/formation_dG.txt'; %
pka_measured_file='../data/pKa.mat';

% Reading the reaction file to obtain thermodynamic constraints
fid = fopen(reaction_file);
allData = textscan(fid, '%s', 'Delimiter', '\n');
reactionStrings = allData{1,1};
fclose(fid);
[S, rids, cids] = parseKeggModel(reactionStrings); % to obtain stoichiometric matrix S 

%reading the thermodynamic reaction and formation Gibbs free energy data
fid=fopen(dG0r_measured_file);
allData = textscan(fid, '%s %f %f %f %f %f', 'headerlines',1);
rids_m=(allData{1}); %RIDS of the corresponding available measurements
dG0r_prime=allData{2};
pH_r=allData{3};
T_r=allData{4};
I_r=allData{5};
fclose(fid);

%
fid=fopen(dG0f_measured_file);
allData = textscan(fid, '%s %f %f %f %f %f', 'headerlines',1);
cids_m=(allData{1}); %RIDS of the corresponding available measurements
dG0f_prime=allData{2};
pH_f=allData{3};
T_f=allData{4};
I_f=allData{5};
var_f=allData{6};
fclose(fid);

%Loading the relevant pKa data
load(pka_measured_file);
for i=1:length(cids)
    for j=1:length(kegg_pKa)
        if isequal(char(cids(i)),kegg_pKa(j).cid)
            pKa(i)=kegg_pKa(j);
        end
    end
end

% Determining the measured and unmeasured variables
[~,index_rm]=intersect(rids,rids_m);
[~,index_fm]=intersect(cids,cids_m);
index_rm=sort(index_rm); %indices of rids that correspond to all available measurements
index_fm=sort(index_fm); %indices of cids that correspond to all available measurements

%Inverse Legendre transform to convert all measurements into standard form
S_constraints=[]; %S_constraints generates a constraint matrix using S for each measurement of dG0'_r
for i=1:length(rids_m)
    j=find(strcmp(rids_m(i),rids));
    S_constraints(:,i)=S(:,j);
end
[dG0r_standard,reverse_ddG0r] = legendretransformR(S_constraints,pKa,pH_r,T_r,I_r,dG0r_prime);
[dG0f_standard,reverse_ddG0f] = legendretransformF(length(cids_m),pKa,cids_m,pH_f,T_f,I_f,dG0f_prime);

% Calculating mean values and variances in the measured data

%Reaction Gibbs free energies
dG0r=[];
var_r=[];
for i=1:length(rids(index_rm))
    dG0r=[dG0r;mean(dG0r_standard(find(strcmp(rids(index_rm(i)),rids_m))))];
    var_r=[var_r;var(dG0r_standard(find(strcmp(rids(index_rm(i)),rids_m))))];
end

%Formation Gibbs free energies
dG0f=[];
var_f=[];
for i=1:length(cids(index_fm))
    dG0f=[dG0f;mean(dG0f_standard(find(strcmp(cids(index_fm(i)),cids_m))))];
    var_f=[var_f;var(dG0f_standard(find(strcmp(cids(index_fm(i)),cids_m))))];
end


%calculating variance matrix
variance=[var_r;var_f];
variance(find(variance==0))=1;
sigma=inv(diag(variance)); % inverse of the variance matrix

%Performing reconciliation
%[m,n]=size(S);
%index_m=[index_rm;n+index_fm];
y=[dG0r;dG0f];
[m,n]=size(S);
[reconciled_var,index_unobservable]=recon_l(S,y,index_rm,index_fm,sigma,true);
index_m=[index_rm;index_fm+n];
result.dG0r_standard=reconciled_var(1:n)';
result.dG0f_standard=reconciled_var(n+1:end)';
unobservableRIDS=[];
unobservableCIDS=[];
for i=1:length(index_unobservable)
    if index_unobservable(i)<=n
        unobservableRIDS=[unobservableRIDS;rids(index_unobservable(i))];
    end
end
for i=1:length(index_unobservable)
    if index_unobservable(i)>=n
        unobservableCIDS=[unobservableCIDS;cids(index_unobservable(i)-n)];
    end
end
%RIDS_GC_NA and CIDS_GC_NA refer to the unobservable RIDS and CIDS whose estimates could
%not be obtained using group contribution
if size(find(isnan(result.dG0r_standard)))>0
    result.RIDS_GC_NA=rids(find(isnan(result.dG0r_standard)));
end
if size(find(isnan(result.dG0f_standard)))>0
    result.CIDS_GC_NA=cids(find(isnan(result.dG0f_standard)));
end
result.unobservableRIDS=unobservableRIDS;
result.unobervableCIDS=unobservableCIDS;
result.rids=rids;
result.cids=cids;
%Displaying the result
result






