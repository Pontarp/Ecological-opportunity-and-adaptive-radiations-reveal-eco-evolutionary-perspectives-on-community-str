function [N0,P0,m,n,V,Z,r_A,mu_A, abun_vec, growth_death_vec, mut_vec,trait_vec, prey_sp_id, pred_sp_id] = uppdate_sys_func(abun_vec,m,n,V,Z,intrin_growth_comp,intrin_death_pred,mut_N,mut_P, prey_sp_id, pred_sp_id)

%This function updates the system
%Recreates N0 and P0 from the abun_vec
%Removes extinct species from NO, P0, V and Z
% set new m and n 
% sets new abun_vec, trait_vec, ....and more


%Uppdate abundance and trait vectors
N0=[]; N0=abun_vec(1:m); %N0 from abunvec
P0=[]; P0=abun_vec(m+1:m+n); %P0 from abunvec

tmp1=find(N0<1);      
N0(tmp1,:)=[]; %Remove N0 element for extinct species
V(tmp1)=[]; %Remove trait of extinct species
prey_sp_id(tmp1)=[];

if length(N0)<1
    disp('WARNING all prey extinct')
end

tmp1=find(P0(:,1)<1);      
P0(tmp1,:)=[]; %Remove P0 element for extinct species
Z(tmp1)=[];  %Remove trait of extinct species
pred_sp_id(tmp1)=[];

%Do not let predator vectors go away, if all predators go extinct keep the
%P0 and Z but with no abundance
if length(P0)<1
    P0=0; Z=0; pred_sp_id=0;
end

m=length(V); n=length(Z); %Update species numbers
r_A=ones(m,1)*intrin_growth_comp; mu_A=ones(n,1)*intrin_death_pred; %Update births and death rate vector 
 
%Update more vectors and parameters
abun_vec=[];
growth_death_vec=[];
mut_vec=[];
trait_vec=[];
 
abun_vec=[abun_vec; [N0; P0]]; %set the abundance vector
growth_death_vec=[growth_death_vec; [r_A; mu_A]]; %set up the growth and death rates accordingly
mut_vec=[mut_vec; ones(m,1)*mut_N; ones(n,1)*mut_P]; %set up mutation vector
trait_vec=[trait_vec; [V'; Z']];
 