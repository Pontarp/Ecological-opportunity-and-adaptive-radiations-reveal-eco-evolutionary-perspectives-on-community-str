function main_ecoevo(sigma_a_ind,rep_ind)

%This is the main script for the eco-evolutionary dynamics analysis if a
%diversifying competitive community. 
%The code builds on the Generalized Lotka-Volterra model and predation can
%thus be included. In this analyzes predator abudnance are, however, always
%set to zero leading to a competitive community only.  

rand('seed',fix(sum(1000*clock)));
randn('seed',fix(sum(2000*clock)+1));

name=['mainout_sig_a' num2str(sigma_a_ind) 'rep' num2str(rep_ind) '.mat']; %Set the name for the focal realization and its output

sigma_a_vec=[0.1 0.2 0.3 0.4 0.5 0.6 0.7]; %Range of parameters analyzed
sigma_a=sigma_a_vec(sigma_a_ind);


%% Set the paramters and innitial conditions of the ecological model

%Set habitat variables
%Parameters for the resource distributions in habitats A and B
K0=[10000]; %Max K in habitat 
U=[0]; %Resource peeks in the habitat
sigma_K=[1]; %Width of the resource distribution

%Set consumer variables
N0=[1]; %Row i in matrix denote abundance for species i 
V=[-0]; %Trait, each element, one species
prey_sp_id=[0]; %Set species id as zero for the first population
prey_sp_id_counter=max(prey_sp_id);
m=length(V); %number of competitors
intrin_growth_comp=1; %intrisic growth rate
r_A=ones(m,1)*intrin_growth_comp; %growth rates for each species in one vector 
mut_N=1*10e-3; %mutation rate
sigma_mut_N = 0.02; % standard deviation of consumer mutations

%Set predator variables (Note: Set P0=0 to run prey only adaptive radiation)
P0=[0]; %Row i in matrix denote abundance for species i 
Z=[0]; %Trait, each element, one species

%Note: Lines 44-51 are not relevant for a prey only simulation but included for generality of the code according to the Generalized Lotka-Volterra form 
n=length(Z); %number of predators
intrin_death_pred=-0.2; %intrisic death rate
mu_A=ones(n,1)*intrin_death_pred; %death rates for each species in one vector 
mut_P=1*10e-2; %mutation rate
sigma_mut_P = 0.02; % standard deviation of predator mutations
bmax=0.0001; %Max attack rate (i.e. attack rate when complete trait match between pred and prey)
sigma_b=0.1; %niche width of the predator 
cP=0.3; %Conversion from competitor to predator
pred_sp_id=[0]; %Set species id as zero for the first population
pred_sp_id_counter=max(pred_sp_id);

%Implementation parameters
t_end = 1000000; %end time step in the population equilibrium simulation
dt=0.5; %set the delta time implementation parameter for the population dynamics simulation
epsi=1; %Extinction threshold
t_evo=0; %Evolutonary time, first

%Variables for data output
prey_dist_data{1,1}=1; %time
prey_dist_data{1,2}=[V; N0; prey_sp_id]; %matrix for trait distribution, abundance and species id

prey_phylo_data(1,1)=0; %species id
prey_phylo_data(1,2)=0; %origin
prey_phylo_data(1,3)=0; %time first registered

pred_dist_data{1,1}=1; %time
pred_dist_data{1,2}=[Z; P0; pred_sp_id]; %matrix for trait distribution, abundance and species id

pred_phylo_data(1,1)=0; %species id
pred_phylo_data(1,2)=0; %origin
pred_phylo_data(1,3)=0; %time first registered

com_info_data{1,1}=1; %time
com_info_data{1,2}=0; %K_A
com_info_data{1,3}=0; %alpha
com_info_data{1,4}=0; %a

prey_fitland_data{1,1}=1; %time
prey_fitland_data{1,2}=0; %landscape

pred_fitland_data{1,1}=1; %time
pred_fitland_data{1,2}=0; %landscape



%%
%Plot the resource landscape
count=0; 
range=-3:0.01:3;
for i=range
    count=count+1;
    K_A(count) = K0*exp(-(i-U).^2/2/sigma_K^2);
end

figure(1)
subplot(4,1,1)
plot(range,K_A);
title('resource distribution')

%%
%Set up the abundances, births and deaths, and mutation rates
%according to the ODE system, which will be solved/ simulated
%later
abun_vec=[];
growth_death_vec=[];
mut_vec=[];
trait_vec=[];

abun_vec=[abun_vec; [N0; P0]]; %set the abundance vector
growth_death_vec=[growth_death_vec; [r_A; mu_A]]; %set up the growth and death rates accordingly
mut_vec=[mut_vec; ones(m,1)*mut_N; ones(n,1)*mut_P]; %set up mutation vector
trait_vec=[trait_vec; [V'; Z']];

%Save the parameters in a struct that will be piped into functions later
parameters.N0=N0; 
parameters.P0=P0;
parameters.V=V; %Initial conditions
parameters.Z=Z; %Initial conditions
parameters.m=m; %No consumers
parameters.n=n; %No predators
parameters.r_A=r_A; %Consumer growth
parameters.mu_A=mu_A; %predator death
parameters.K0=K0; %Max K
parameters.U=U; %Resource peak
parameters.sigma_K=sigma_K; %Resource width
parameters.sigma_a=sigma_a; %Width consumption kernel
parameters.bmax=bmax; %Max attack rate
parameters.sigma_b=sigma_b; %niche width of predator
parameters.cP=cP; %conversion rate
parameters.t_end=t_end; %end time 
parameters.abun_vec=abun_vec; %Abundances collected in one vector for the ODE
parameters.growth_death_vec=growth_death_vec; %growth and death rates for the ODE
parameters.intrin_growth_comp=intrin_growth_comp; 
parameters.intrin_death_pred=intrin_death_pred;
parameters.mut_N=mut_N;
parameters.mut_P=mut_P;
parameters.dt=dt; 
%%

mut_range=-3:0.01:3; %a range of trait values that will be used to compute fitness lanscape
 

%% Pipe the parameters into the function that computes equilibrium community
mutant_trait=1; 
mut_type_flag=1; 


dummy=0;
 simsteps=t_end/dt;
[t,y,A,mut_fit, K_A, alpha, a] = mutfit_and_popequi_func(parameters, mutant_trait, mut_type_flag, simsteps, 0);

%Update the abundances according to the computed equilibrium
abun_vec=y(end,:)'; %mean(y(end-100:end,:))'; %take the mean of the last few points in the time series

%Remove extinct species
[mm nn]=size(y); 
for i=1:nn
    if length(find(y(:,i)<=0))>1 %if any of the species were below 1 for the last few time steps
        abun_vec(i)=0;
        disp('Extinction already before evolution starts')
    end
end


%Update the system according to the latest simulated equilibrium
[N0,P0,m,n,V,Z,r_A,mu_A, abun_vec, growth_death_vec, mut_vec,trait_vec, prey_sp_id, pred_sp_id] = uppdate_sys_func(abun_vec,m,n,V,Z,intrin_growth_comp,intrin_death_pred,mut_N,mut_P, prey_sp_id, pred_sp_id);

%Pipe the new variables into the parameters struct for subsequent fitness
%and equilibrium calculations
parameters.N0=N0; 
parameters.P0=P0; 
parameters.V=V;
parameters.Z=Z;
parameters.m=m; %No consumers
parameters.n=n; %No predators
parameters.r_A=r_A; %Consumer growth
parameters.mu_A=mu_A; %Predator death
parameters.abun_vec=abun_vec; %Abundances collected in one vector for the ODE
parameters.growth_death_vec=growth_death_vec; %growth and death rates for the ODE


%%HERE WE START THE EVOLUTOINARY TIMELINE, LOOPING OVER EVOLUTIONARY TIME
%% Loop over evolutionary time
while t_evo < 5000

 %Compute fitness lanscapes
 simsteps=3; 
     count=1;
     for j=mut_range
         mutant_trait=j; 
         
         %Consumer landscape
         mut_type_flag=1;
         
         [t,y,A,mut_fit, K_A, alpha, a] = mutfit_and_popequi_func( parameters, mutant_trait, mut_type_flag, simsteps, 1);
         N_fitland(count)=max(mut_fit,0);
         
         %Predator lanscape
         mut_type_flag=2;
         [t,y,A,mut_fit, K_A, alpha, a] = mutfit_and_popequi_func( parameters, mutant_trait, mut_type_flag, simsteps, 1);
         P_fitland(count)=max(mut_fit,0);
         
         count=count+1; 
     end
  
 if t_evo==0; 
     %Save community info 
    com_info_data{1,1}=t_evo; %time
    com_info_data{1,2}=K_A;
    com_info_data{1,3}=alpha;
    com_info_data{1,4}=a;

    %Save predator and prey fitness landscape
    prey_fitland_data{1,1}=t_evo; %time
    prey_fitland_data{1,2}=N_fitland; %landscape

    pred_fitland_data{1,1}=t_evo; %time
    pred_fitland_data{1,2}=P_fitland; %landscape
 end
    

w=abun_vec.*mut_vec;

w_tot=sum(w);

%Pick the species that should mutate according to w/w_tot
tmp=1:length(w);
mut_sp=tmp(find(rand<cumsum(w/w_tot),1,'first'));
mut_fit=[];

%Find out if the mutating population is a consumer or a predator
%This can be done by looking at the groth_death_vector

if growth_death_vec(mut_sp)>0 %if the mutating species is a competitior
    
    %Get the mutant trait
    mutant_trait=trait_vec(mut_sp) + sigma_mut_N*randn;
    
    %Set the mut_type_flag to 1
    mut_type_flag=1;
    
    %Update evolutionary time
    t_evo=t_evo+1; %OBS IN ITO & DIECKMANN THEY DO (-(1/w_tot)*log(rand));
    
    %Compute invasion fitness
    simsteps=3; 
    [t,y,A,mut_fit, K_A, alpha, a] = mutfit_and_popequi_func( parameters, mutant_trait, mut_type_flag, simsteps, 1);

    %Invade with probability that is based on fitness
    tmp=rand; %draw a random number between 0-1
        
    if mut_fit>0; %tmp < mut_fit/birth_N
        
        disp('consumer mutated with positive fitness')

        V_tmp=V; V_tmp(mut_sp)=mutant_trait; %Switch mutating species trait for mutant trait
        Z_tmp=Z; %create a Z_tmp (not used untill the next update of the system)
                                 
        parameters.V=V_tmp; %pipe the new trait vector into community_equilibrium 
        
         simsteps=t_end/dt;
        [t,y,A,mut_fit, K_A, alpha, a] = mutfit_and_popequi_func( parameters, mutant_trait, mut_type_flag, simsteps, 0);
       
        
        %Save the abundances according to the computed equilibrium
        abun_vec_tmp=y(end,:)'; %mean(y(end-100:end,:))'; %take the mean of the last few points in the time series

        %Remove extinct species
        [mm nn]=size(y); 
        for i=1:nn
            if length(find(y(:,i)<=0))>1 %if any of the species were below 1 for the last few time steps
                abun_vec(i)=0;
            end
        end
        
        %Update the system according to the latest simulated equilibrium
        [N0_tmp,P0_tmp,m_tmp,n_tmp,V_tmp,Z_tmp,r_A_tmp,mu_A_tmp, abun_vec_tmp, growth_death_vec_tmp, mut_vec_tmp,trait_vec_tmp, prey_sp_id_tmp, pred_sp_id_tmp] = uppdate_sys_func(abun_vec_tmp,m,n,V_tmp,Z_tmp,intrin_growth_comp,intrin_death_pred,mut_N,mut_P, prey_sp_id, pred_sp_id);

        parameters.V=V_tmp;
        parameters.Z=Z_tmp;
        parameters.m=m_tmp; %No consumers
        parameters.n=n_tmp; %No predators
        parameters.r_A=r_A_tmp; %Consumer growth
        parameters.mu_A=mu_A_tmp; %Predator death
        parameters.abun_vec=abun_vec_tmp; %Abundances collected in one vector for the ODE
        parameters.growth_death_vec=growth_death_vec_tmp; %growth and death rates for the ODE
        parameters.N0=N0_tmp;
        parameters.P0=P0_tmp;
     
        %Test for mutual invasivability
        %compute mutating species invasion fitness given equilibrium above
        simsteps=3; 
        [t,y,A,mutual_inv, K_A, alpha, a] = mutfit_and_popequi_func( parameters, trait_vec(mut_sp), mut_type_flag, simsteps, 1);
        
        if mutual_inv > 0 %if mutual invasibility exist
            
            disp('consumer mutated with mutual invasibility')
            
            %Include the mutant and keep the mutating species
            N0_tmp=[1; N0]; V_tmp=[mutant_trait V]; m_tmp=length(V_tmp); r_A_tmp=ones(m_tmp,1)*intrin_growth_comp;
            n_tmp=length(Z); prey_sp_id_tmp=[prey_sp_id(mut_sp) prey_sp_id]; pred_sp_id_tmp=pred_sp_id; 
            
            abun_vec_tmp=[];
            growth_death_vec_tmp=[];
            mut_vec_tmp=[];
            trait_vec_tmp=[];
            abun_vec_tmp=[abun_vec_tmp; [N0_tmp; P0]]; %set the abundance vector
            growth_death_vec_tmp=[growth_death_vec_tmp; [r_A_tmp; mu_A]]; %set up the growth and death rates accordingly
            mut_vec_tmp=[mut_vec_tmp; ones(m_tmp,1)*mut_N; ones(n,1)*mut_P]; %set up mutation vector
            trait_vec_tmp=[trait_vec_tmp; [V_tmp'; Z']];
            
            parameters.V=V_tmp;
            parameters.Z=Z;
            parameters.m=m_tmp; %No consumers
            parameters.n=n_tmp; %No predators
            parameters.r_A=r_A_tmp; %Consumer growth
            parameters.mu_A=mu_A; %Predator death
            parameters.abun_vec=abun_vec_tmp; %Abundances collected in one vector for the ODE
            parameters.growth_death_vec=growth_death_vec_tmp; %growth and death rates for the ODE
            parameters.N0=N0_tmp;
            parameters.P0=P0;
            
            %Compute new equilibrium
             simsteps=(t_end/dt)*100; %run simulation 100 times longer as it can take time to find eq
           [t,y,A,mutual_inv, K_A, alpha, a] = mutfit_and_popequi_func( parameters, trait_vec(mut_sp), mut_type_flag, simsteps, 0);
           

           %Save the abundances according to the computed equilibrium
           abun_vec_tmp=y(end,:)'; %mean(y(end-100:end,:))'; %take the mean of the last few points in the time series
           
           %Remove extinct species
            [mm nn]=size(y); 
            for i=1:nn
                if length(find(y(:,i)<=0))>1 %if any of the species were below 1 for the last few time steps
                    abun_vec(i)=0;
                end
            end
            Z_tmp=Z; %reset Z_tmp to Z before updating the system

            %Update the system according to the latest simulated equilibrium            
           [N0_tmp,P0_tmp,m_tmp,n_tmp,V_tmp,Z_tmp,r_A_tmp,mu_A_tmp, abun_vec_tmp, growth_death_vec_tmp, mut_vec_tmp,trait_vec_tmp,prey_sp_id_tmp, pred_sp_id_tmp] = uppdate_sys_func(abun_vec_tmp,m_tmp,n_tmp,V_tmp,Z_tmp,intrin_growth_comp,intrin_death_pred,mut_N,mut_P,prey_sp_id_tmp, pred_sp_id_tmp);
            
        end %end of mutual invasivability check
        
        %Update the system
        N0=N0_tmp; P0=P0_tmp; 
        V=V_tmp; Z=Z_tmp;
        abun_vec=abun_vec_tmp; trait_vec=trait_vec_tmp; 
        m=m_tmp; n=n_tmp; 
        r_A=r_A_tmp; mu_A=mu_A_tmp; 
        growth_death_vec=growth_death_vec_tmp; 
        mut_vec=mut_vec_tmp;
        prey_sp_id=prey_sp_id_tmp; pred_sp_id=pred_sp_id_tmp;

            %Update the species id and phylogenetic information

            V_id_combo=[V; prey_sp_id]; %put trait distribution and species id together
            V_id_combo_sort=sortrows(V_id_combo')'; %Sort the matrix according to the trait distribuiton
            [dummy,p]=sort(V); %p provides the element indeces of V before sorted, will help us turn  the data back to original form

            sp_pres=unique(V_id_combo(2,:)); %Get present species to loop over

            for j=sp_pres %loop over species to check for gaps in them
                sp_j_ind=find(V_id_combo_sort(2,:)==j); %get index of specis j

                V_diff=abs(diff(V_id_combo_sort(1,sp_j_ind))); %compute the distance between each trait value in the j'th species trait distribution
                gap=find(V_diff > sigma_mut_N*3); %find gaps in V that are > x times the size of the mutations

                if length(gap)>0 %if gaps exist in the trait distribution vector

                    tmp1=sp_j_ind(1:gap); %Get elements up to the gap, for the given species j
                    tmp2=sp_j_ind(gap+1:end); %Get elements past the gap, for the given species j

                    V_id_combo_sort(2,tmp1)=prey_sp_id_counter+1; %give the elements up to the gap a new id
                    V_id_combo_sort(2,tmp2)=prey_sp_id_counter+2; %give the elements beyond the gap a new id

                    prey_sp_id_counter=prey_sp_id_counter+2; %update the counter
                    
                    %Fill in phylo data 
                    prey_phylo_data(end+1,1)=prey_sp_id_counter-1; %species id
                    prey_phylo_data(end,2)=j; %origin
                    prey_phylo_data(end,3)=t_evo; %time first registered
                    
                    prey_phylo_data(end+1,1)=prey_sp_id_counter; %species id
                    prey_phylo_data(end,2)=j; %origin
                    prey_phylo_data(end,3)=t_evo; %time first registered
                end
            end

            %Change the data back to original form
            tmp=[p; V_id_combo_sort];
            tmp2=sortrows(tmp',1)';
            prey_sp_id=tmp2(3,:); %prey_sp_id is the only vector that needs updating, N0 and V should be the same

     end %end for check if the mutant should be allowed to invade
               
elseif growth_death_vec(mut_sp)<0 %if the mutating species is a predator and the predator exist
     
    %Get the mutant trait
    mutant_trait=trait_vec(mut_sp) + sigma_mut_P*randn;
    
    %Set the mut_type_flag to 2
    mut_type_flag=2;
    
    %Update evolutionary time
    t_evo=t_evo+1; %OBS IN ITO & DIECKMANN THEY DO (-(1/w_tot)*log(rand));
    
    %Compute mutant fitness
    simtep=3;
    [t,y,A,mut_fit, K_A, alpha, a] = mutfit_and_popequi_func( parameters, mutant_trait, mut_type_flag,simtep, 1);

    %Invade with probability that is based on fitness
    tmp=rand; %draw a random number between 0-1
         
    if mut_fit>0; %tmp < mut_fit/birth_N
        
        disp('predator invaded with positive fitness')
        
        %compute eqilibrium with mutating species switched for the mutant
        mut_sp_tmp=mut_sp-m;
               
        Z_tmp=Z; Z_tmp(mut_sp_tmp)=mutant_trait; %Switch mutating species trait for mutant trait
        V_tmp=V; %create V_tmp (not used untill later)
        
        parameters.Z=Z_tmp; %pipe the new trait vector into community_equilibrium 
         
         simsteps=t_end/dt;
        [t,y,A,mut_fit, K_A, alpha, a] = mutfit_and_popequi_func( parameters, mutant_trait, mut_type_flag, simsteps, 0);
         
        
        %Save the abundances according to the computed equilibrium
        abun_vec_tmp=y(end,:)'; %mean(y(end-100:end,:))'; %take the mean of the last few points in the time series

        %Remove extinct species
        [mm nn]=size(y); 
        for i=1:nn
            if length(find(y(:,i)<=0))>1 %if any of the species were below 1 for the last few time steps
                abun_vec(i)=0;
            end
        end
        
        %Update the system according to the latest simulated equilibrium
        [N0_tmp,P0_tmp,m_tmp,n_tmp,V_tmp,Z_tmp,r_A_tmp,mu_A_tmp, abun_vec_tmp, growth_death_vec_tmp, mut_vec_tmp,trait_vec_tmp, prey_sp_id_tmp, pred_sp_id_tmp] = uppdate_sys_func(abun_vec_tmp,m,n,V_tmp,Z_tmp,intrin_growth_comp,intrin_death_pred,mut_N,mut_P, prey_sp_id, pred_sp_id);


        parameters.V=V_tmp;
        parameters.Z=Z_tmp;
        parameters.m=m_tmp; %No consumers
        parameters.n=n_tmp; %No predators
        parameters.r_A=r_A_tmp; %Consumer growth
        parameters.mu_A=mu_A_tmp; %Predator death        
        parameters.abun_vec = abun_vec_tmp; %Abundances collected in one vector for the ODE
        parameters.growth_death_vec=growth_death_vec_tmp; %growth and death rates for the ODE
        parameters.N0=N0_tmp;
        parameters.P0=P0_tmp;
        
        %Test for mutual invasivability
        %compute mutating species invasion fitness given equilibrium above
        simsteps=3; 
        [t,y,A,mutual_inv, K_A, alpha, a] = mutfit_and_popequi_func( parameters, trait_vec(mut_sp), mut_type_flag, simsteps, 1);
     
               
        if mutual_inv > 0 %if mutual invasibility exist
            
            disp('predator mutated with mutual inv.')
            
            %Include the mutant and keep the mutating species
            P0_tmp=[1; P0]; Z_tmp=[mutant_trait Z]; n_tmp=length(Z_tmp); mu_A_tmp=ones(n_tmp,1)*intrin_death_pred;
            m_tmp=length(V); pred_sp_id_tmp=[pred_sp_id(mut_sp_tmp) pred_sp_id]; prey_sp_id_tmp=prey_sp_id; 
            
            abun_vec_tmp=[];
            growth_death_vec_tmp=[];
            mut_vec_tmp=[];
            trait_vec_tmp=[];
            abun_vec_tmp=[abun_vec_tmp; [N0; P0_tmp]]; %set the abundance vector
            growth_death_vec_tmp=[growth_death_vec_tmp; [r_A; mu_A_tmp]]; %set up the growth and death rates accordingly
            mut_vec_tmp=[mut_vec_tmp; ones(m,1)*mut_N; ones(n_tmp,1)*mut_P]; %set up mutation vector
            trait_vec_tmp=[trait_vec_tmp; [V'; Z_tmp']];
            
            parameters.V=V;
            parameters.Z=Z_tmp;
            parameters.m=m_tmp; %No consumers
            parameters.n=n_tmp; %No predators
            parameters.r_A=r_A; %Consumer growth
            parameters.mu_A=mu_A_tmp; %Predator death
            parameters.abun_vec=abun_vec_tmp; %Abundances collected in one vector for the ODE
            parameters.growth_death_vec=growth_death_vec_tmp; %growth and death rates for the ODE
            parameters.P0=P0_tmp;
            parameters.N0=N0;
            
            %Compute new equilibrium
             simsteps=(t_end/dt)*100; %Run simulation 100 times longer as it can take time to find eq
          [t,y,A,mutual_inv, K_A, alpha, a] = mutfit_and_popequi_func( parameters, trait_vec(mut_sp), mut_type_flag, simsteps, 0);
           
        
          
           %Save the abundances according to the computed equilibrium
           abun_vec_tmp=y(end,:)'; %mean(y(end-100:end,:))'; %take the mean of the last few points in the time series
           
           %Remove extinct species
            [mm nn]=size(y); 
            for i=1:nn
                if length(find(y(:,i)<=0))>1 %if any of the species were below 1 for the last few time steps
                    abun_vec(i)=0;
                end
            end
            V_tmp=V; %Reset V_tmp to V before updating

            %Update the system according to the latest simulated equilibrium            
            [N0_tmp,P0_tmp,m_tmp,n_tmp,V_tmp,Z_tmp,r_A_tmp,mu_A_tmp, abun_vec_tmp, growth_death_vec_tmp, mut_vec_tmp,trait_vec_tmp, prey_sp_id_tmp, pred_sp_id_tmp] = uppdate_sys_func(abun_vec_tmp,m_tmp,n_tmp,V_tmp,Z_tmp,intrin_growth_comp,intrin_death_pred,mut_N,mut_P, prey_sp_id_tmp, pred_sp_id_tmp);
           
        end %end of mutual invasibility check
        
        %Update the system
        N0=N0_tmp; P0=P0_tmp; 
        V=V_tmp; Z=Z_tmp;
        abun_vec=abun_vec_tmp; trait_vec=trait_vec_tmp; 
        m=m_tmp; n=n_tmp; 
        r_A=r_A_tmp; mu_A=mu_A_tmp; 
        growth_death_vec=growth_death_vec_tmp; 
        mut_vec=mut_vec_tmp;
        prey_sp_id=prey_sp_id_tmp; pred_sp_id=pred_sp_id_tmp;
        
        
            %Update the species id and phylogenetic inforamation

            Z_id_combo=[Z; pred_sp_id]; %put trait distribution and species id together
            Z_id_combo_sort=sortrows(Z_id_combo')'; %Sort the matrix according to the trait distribuiton
            [dummy,p]=sort(Z); %p provides the element indeces of V before sorted, will help us turn  the data back to original form

            sp_pres=unique(Z_id_combo(2,:)); %Get present species to loop over

            for j=sp_pres %loop over species to check for gaps in them
                sp_j_ind=find(Z_id_combo_sort(2,:)==j); %get index of specis j

                Z_diff=abs(diff(Z_id_combo_sort(1,sp_j_ind))); %compute the distance between each trait value in the j'th species trait distribution
                gap=find(Z_diff > sigma_mut_P*3); %find gaps in V that are > x times the size of the mutations

                if length(gap)>0 %if gaps exist in the trait distribution vector

                    tmp1=sp_j_ind(1:gap); %Get elements up to the gap, for the given species j
                    tmp2=sp_j_ind(gap+1:end); %Get elements past the gap, for the given species j

                    Z_id_combo_sort(2,tmp1)=pred_sp_id_counter+1; %give the elements up to the gap a new id
                    Z_id_combo_sort(2,tmp2)=pred_sp_id_counter+2; %give the elements beyond the gap a new id

                    pred_sp_id_counter=pred_sp_id_counter+2; %update the counter
                    
                    %Fill in phylo data 
                    pred_phylo_data(end+1,1)=pred_sp_id_counter-1; %species id
                    pred_phylo_data(end,2)=j; %origin
                    pred_phylo_data(end,3)=t_evo; %time first registered
                    
                    pred_phylo_data(end+1,1)=pred_sp_id_counter; %species id
                    pred_phylo_data(end,2)=j; %origin
                    pred_phylo_data(end,3)=t_evo; %time first registered
                end
            end

            %Change the data back to original form
            tmp=[p; Z_id_combo_sort];
            tmp2=sortrows(tmp',1)';
            pred_sp_id=tmp2(3,:); %prey_sp_id is the only vector that needs updating, N0 and V should be the same

    end %end for check if the mutant should be allowed to invade
             
end %end of if/ elseif loop

%Save data 
prey_dist_data{end+1,1}=t_evo; %time
prey_dist_data{end,2}=[V; N0'; prey_sp_id]; %matrix for trait distribution, abundance and species id
% pred_dist_data{end+1,1}=t_evo; %time
% pred_dist_data{end,2}=[Z; P0'; pred_sp_id]; %matrix for trait distribution, abundance and species id

com_info_data{end+1,1}=t_evo; %time
com_info_data{end,2}=K_A;
com_info_data{end,3}=alpha;
com_info_data{end,4}=a;

prey_fitland_data{end+1,1}=t_evo; %time
prey_fitland_data{end,2}=N_fitland; %landscape
% pred_fitland_data{end+1,1}=t_evo; %time
% pred_fitland_data{end,2}=P_fitland; %landscape


%Update parameters
parameters.V=V;
parameters.Z=Z;
parameters.m=m; %No consumers
parameters.n=n; %No predators
parameters.r_A=r_A; %Consumer growth
parameters.mu_A=mu_A; %Predator death
parameters.abun_vec=abun_vec; %Abundances collected in one vector for the ODE
parameters.growth_death_vec=growth_death_vec; %growth and death rates for the ODE
parameters.N0=N0;
parameters.P0=P0;


%Do some plotting every n time steps
if mod(t_evo,50)==1
    figure(1)
    subplot(4,1,2)
    cla
    hold on
    plot(mut_range,N_fitland,'color',[0.4 0.4 0.4])
    %xlim([-1 3])
    ylabel('Cons. Fitness')

%     subplot(4,1,2)
%     hold on
%     plot(mut_range,P_fitland,'r')
%     %xlim([-1 3])
%     %ylim([-0.2 1])
%     ylabel('Pred. Fitness')
%     %axis('square')

    subplot(4,1,3)
    hold on
    plot(parameters.V,t_evo,'.','color',[0.4 0.4 0.4])
    xlim([-3 3])
    ylabel('Evo. time')

%     subplot(4,1,3)
%     plot(parameters.Z,t_evo,'.r')
%     xlim([-3 3])
%     %axis('square')
%     xlabel('Trait value')

    subplot(4,1,4)
    cla
    hold on
    for i=1:length(V)
        bar(V(i),N0(i),'facecolor',[0.4 0.4 0.4],'barwidth',0.1)
    end
    for i=1:length(Z)
        bar(Z(i),P0(i),'facecolor','r','barwidth',0.1)
    end
    xlim([-3 3])
    ylabel('N*')
    xlabel('Trait')
    title('Trait distribution')
end
    
end %end of evolutonary time while loop

%Save output
save(name,'prey_dist_data','prey_phylo_data','pred_dist_data','pred_phylo_data', 'com_info_data', 'prey_fitland_data', 'pred_fitland_data');

%Save figures
for i=[1]
   saveas(i,['mainout_sig_a' num2str(sigma_a_ind) 'rep' num2str(rep_ind) '_fig_' num2str(i)])
end

