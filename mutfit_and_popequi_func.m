function [t,y,A,mut_fit, K_A, alpha, a] = mutfit_and_popequi_func( parameters, mutant_trait, mut_type_flag,simsteps, fitcompflag)

%This function computes the equilibrium population sizes by solving the ODE
%system defined by parameters in main and piped into this function

%Extract the parameters that were piped into this script
N0=parameters.N0; 
P0=parameters.P0;
V=parameters.V; %Initial conditions
Z=parameters.Z; %Initial conditions
m=parameters.m; %No consumers
n=parameters.n; %No predators
r_A=parameters.r_A; %Consumer growth
mu_A=parameters.mu_A; %predator death
K0=parameters.K0; %Max K
U=parameters.U; %Resource peaks
sigma_K=parameters.sigma_K; %Resource width
sigma_a=parameters.sigma_a; %Width consumption kernel
bmax=parameters.bmax; %Max attack rate
sigma_b=parameters.sigma_b; %niche width of predator
cP=parameters.cP; %conversion rate
t_end=parameters.t_end; %end time 
abun_vec=parameters.abun_vec; %Abundances in vector for the ODE
growth_death_vec=parameters.growth_death_vec; %Growth and death rates
intrin_growth_comp=parameters.intrin_growth_comp; 
intrin_death_pred=parameters.intrin_death_pred;
mut_N=parameters.mut_N;
mut_P=parameters.mut_P;
dt=parameters.dt;

%%
if fitcompflag==1  %If the function should be used to compute fitness for a given mutant

    %Include mutant
    if mut_type_flag==1 %if the mutant is a consumer
        N0prim=1; 
        N0=[N0prim; N0]; V=[mutant_trait V]; m=length(V); r_A=ones(m,1)*intrin_growth_comp;
        
    elseif mut_type_flag==2 %If the mutant is a predator
            P0prim=1;
            P0=[P0prim; P0]; Z=[mutant_trait Z]; n=length(Z); mu_A=ones(n,1)*intrin_death_pred;
    end
    
    abun_vec=[];
    growth_death_vec=[];
    mut_vec=[];
    trait_vec=[];
    
    abun_vec=[abun_vec; [N0; P0]]; %set the abundance vector
    growth_death_vec=[growth_death_vec; [r_A; mu_A]]; %set up the growth and death rates accordingly
    mut_vec=[mut_vec; ones(m,1)*mut_N; ones(n,1)*mut_P]; %set up mutation vector
    trait_vec=[trait_vec; [V'; Z']];   
    
    
end 
%% Compute and set up interaction and dispersal matrises 

%Compute K for the species 
K_A = K0*exp(-(V-U).^2/2/sigma_K^2);

%Set up the consumer interaction matrix
alpha=[];
for j=1:m
    alpha(j,:)=exp(-(V-V(j)).^2/2/sigma_a^2);  
end

%Set up the matrix for the predator attack rate
a=[]; 
for i=1:m
    for j=1:n
        a(i,j)=-bmax*exp(-(V(i)-Z(j)).^2/2/sigma_b^2); 
    end
end

%Conversion from prey to pred, one per habitat but here we use same
c=zeros(n,n);
c(1:n+1:n*n) = cP; %Convert the diagonal to what ever conversion rate

%% Set up the matrices of the GLV (see pp 348 in Case)
%See also Case pp 387 for an example of the full model with dispersal
%Note that each of the matrices will be devided in subblocks here as we
%have two habitats. 

%Set up the r-alpha-K habitat A subblock for the A matrix
for j=1:m
    A_r_alpha_K(j,:)=-(r_A(j)/K_A(j)).*alpha(j,:);    
end


%Put the whole block matrix together and set up the growth and death rate
%vector
A=[];

%Put the two subblocks together in one block which will be the upper left
%part of the A matrix
r_alpha_K=A_r_alpha_K;


%Set up the ka block in the A matric
a_tmp=a.*-1; %change signe on a
a_tmp=a_tmp'; %use the transpose of a

c=[c]; %put togheter the conversion rate matrices
ka=c*a_tmp; 

%Set up the complete A matrix
A=blkdiag(A,[r_alpha_K a; ka zeros(n,n)]);
    

    %Set parameters and community matrix in a struct that will be piped into
    %the ODE solver
    
    p.g_d_vec=growth_death_vec; %Growth and death rates
    p.A=A; %Pipe the A matrix into the ODE solver

    y0 = abun_vec;
    ode_opts = odeset( 'nonnegative' , 1:length(y0) ); %options for the ode

     if fitcompflag == 0; %Use the ode e solver for the equilibrium comuptation
             %Call the ode45 or ode15s solver and input ode_sys function which contains the diff.
            %equations
            [t,y] = ode15s( @(t,y) ode_sys( t , y , p) , [0,t_end] , y0 , ode_opts );
        elseif fitcompflag == 1;
            %Simulate the equilibrium 
            [t,y]= sim_equilibrium(simsteps, y0, p, dt);
            y=y';
     end
    
    if mut_type_flag==1 %if the mutant is a consumer
        mut_fit= (y(2,1)-y(1,1))/dt;
    elseif mut_type_flag==2 %If the mutant is a predator
        mut_fit= (y(2,m+1)-y(1,m+1))/dt;
    end
  






