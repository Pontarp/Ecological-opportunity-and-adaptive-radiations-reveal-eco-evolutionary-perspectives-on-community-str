function [fy] = ode_sys( t , y , p )

%Competitors, GLV
% Set up the matrices of the GLV (see pp 348 in Case)

% k=[p.r_A; p.mu_A];
k=p.g_d_vec;
A=p.A;
% disp_mat=p.disp_mat; 

N=y;

fN=k.*N+(A*N).*N;

fy=fN;
