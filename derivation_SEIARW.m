close all; clear all; clc;

%% Define all symbolic variables and parameters involved in SEIARW model
syms beta beta_w p omega1 omega2 mu1 mu2 gamma1 gamma2 epsilon k dr;
syms S E I A R W N;

%% Define infected variables
variables = [E, I, A, W]; % all infected cells
dim = numel(variables);  
F = sym(zeros(dim,1));  % initialize the vector of newly infections
V = sym(zeros(dim,1));  % initialize the vector of transitions (exluding newly infections)

% Define vector F and V
% Note that the original ODEs can be written as  dx/dt = F - V
F(1) = beta*S*(I+k*A) + beta_w*S*W;
V(1) = (1-p)*omega1*E + p*omega2*E + dr*E;
V(2) = -(1-p)*omega1*E + gamma1*I + dr*I;
V(3) = -p*omega2*E + gamma2*A + dr*A;
V(4) = -mu1*I - mu2*A + epsilon*W;



%% construct the next generation matrix and compute its leading eigenvalue.
% compute the jacobian matrix 
JF = sym(zeros(dim)); 
JV = sym(zeros(dim)); 
for i = 1:dim
    for j = 1:dim
        JF(i,j) = diff(F(i),variables(j));
        JV(i,j) = diff(V(i),variables(j));
    end
end
JF
JV

%% the next-generation matrix F*V^(-1)
% the matrix is in fact the reproduction matrix for interactive
% transmissibility
Mat = JF*inv(JV)

% all eigenvalues of the Next-Generation matrix
eigenvalues = eig(Mat)

% maximum real part of eigenvalues 
Reff = eigenvalues(end)

% R0 is obtained at the disease free equalibrium (N, 0, 0, 0, 0, 0, 0)
R0 = subs(Reff,S,N)

% latex code
latex(R0)




