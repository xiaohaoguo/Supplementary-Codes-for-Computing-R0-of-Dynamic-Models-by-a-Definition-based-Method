close all; clear; clc;
format long;

%% Define all symbolic variables and parameters involved in (SEIAR)_i^n
n = 4;  % number of groups
syms S1 S2 E1 E2 I1 I2 A1 A2 R1 R2;
syms S3 E3 I3 A3 R3;
syms S4 E4 I4 A4 R4;
syms beta11 beta21 beta12 beta22;
syms beta13 beta23 beta33 beta31 beta32;
syms beta14 beta24 beta34 beta44 beta41 beta42 beta43;
syms omega1 omega2 gamma1 gamma2 p kappa f1 f2 f3 f4;


%% Define infected variables
variables = [E1,A1,I1,E2,A2,I2,E3,A3,I3,E4,A4,I4]; % all infected cells
dim = numel(variables);   % dimension, number of infected variables
F = sym(zeros(dim,1));  % initialize the vector of newly infections
V = sym(zeros(dim,1));  % initialize the vector of transitions (exluding newly infections)

% Define vector F and V
% Note that the original ODEs can be written as  dx/dt = F - V
F(1) = beta11*S1*(I1+kappa*A1) + beta21*S1*(I2+kappa*A2) + beta31*S1*(I3+kappa*A3) + beta41*S1*(I4+kappa*A4);
V(1) = p*omega2*E1 + (1-p)*omega1*E1;
V(2) = -p*omega2*E1 + gamma2*A1;
V(3) = -(1-p)*omega1*E1 + (f1+gamma1)*I1;

F(4) = beta22*S2*(I2+kappa*A2) + beta12*S2*(I1+kappa*A1) + beta32*S2*(I3+kappa*A3) + beta42*S2*(I4+kappa*A4);
V(4) = p*omega2*E2 + (1-p)*omega1*E2;
V(5) = -p*omega2*E2 + gamma2*A2;
V(6) = -(1-p)*omega1*E2 + (f2+gamma1)*I2;

F(7) = beta33*S3*(I3+kappa*A3) + beta13*S3*(I1+kappa*A1) + beta23*S3*(I2+kappa*A2) + beta43*S3*(I4+kappa*A4);
V(7) = p*omega2*E3 + (1-p)*omega1*E3;
V(8) = -p*omega2*E3 + gamma2*A3;
V(9) = -(1-p)*omega1*E3 + (f3+gamma1)*I3;

F(10) = beta44*S4*(I4+kappa*A4) + beta14*S4*(I1+kappa*A1) + beta24*S4*(I2+kappa*A2) + beta34*S4*(I3+kappa*A3);
V(10) = p*omega2*E4 + (1-p)*omega1*E4;
V(11) = -p*omega2*E4 + gamma2*A4;
V(12) = -(1-p)*omega1*E4 + (f4+gamma1)*I4;

%% construct the next generation matrix
% compute the jacobian matrices JF, JV via differentiating F, V on variables
JF = sym(zeros(dim)); 
JV = sym(zeros(dim)); 
for i = 1:dim
    for j = 1:dim
        JF(i,j) = diff(F(i),variables(j));
        JV(i,j) = diff(V(i),variables(j));
    end
end

% construct the next-generation matrix Mat = F*V^(-1)
invJV = inv(JV)
Mat = JF * invJV


%% Compute the leading eigenvalues of sub-matrices in the next-generation matrix
% assume that {omega1 omega2 gamma1 gamma2 p kappa} are group irrelevent, and f1 = f2 = f3 = f4 
C = sym(zeros(3));
C(1,1) = (kappa*omega2*p)/(gamma2*(omega1 - omega1*p + omega2*p)) - (omega1*(p - 1))/((f1 + gamma1)*(omega1 - omega1*p + omega2*p));
C(1,2) = kappa/gamma2;
C(1,3) = 1/(gamma1+f1);

% by defining the following symbolic Beta matrix with B(i,j) = betaji * Si, 
% and denote B(i,j) as betaij)
for i = 1:n     
    for j = 1:n
        B(i,j) = eval(['sym(beta',num2str(j),num2str(i),')']);
    end
end

% then the next-generation matrix can be formularized as: Mat = kron(B, C)
% thus the interactive Rij given by eigenvalues of sub-matrices is
% B(i,j)*max(eig(C)), which leads to formula (6) in paper.
%

% since matrix C has only one non-zeros row, i.e. the first row, then C has
% only one non-zero eigenvalue: C(1,1)
% 

% once values of all parameters are obtained, one can use formula (8) in paper:
% R0 = max(eig(eval(B))) * max(eig(eval(C)))
% (note that for R0, we shall set S1 = N1,...,Sn = Nn in substitution)
%


%% Numerical Experiment proves the derivation is correct
% The derivation may be complicate to check from one cell to another,
% therefore, a numerical experiment is designed based on random-generated
% parameters. If the derivation is correct, then the result from our derivation 
% should identical to the result of implicit form.
%

% random generated Beta matrix and parameters
for i = 1:n
    for j = 1:n
        eval(['beta',num2str(j),num2str(i),'= rand(1);']);
    end
end

omega1 = rand(1);
omega2 = rand(1);
gamma1 = rand(1);
gamma2 = rand(1);
p = rand(1);
kappa = rand(1);
f1 = rand(1);
f2 = f1;
f3 = f1;
f4 = f1;
S1 = 1; 
S2 = 1;
S3 = 1;
S4 = 1;



% substitute parameteres (Note that the population size is included in
% betaij in Line 66-67. ) 
eMat = eval(Mat);   
eB = eval(B);   
eC = eval(C);   
eJF = eval(JF);
eJV = eval(JV);
eJV1 = eJV(1:3,1:3);


% differences between the implicit R0 via NGM and the derived explicit R0 via NGM
err = eMat - kron(eB,eC);    
normErr = norm(err)



% all blocks of Next-Generation matrix are identtial
for i = 1:n
    for j = 1:n
        eval(['err_eMat',num2str(i),num2str(j),' = norm(eJF( 3*i-2:3*i , 3*j-2:3*j )*inv(eJV1) - eMat( 3*i-2:3*i , 3*j-2:3*j ))']);
    end
end


% R0 are identical
R1 = max(eig(eMat))
R2 = max(eval(eig(B)))*max(eval(eig(C)))


% Notes: the Matlab function latex() is a useful tool that transform the
% symbolic expression into latex codes, and combine with Mathpix, it will
% be convienient to paste the complecated formula into a Microsoft Word
% document. 
