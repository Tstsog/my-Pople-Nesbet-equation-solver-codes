% This Matlab code computes the ground state energy for nitrogen (N) atom by solving the Pople-Nesbet equation using the
% 6-31G basis set ([3s2p]) (unrestricted Hartree-Fock (uhf) calculation). 
%
% The core Hamiltonian matrix (H_core), overlap matrix (S_ov) and two-electron integrals (tei) (Li_6_31G_tei.txt) are computed 
% by my own developing code. An obtained total energy is compared with
% those of https://cccbdb.nist.gov/. 
% 
%
% Ref: A. Szabo and N. S. Ostlund "Modern Quantum Chemistry" book.  
%
% Written by Tsogbayar Tsednee (PhD)
% Email: tsog215@gmail.com
% March 13, 2025 & University of North Dakota 
%
function [] = N_uhf_6_31G
clear N_uhf_6_31G; clc
%
format short 
%
S_ov = [
  1.         0.22221582 0.         0.         0.         0.17186917  0.         0.         0.        ;
  0.22221582 1.         0.         0.         0.         0.77550678  0.         0.         0.        ;
  0.         0.         1.         0.         0.         0.          0.51929556 0.         0.        ;
  0.         0.         0.         1.         0.         0.          0.         0.51929556 0.        ;
  0.         0.         0.         0.         1.         0.          0.         0.         0.51929556;
  0.17186917 0.77550678 0.         0.         0.         1.          0.         0.         0.        ;
  0.         0.         0.51929556 0.         0.         0.          1.         0.         0.        ;
  0.         0.         0.         0.51929556 0.         0.          0.         1.         0.        ;
  0.         0.         0.         0.         0.51929556 0.          0.         0.         1.        ];
%
H_core = [
  -24.44343651  -5.12149184   0.           0.           0.   -4.05443297   0.           0.           0.        ;
   -5.12149184  -6.61203396   0.           0.           0.   -5.03549926   0.           0.           0.        ;
    0.           0.          -5.60525108   0.           0.    0.          -2.5397417    0.           0.        ;
    0.           0.           0.          -5.60525108   0.    0.           0.          -2.5397417    0.        ;
    0.           0.           0.           0.          -5.60525108         0.           0.           0.          -2.5397417 ;
   -4.05443297  -5.03549926   0.           0.           0.                -4.82556578   0.           0.           0.        ;
    0.           0.          -2.5397417    0.           0.                 0.          -2.89899661   0.           0.        ;
    0.           0.           0.          -2.5397417    0.                 0.           0.          -2.89899661   0.        ;
    0.           0.           0.           0.          -2.5397417          0.           0.           0.          -2.89899661];
%
%
N_el = 7.;
N_a = 5;
N_b = N_el - N_a;
dim = 9;      % size of basis sets
itermax = 300; tol = 1e-12;
%
tei_n = 6561;
%
read_tei_data = fopen('N_6_31G_tei.txt', 'r');
tei_data_n5 = textscan(read_tei_data, '%d %d %d %d %f');
%
p = zeros(tei_n,1); q = zeros(tei_n,1); r = zeros(tei_n,1); s = zeros(tei_n,1); vals = zeros(tei_n,1);
p(1:tei_n) = tei_data_n5{1};
q(1:tei_n) = tei_data_n5{2};
r(1:tei_n) = tei_data_n5{3};
s(1:tei_n) = tei_data_n5{4};
vals(1:tei_n) = tei_data_n5{5};
for i = 1:tei_n
    tei(p(i),q(i),r(i),s(i)) = vals(i);
%    tei(q(i),p(i),r(i),s(i)) = vals(i);    
%    tei(p(i),q(i),s(i),r(i)) = vals(i);    
%    tei(q(i),p(i),s(i),r(i)) = vals(i);   
    %
%    tei(r(i),s(i),p(i),q(i)) = vals(i);    
%    tei(s(i),r(i),p(i),q(i)) = vals(i);        
%    tei(r(i),s(i),q(i),p(i)) = vals(i);        
%    tei(s(i),r(i),q(i),p(i)) = vals(i);            
end
%
Q_tei = tei;
%
P_old_a = 0.5 * ones(dim,dim); % initial charge population for alpha-spin electron
P_old_b = 0.5 * ones(dim,dim); % initial charge population for beta-spin electron
P_T = 0.5 * ones(dim,dim); 
%
for iter = 1:itermax
    iter
    P_a = P_old_a;
    P_b = P_old_b;
    %
    F_a = H_core;
    F_b = H_core;
    for p = 1:dim
        for q = 1:dim
            for r = 1:dim
                for s = 1:dim
                    F_a(p,q) = F_a(p,q) + (P_T(r,s) * Q_tei(p,q,r,s) - P_a(r,s) * Q_tei(p,r,q,s));  % Ref: Eq. (3.348)
                    F_b(p,q) = F_b(p,q) + (P_T(r,s) * Q_tei(p,q,r,s) - P_b(r,s) * Q_tei(p,r,q,s));  % Ref: Eq. (3.349)                    
                end
    
            end
    
        end
    end
    Ham_fock_a = F_a ;     % Fock matrix
    S_mat_fock = S_ov;

    [Vec_a,En_a] = eig(Ham_fock_a,S_mat_fock);                                     % % Ref: Eq. (3.351)
    En_a = diag(En_a);
    [foo, ij] = sort(En_a);
    En_a = En_a(ij);
    [En_a(1), En_a(2)];  % orbital energies for alpha-spin electron
    %
    Vec_a = Vec_a(:,ij);                       % expansion coefficients 
    %
    for i = 1:dim
        norm = 0.;
        for p = 1:dim
            for q = 1:dim
                norm = norm + Vec_a(p,i) * Vec_a(q,i) * S_ov(p,q);
            end
        end
        Vec_a(:,i) = Vec_a(:,i)/sqrt(norm);
    end
    %
    P_new_a = zeros(dim,dim);
    for i = 1:N_a
        for pp = 1:dim
            for qq = 1:dim
                P_new_a(pp,qq) = P_new_a(pp,qq) + Vec_a(pp,i)*Vec_a(qq,i);
            end
        end
    end
    %%% beta spin
    Ham_fock_b = F_b ;     % Fock matrix
    S_mat_fock = S_ov;

    [Vec_b,En_b] = eig(Ham_fock_b,S_mat_fock);                                     % Eigenvalue problem: F*c = En*S*c % Ref: Eq. (3.352)
     En_b = diag(En_b);
    [foo, ij] = sort(En_b);
     En_b = En_b(ij);
    [En_b(1), En_b(2)];  % orbital energies for beta-spin electron
    %
    Vec_b = Vec_b(:,ij);                       % expansion coefficients 
    %
    for i = 1:dim
        norm = 0.;
        for p = 1:dim
            for q = 1:dim
                norm = norm + Vec_b(p,i) * Vec_b(q,i) * S_ov(p,q);
            end
        end
        Vec_b(:,i) = Vec_b(:,i)/sqrt(norm);
    end
    %
    P_new_b = zeros(dim,dim);
    for i = 1:N_b
        for pp = 1:dim
            for qq = 1:dim
                P_new_b(pp,qq) = P_new_b(pp,qq) + Vec_b(pp,i)*Vec_b(qq,i);
            end
        end
    end
    %
    P_T = P_new_a + P_new_b;

    %
     if ( (abs(sum(sum(P_new_a-P_old_a)))) && (abs(sum(sum(P_new_b-P_old_b)))) < tol)
            break 
     end
    %        
    P_old_a = P_new_a;
    P_old_b = P_new_b;

end
%%%

[En_a';
 En_b']
En_0 = sum(0.5*diag( P_T(:,:)*H_core(:,:) + P_a(:,:)*F_a(:,:) + P_b(:,:)*F_b(:,:) ) ) % -54.385007714099110 ground state energy in atomic unit

%%%
%[En_a';
% En_b']
% -15.6766   -1.1587   -0.5674   -0.5674   -0.5674    0.8449    0.8449    0.8449    0.9262
% -15.5858   -0.7257    0.1620    0.1620    0.1620    1.0571    1.0591    1.0591    1.0591

%
%En_0 = -54.385007714099110 vs -54.385008 = https://cccbdb.nist.gov/

%%%
return
end
