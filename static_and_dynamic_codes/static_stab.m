
%This is the function file used in the main file named "Call_static_stab.m".
%Developed by: Mojtaba Kheiri and Shahrzad Tabatabaei
%Affiliation: Concordia University, Mechanical, Industrial & Aerospace Engineering, Montreal, Canada
%Last Modified: May 12, 2021
%For inquiries, please contact us at [sheze.tb@gmail.com]



function omega = static_stab ( typ , N , u )
global stvar

%%%%%%%%%%%%%%%%%%%   System parameters   %%%%%%%%%%%%%%%%%%%


betta = stvar.betta ;
gamma = stvar.gamma ;
f = stvar.f ;
khi_e = stvar.khi_e ;
khi_e_bar = stvar.khi_e_bar ;
epsilon_cn = stvar.epsilon_cn ;
epsilon_ct = stvar.epsilon_ct ;
cb = stvar.cb ;
epsilon_c = stvar.epsilon_c ;
xi = stvar.xi ;
h = stvar.h ;
k0 = stvar.k0;
k0_star = stvar.k0_star;
N = stvar.N ;
theta_knot = stvar.theta_knot ;


%%%%%%%%%%%%%%%%%%%   Pre-setting matrices and arrays   %%%%%%%%%%%%%%%%%%%


Mmat=zeros(N);
Cmat=zeros(N);
Kmat=zeros(N);

phi_1=zeros(1,N);
phi_0=zeros(1,N);
phi_1_prime=zeros(1,N);
phi_0_prime=zeros(1,N);

A = zeros ( 2 * N , 2 * N );
B = zeros ( N , N );
C = zeros ( N , N );
D = zeros ( N , N );
kronecker_delta=eye(N);


    % we use free-free eigenfunctions
    
    lam = [0 0 4.73004 7.85320 10.9956 14.1372 17.2788]; % obtained from the book "INTRODUCTION TO STRUCTURAL DYNAMICS AND AEROELASTICITY, SECOND EDITION"  by Hodges', Table 3.2
    for i = 8 : N
        lam ( i ) = ( 2 * i + 1 ) * pi / 2;
    end
    % free-free euler-bernoulli beam eigenfunctions, obtained from PhD dissertation "Dynamics of flexible slender structures in contact with axial flow and in presence of flexible end-restraints" by Mojtaba Kheiri, eq. 4.32
    sigm(1:2) = 0;
    sigm(3:N) = ( cos ( lam(3:N) ) - cosh ( lam(3:N) ) ) ./ ( sin ( lam(3:N) ) -sinh ( lam(3:N) )  );  %Mojtaba's dissertation....
    phi_1(1,1:2)=[1 sqrt(3)];
    phi_1_prime(1,1:2)=[0 2*sqrt(3)];
    phi_0(1,1:2)=[1 -sqrt(3)];
    phi_0_prime(1,1:2)=[0 2*sqrt(3)];
    
    for r = 3:N
        
        phi_1(1,r) = 2*(-1)^(r-3);
        phi_1_prime(1,r) = 2*(-1)^(r-3)*sigm(r)*lam(r);
        phi_0(1,r)=2;
        phi_0_prime(1,r)= -2*lam(r)*sigm(r);
        
    end
    
    
%%%%%%%%%%%%%%%%%%%   calculation of B C D matrices   %%%%%%%%%%%%%%%%%%%
    
    %B and C (for free-free) from the aforementioned dissertation
    B ( 1:2 , 1:2 )= [0 , 2*sqrt(3) ;0 , 0];
    C ( 1:2 , 1:2 )= [0 , 0 ;0 , 0];
    D ( 1:2 , 1:2 )= [0 , 0 ;0 , 0];
    %vertical boxes
    
    B(3:N,1) = 0;
    B(3:N,2) = 0;
    C(3:N,1) = 0;
    C(3:N,2) = 0;
    D(3:N,1) = 0;
    D(3:N,2) = 0;
    
    %horizontal boxes
    
    for r=3:N
        B ( 1 , r) = 2 * ( ( -1 )^( r - 3 ) -1 );
        B ( 2 , r) = (2 * sqrt(3)) * ( 1 + ( -1 )^( r - 3 ) ) ;
        C ( 1 , r) = 2 *( lam ( r )) * ( sigm ( r ) ) * ( 1 + ( -1 )^( r - 3 ) ) ;
        C ( 2 , r) = (2 * sqrt(3)) * ( ( -1 )^( r - 3 ) -1 ) * ( (sigm(r)) * lam( r ) - 2 );
        D ( 1 , r) = 2 * (( -1 )^( r - 3 )) * ( (sigm ( r )) * lam ( r ) - 1) + 2;
        D ( 2 , r) = (2 * sqrt(3))* ( ( ( -1 )^( r - 3 ) ) * ( (sigm ( r )) * lam ( r ) - 3 ) - 1 );
    end
    
    for s = 3 : N
        for r = 3 : N
            if r==s
                B ( s , r ) = 0;
                C ( s , r ) = ( lam ( r ) * sigm (r) )*( 2 - lam ( r ) * sigm (r) );
                D ( s , r ) = lam(r)* sigm (r) - 0.5 * ( (lam ( r ))^2 ) * ( (sigm (r))^2 );
            else
                B ( s , r ) = 4*( ( (( -1 )^(s+r-6)) - 1) / ( 1 - ( lam ( s ) / lam ( r ) )^4) );
                C ( s , r ) = 4 * (lam ( r ) * sigm (r) - lam ( s ) * sigm (s)) *( ( (( -1 )^(s+r-6)) + 1) / ( 1 - ( lam ( s ) / lam ( r ) )^4) );
                D ( s , r )= ( (phi_1(1,s)) * (phi_1_prime(1,r)) - (phi_1(1,s))* (phi_1(1,r)) + (phi_0(1,r)) * (phi_0(1,s)) - (phi_1_prime(1,s)) * (phi_1(1,r)) - 4 * B ( s , r )* ( ( lam ( s ) / lam ( r ) ) ^ 4 )  ) / ( 1- (lam ( s )/lam ( r )) ^ 4);
            end
        end
    end
    
  
  
  
   
     for s = 1 : N
            for r = 1 : N
               
                    Mmat(s,r) = ( 1 + betta * ( xi - 1 ) ) * kronecker_delta(s,r) + ( 1 + betta * ( f * xi - 1 ) ) * khi_e * phi_1(1,r) * phi_1(1,s);
                    Cmat(s,r) = 2 * xi * u * ( cos( theta_knot ) ) * sqrt ( betta ) * B(s,r) ...
                        + 0.5 * ( kronecker_delta(s,r) ) * ( epsilon_cn * u * ( 1 / ( cos ( theta_knot ) ) ) * sqrt ( betta ) ...
                        + epsilon_c * ( 1 / ( ( cos ( theta_knot ) ) ^ 2 ) ) * sqrt ( betta ) ) ...
                        - ( f * xi * u * ( cos ( theta_knot ) ) * sqrt ( betta ) ...
                        - 0.5 * khi_e_bar * epsilon_cn * u * ( 1 / ( cos ( theta_knot ) ) ) * sqrt ( betta )...
                        - 0.5 * epsilon_c * ( 1 / ( ( cos ( theta_knot ) ) ^ 2 ) ) * khi_e_bar * sqrt ( betta ) ) * phi_1(1,r) * phi_1(1,s)...
                        + xi * khi_e * f * u * ( cos ( theta_knot ) ) * sqrt ( betta ) * phi_1_prime(1,r) * phi_1(1,s);
                    
                    
                    
                    Kmat(s,r) = ( kronecker_delta ( s,r ) ) * ( lam ( r ) ) ^ 4 + xi * ( u ^ 2 ) * ( ( cos( theta_knot ) ) ^ 2 ) * C(s,r)...
                        - ( 0.5 * epsilon_ct * ( u ^ 2 ) * ( 1 + h ) + gamma * ( cos ( theta_knot ) ) ) *( C(s,r) - D(s,r) )...
                        - ( 0.5 * cb * ( u ^ 2 ) + 0.5 * epsilon_ct * ( u ^ 2 ) * khi_e_bar + gamma * khi_e * ( cos ( theta_knot ) ) + 0.5 * epsilon_ct * ( u ^ 2 ) * h * khi_e ) * C(s,r)...
                        + ( 0.5 * epsilon_cn * ( u ^ 2 ) + gamma * ( cos ( theta_knot ) ) + 0.5 * epsilon_ct * u ^ 2 * h ) * B (s,r) ...
                        - ( - 0.5 * khi_e_bar * epsilon_cn * ( u ^ 2 ) ...
                        + f * xi * ( u ^ 2 ) * ( cos ( 2 * ( theta_knot ) ) ) ...
                        - gamma * khi_e * ( cos ( theta_knot ) ) ...
                        - 0.5 * epsilon_ct * ( u ^ 2 ) * h * khi_e ) * phi_1_prime(1,r) * phi_1(1,s)...
                        - ( 0.5 * epsilon_ct * ( u ^ 2 ) * ( 1 + khi_e_bar ) + 0.5 * cb * ( u ^ 2 ) ...
                        + gamma * khi_e * ( cos ( theta_knot ) ) + gamma * ( cos ( theta_knot ) ) ...
                        + 0.5 * epsilon_ct * ( u ^ 2 ) * h ...
                        + 0.5 * epsilon_ct * ( u ^ 2 ) * h * khi_e ) * phi_0(1,s) * phi_0_prime(1,r)...
                        + k0 * phi_0(1,r) * phi_0(1,s) + k0_star * phi_0_prime(1,r) * phi_0_prime(1,s);
     
            end
     end
     
A ( 1 : N , 1 : N ) = 0;
A ( 1 : N , N + 1 : 2  * N ) = eye ( N );
A ( N + 1 : 2 * N , 1 : N ) = - inv ( Mmat ) * Kmat;
A ( N + 1 : 2 * N , N + 1 : 2 * N ) = - inv ( Mmat ) * Cmat;

% V is the matrix of eigenvectors, and D is a diagonal matrix whose
% elements are eigenvalues.
[V,Dmat] = eig ( A );

eigenvals = diag ( Dmat );
eigenfreq = - 1i * eigenvals;
% omega = eigenfreq;
omega = esort ( eigenfreq );
