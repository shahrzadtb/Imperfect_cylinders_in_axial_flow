%This code solves for the static shape of a flexible cylinder in contact with axial flow
%The cylinder may have an imperfect upstream support and initial inclination with respect to the oncoming flow
%The solution is obtained via a central finite difference method (FDM)
%Developed by: Shahrzad Tabatabaei 
%Affiliation: Concordia University, Mechanical, Industrial & Aerospace Engineering, Montreal, Canada
%Last modified: May 12, 2021
%For inquiries, please contact us at [sheze.tb@gmail.com]
%Copyright of Shahrzad Tabatabaei 2021

clc;
clear;
close all;


%%%%%%%%%%%%%%%%%%%   System parameters   %%%%%%%%%%%%%%%%%%%

L = 0.4 ;    %length of the cylinder, m
gamma = 17.6 ; %dimensionless gravitational parameter
f = 0.8 ; %end-piece shape factor 1
cb = 0.1 ; %end-piece form drag coefficient
khi_e = 0.00792 ; %end-piece shape factor 2
khi_e_bar = 0.01056 ; %end-piece shape factor 3
epsil = 25.3; %slenderness ratio
cn = 0.001; %normal frictional drag coefficient
ct = 0.0125; %tangential frictional drag coefficient
epsilon_cn = epsil * cn ; 
epsilon_ct = epsil * ct ;
xi = 1.22 ; %confinment parameter 1
h = 0.455 ; %confinment parameter 2
theta0 = ( 5 * pi ) / 180 ; %initial inclination angle, rad
k0 = 10 ^ 10 ; %dimensionless translational spring constant
k0_star = 10 ^ 10 ; %dimensionless rotational spring constant

%%%%%%%%%%%%%%%%%%%   Varying parameters   %%%%%%%%%%%%%%%%%%%

N = 35 ;     %number of elements used in the FDM solution
u = 0.001 : 0.01 : 2.97 ;    %nondimensional flow velocity

%%%%%%%%%%%%%%%%%%%   Matrices   %%%%%%%%%%%%%%%%%%%

d = zeros ( N , 1 ) ;
AAmatrix = zeros ( N+4, N+4 ) ;
BBmatrix = zeros (N+4,1) ;
x = zeros (N,1);
c1 = zeros( N , 1 ) ;
c2 = zeros( N , 1 ) ;
c3 = zeros( N , 1 ) ;
displacement = zeros ( 10 , N ) ;

%%%%%%%%%%%%%  discretization along the length  %%%%%%%%%%%%%  

for r = 1 : N
        
        x ( r , 1 ) = ( L / ( 2 * N ) ) + ( r - 1 ) * ( L / N ) ;
        
end
 
zeta = x / L ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

average_theta_s = zeros ( 1 , length ( u ) ) ;
Slope_of_the_line = zeros ( 1 , length ( u ) ) ;
Slope_of_the_line_in_degrees = zeros ( 1 , length ( u ) ) ;
Slope = zeros ( length ( u ) , N ) ;
theta_s_of_tip = zeros ( length ( u ) ) ;

for mn = 1 : length ( u )
  
          a = ( gamma * ( cos(theta0) ) ) + ( 0.5 * epsilon_ct * ( ( u ( mn ) ) ^ 2 ) * h ) + ( 0.5 * epsilon_cn * ( ( u ( mn ) ) ^ 2 ) ) ;
    
    
    for ii = 1 : N 
        
        d ( ii , 1 ) = - ( gamma * ( cos(theta0) ) * ( 1 - ( 1 / ( 2 * N ) ) - ( ( 1 * ( ii - 1 ) ) / N ) ) ) - ( 0.5 * epsilon_ct * ( ( u ( mn ) ) ^ 2 ) * ( 1 - ( 1 / ( 2 * N ) ) - ( ( 1 * ( ii - 1 ) ) / N ) ) ) ...
            + ( xi * ( ( u ( mn ) ) ^ 2 ) * ( ( cos(theta0) ) ^ 2 ) ) -( gamma * ( cos(theta0) ) * khi_e ) - ( 0.5 * cb * ( ( u ( mn ) ) ^ 2 ) )...
            - ( 0.5 * epsilon_ct * ( ( u ( mn ) ) ^ 2 ) * khi_e * h ) - ( 0.5 * epsilon_ct * ( ( u ( mn ) ) ^ 2 ) * h * ( 1 - ( 1 / ( 2 * N ) ) - ( ( 1 * ( ii - 1 ) ) / N ) ) ) - ( 0.5 * epsilon_ct * ( ( u ( mn ) ) ^ 2 ) * khi_e_bar ) ;
       
    end
    
   e = - ( gamma * ( sin(theta0) ) )  - ( 0.5 * epsilon_cn * ( ( u ( mn ) ) ^ 2 ) * ( ( sin(theta0) ) / ( cos(theta0) ) ) )  ;
       
    for j = 1 : N
    
       c1 (j,1) =  ( a / ( 2 * ( N  ^ 3 ) ) ) - 4 + ( d(j,1) / ( N ^ 2 ) ) ;
       c2 (j,1) =  6 - ( ( 2 * d(j,1) ) / ( N ^ 2 ) ) ;
       c3 (j,1) =  - ( ( a / ( 2 * ( N ^ 3 ) ) ) + 4 - ( ( d(j,1) ) / ( N ^ 2 ) ) ) ;
    
    end
    
    c4 =  - ( e / ( N ^ 4 ) ) ;
   
      
    dprime = 3 - ( k0 / ( 2 * ( N ^ 3 ) ) ) + ( 1 / ( N ^ 2 ) ) * ( gamma * ( cos(theta0) ) + 0.5 * epsilon_ct * ( ( u ( mn ) ) ^ 2 ) * ( 1 + khi_e_bar ) ...
             + gamma * ( cos(theta0) ) * khi_e + 0.5 * epsilon_ct * ( ( u ( mn ) ) ^ 2 ) * khi_e * h + 0.5 * epsilon_ct * ( ( u ( mn ) ) ^ 2 ) * h ...
             + 0.5 * cb * ( ( u ( mn ) ) ^ 2 ) ) ;
       
          
    eprime = 3 + ( k0 / ( 2 * ( N ^ 3 ) ) ) + ( 1 / ( N ^ 2 ) ) * ( gamma * ( cos(theta0) ) + 0.5 * epsilon_ct * ( ( u ( mn ) ) ^ 2 ) * ( 1 + khi_e_bar ) ...
             + gamma * ( cos(theta0) ) * khi_e + 0.5 * epsilon_ct * ( ( u ( mn ) ) ^ 2 ) * khi_e * h + 0.5 * epsilon_ct * ( ( u ( mn ) ) ^ 2 ) * h ...
             + 0.5 * cb * ( ( u ( mn ) ) ^ 2 ) ) ;
       
       
    bprime = ( ( f * xi * ( ( u ( mn ) ) ^ 2 ) * ( cos( 2 * theta0) ) ) / ( N ^ 2 ) ) ...
           - ( ( epsilon_cn * ( ( u ( mn ) ) ^ 2 ) * khi_e_bar ) / ( 2 * ( N ^ 2 ) ) ) ...
           - ( ( epsilon_ct * ( ( u ( mn ) ) ^ 2 ) * khi_e * h ) / ( 2 * ( N ^ 2 ) ) ) ...
           - ( ( gamma * ( cos(theta0) ) * khi_e ) / ( N ^ 2 ) ) ;
       
       
       cprime = ( ( gamma * ( sin (theta0) ) * khi_e ) / ( N ^ 3 ) ) - ( ( f * xi * ( ( u ( mn ) ) ^ 2 ) * ( cos( theta0) ) * ( sin( theta0) ) ) / ( N ^ 3 ) ) ...
           + ( ( epsilon_cn * ( ( u ( mn ) ) ^ 2 ) * khi_e_bar * ( ( sin(theta0) ) / ( cos(theta0) ) ) ) / ( 2 * ( N ^ 3 ) ) ) ;
  
    
    c1prime = - ( 1 + ( ( 2 * k0_star ) / N ) ) ;
    c2prime = - ( 1 - ( ( 2 * k0_star ) / N ) ) ;
    c3prime = - dprime ;
    c4prime =  eprime ;
    c5prime = - ( 3 - bprime ) ;
    c6prime = 3 - bprime ;
    c7prime = - cprime ; 
  
  
    AAmatrix ( 1:2 , 1:4 ) = [ 1 , c2prime , c1prime , 1 ; -1 , c4prime , c3prime , 1 ] ;
    AAmatrix ( 3 , 1 ) = 1 ;
    AAmatrix ( 3 , 2 ) = c3 (1,1) ;
    AAmatrix ( 4 , 2 ) = 1 ;
    AAmatrix ( N+2 , N+3 ) = c1 (N,1) ;
    AAmatrix ( N+1 , N+3 ) = 1 ;
    AAmatrix ( N+2 , N+4 ) = 1 ;
    
    
    Amatrix = diag( ones(1,N-2) , 2 ) + diag( ones(1,N-2) , -2 );
    
    for p2 = 1 : N
        
        Amatrix ( p2 , p2 ) = c2 ( p2 , 1) ;
        
    end
    
    for pp = 1: N-1
        
         Amatrix ( pp , pp + 1 ) = c1 ( pp , 1) ;
         Amatrix ( pp + 1 , pp ) = c3 ( pp + 1 , 1) ;
         
    end
    

    
    AAmatrix ( 3:N+2 , 3:N+2) = Amatrix ;

    AAmatrix ( N+3:N+4 , N+1:N+4 ) = [ -1 , c6prime , c5prime , 1 ; 1 , -1 , -1 , 1 ] ;
    Bmatrix = c4 * ones (N,1) ;
    BBmatrix ( 3:N+2 , 1 ) = Bmatrix ;
    BBmatrix ( N+3 , 1 ) = c7prime ;
    eta = AAmatrix \ BBmatrix;   % since nondimensional equation is solved, the obtained value is eta
    eta1 = eta ( 3: N + 2 , 1 ) ;  %The first and second rows are eta(-1) and eta(0).
    y = ( eta1 ) * L ;
    
    %%%%%%%%%%%%%%%%%%    Static shape before rotation    %%%%%%%%%%%%%%%%%%
    
    figure(1)
    displacement(mn,1:N) = y ;
    plot ( displacement(mn,1:N) , x, 'color', [0 0 0])
    ylabel('Y') 
    xlabel('X') 
    axis equal
    view([0 -90])
    
    hold on
    
    %%%% Theta_s of tip is calculated using a straight line, curvefitted to
    %%%% the 1/3 of the end portion of the beam %%%%%%%
    
    q = round ( ( 2 / 3 ) * N ) ;
    displacement1 ( mn , 1 : N - q + 1 ) = y ( q : N ) ;   
    x1 = x ( q : N ) ;
    x11 = x1' ;
    p_line = polyfit ( x11 , displacement1 ( mn , 1 : N - q + 1 ) , 1 ) ;
    Slope_of_the_line( 1, mn ) = p_line ( 1 ) ;
    Slope_of_the_line_in_degrees ( 1 , mn ) = ( ( atan ( Slope_of_the_line( 1 , mn ) ) ) * 180 ) / pi ;
    function_line = polyval ( p_line , x11 ) ; 
    axis equal
    view([0 -90])
    hold on
    
    ax = gca;
    ax.XDir = 'reverse';
    Slope ( mn , 1 ) = ( eta1 ( 1 , 1 ) ) ./ ( zeta ( 1 , 1 ) ) ;
    
    for w = 1 : N - 1
        Slope ( mn , w + 1 ) = ( eta1 ( w + 1 , 1 ) - eta1 ( w , 1 ) ) / ( zeta ( w+1 , 1 ) - zeta ( w , 1 ) ) ;
    end
    
    %Slope = eta1 ./ zeta ;
    theta_s_of_tip ( mn ) = atan ( Slope ( mn , N ) ) ;
    theta_s_of_tip_in_degrees = theta_s_of_tip * (180/pi) ;
    
    figure (2)
    plot ( u(1,mn), displacement(mn,N)/L,'k.','MarkerSize',10 )
    xlabel('u') 
    ylabel('displacement') 
    set(gca,'xtick', (0:0.25:u(mn)))
    set(gca,'XMinorTick','on')
    hold on
    average_theta_s ( 1 , mn ) = mean ( atan ( Slope( mn , 1 : N ) ) ) ;
    average_theta_s_in_degrees = average_theta_s * (180/pi) ;
end

