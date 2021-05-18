%This is the main file 
%This code uses an eigenvalue problem solver to find the critical flow velocity for the divergence of a flexible cylinder in contact with axial flow
%The cylinder may have an imperfect upstream support and initial inclination with respect to the oncoming flow
%Developed by: Mojtaba Kheiri and Shahrzad Tabatabaei
%Affiliation: Concordia University, Mechanical, Industrial & Aerospace Engineering, Montreal, Canada
%Last Modified: May 12, 2021
%For inquiries, please contact us at [sheze.tb@gmail.com]


clc;
clear;
close all;
hold on;

global stvar

%%%%%%%%%%%%%%%%%%%   System parameters   %%%%%%%%%%%%%%%%%%%

stvar.betta = 0.5 ; %mass ratio
stvar.gamma = 17.6 ; %dimensionless gravitational parameter
stvar.khi_e = 0.00792 ;%end-piece shape factor 2
stvar.khi_e_bar = 0.01056 ; %end-piece shape factor 3
stvar.typ = 6 ; %type of the configuration, 6 indicates geometrically imperfect with imperfect support
epsil = 25.3; %slenderness ratio
cn = 0.001; %normal frictional drag coefficient
ct = 0.0125; %tangential frictional drag coefficient
stvar.epsilon_cn = epsil * cn ; 
stvar.epsilon_ct = epsil * ct ;
stvar.xi = 1.22 ; %confinment parameter 1
stvar.alpha_star = 0.0 ; %viscoelastic damping
stvar.mui_star = 0.0 ; %hysteretic damping
stvar.h = 0.455 ; %confinment parameter 2
stvar.epsilon_c = 0.0 ; %viscous forces effect at zero flow velocity   
stvar.f = 0.8 ; %fairly streamlined
stvar.cb = 0.1 ; %end-piece form drag coefficient
stvar.k0_star = 10 ^ 10 ; %dimensionless rotational spring constant
stvar.k0 = 10 ^ 10 ; %dimensionless translational spring constant
stvar.theta_knot = ( 5 * pi ) / 180 ; %initial inclination angle

%%%%%%%%%%%%%%%%%%%   varying parameters   %%%%%%%%%%%%%%%%%%%


u = 0.001 : 0.01 : 20.0;    %nondimensional flow velocity
stvar.N = 10;    %number of eignmodes
stvar.N_S = 35 ;     %number of elements used in the FDM solution

stvar.fformat = 1; %txt file output in plotter
stvar.m = 5; %number of modes considered in Tecplot output



figure ( 1 );
omegamat = zeros ( 2 * stvar.N , length ( u ) );


%%%%%%%%%%%%%%%%%%%   calling static_stab   %%%%%%%%%%%%%%%%%%%


for i = 1 : length ( u )
    omega = static_stab ( stvar.typ , stvar.N , u ( i ) );
    omegamat ( : , i ) = omega;
    plot(real(omega(stvar.N/2:3*stvar.N/2+1,:)),imag(omega(stvar.N/2:3*stvar.N/2+1,:)), 'go');
    hold on;
end

grid on


%%%%%%%%%%%%%%%%%%%   prepare for Tecplot   %%%%%%%%%%%%%%%%%%%


Pstp = u;
plotter;




