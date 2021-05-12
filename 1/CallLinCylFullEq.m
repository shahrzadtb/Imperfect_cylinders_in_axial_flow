

%This code solves the dynamic equation of the mentioned case, to obtain the
%divergence speed
%Developed by: Mojtaba Kheiri & later used and modified by Shahrzad Tabatabaei
%Last Modified: May 12, 2021

clc;
clear;
close all;
hold on;

global stvar


u = 0.001 : 0.01 : 20.0; %flow velocity range
stvar.N = 10; %number of eignmodes
stvar.typ = 6 ; %type of the configuration

stvar.fformat = 1; %txt file output in plotter
stvar.m = 5; %number of modes considered in Tecplot output

%systme parameters

stvar.betta = 0.5 ;   %mass ratio
stvar.gamma = 17.6 ;     %gravity effect
stvar.khi_e = 0.00792 ;
stvar.khi_e_bar = 0.01056 ; 
stvar.epsilon_cn = 25.3 * 0.0100 ;
stvar.epsilon_ct = 25.3 * 0.0125 ;
stvar.xi = 1.22 ; % confinement parameter
stvar.alpha_star = 0.0 ; % viscoelastic damping
stvar.mui_star = 0.0 ; % hysteretic damping
stvar.h = 0.455 ; % diameter of cylinder/Dh
stvar.epsilon_c = 0.0 ; %viscous forces effect at zero flow velocity   
stvar.f = 0.8 ; %fairly streamlined
stvar.cb = 0.1 ;
stvar.k0_star = 10 ^ 10 ;
stvar.k0 = 10 ^ 10 ;



% Static equation values %

stvar.N_S = 35 ;  % number of nodes on the beam
stvar.theta_knot = ( 5 * pi ) / 180 ;   % initial inclination angle



figure ( 1 );
omegamat = zeros ( 2 * stvar.N , length ( u ) );


for i = 1 : length ( u )
    omega = LinCylinStaticFullEqu( stvar.typ , stvar.N , u ( i ) );
    omegamat ( : , i ) = omega;
    plot(real(omega(stvar.N/2:3*stvar.N/2+1,:)),imag(omega(stvar.N/2:3*stvar.N/2+1,:)), 'go');
    hold on;
end

grid on
%prepare for tecplot
Pstp = u;
plotter;




