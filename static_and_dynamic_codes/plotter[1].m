
%This is the plotter file used in the main file named "CallLinCylFullEq.m" which plots the Argand diagrams for the case under study.
%Developed by: Mojtaba Kheiri and Shahrzad Tabatabaei
%Affiliation: Concordia University, Mechanical, Industrial & Aerospace Engineering, Montreal, Canada
%Last Modified: May 12, 2021
%For inquires, please contact: [sheze.tb@gmail.com]



clc;
clear;
close all;
hold on;

global stvar

N = stvar.N;
m = stvar.m;
fileformat = stvar.fformat;

if ( fileformat == 1 )%txt format
    Name1 = 'Argand.txt';
    Name2 = 'Damping.txt';
    Name3 = 'Frequency.txt';
else
    Name1 = 'Argand.plt';
    Name2 = 'Damping.plt';
    Name3 = 'Frequency.plt';
end


fprintf ( 'Argand Diagram\n' );
fid1 = fopen ( Name1 , 'wt' );
fprintf ( fid1 , 'X\t Y\n' );
for i = N : -1 : N - m + 1
    fprintf ( fid1 , '%12.16f\t %12.16f\n', [real( omegamat ( i , : ) )' imag( omegamat ( i , : ) )']' );
end
fclose ( fid1 );


fprintf ( 'Damping Diagram\n' );
fid1 = fopen ( Name2 , 'wt' );
fprintf ( fid1 , 'X\t Y\n' );
for i = N : -1 : N - m + 1
    fprintf ( fid1 , '%12.16f\t %12.16f\n', [Pstp' -imag( omegamat ( i , : ) )']' );
end
fclose ( fid1 );


fprintf ( 'Frequency Diagram\n' );
fid1 = fopen ( Name3 , 'wt' );
fprintf ( fid1 , 'X\t Y\n' );
for i = N : -1 : N - m + 1
    fprintf ( fid1 , '%12.16f\t %12.16f\n', [Pstp' real( omegamat ( i , : ) )']' );
end
fclose ( fid1 );

