

%Plotter

%Developed by: Mojtaba Kheiri

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

