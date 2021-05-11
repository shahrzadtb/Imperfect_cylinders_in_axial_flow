
clc;
clear;
close all;

N_S = 35 ;
L = 0.4 ;
epsilon_cn = 25.3 * 10 * 0.0100 ;
theta_knot = ( 5 * pi ) / 180 ;



syms y(x)

for mn = 1:50
    u (mn) = 0.01 + (mn-1) * 0.05 ;
    a = - 0.5 * epsilon_cn * ( ( u ( mn ) ) ^ 2 ) ;
    b = - 0.5 * epsilon_cn * ( ( u ( mn ) ) ^ 2 ) * ( ( sin(theta_knot) ) / ( cos(theta_knot) ) ) ;
    ode = diff(y,x,4) + a * diff(y,x) == b;
    ySol(x) = dsolve(ode) ;
    
    yy(mn,1:N_S) = ySol(x) ;
    
    for r = 1 : N_S
        
        xx ( r , 1 ) = ( L / ( 2 * N_S ) ) + ( r - 1 ) * ( L / N_S ) ;
        
    end
 
    plot ( yy(mn,1:N_S) , xx , 'DisplayName','Initial static' )
    %xlabel('displacement y') 
    ylabel('x') 
    axis equal
    view([0 -90])
    xlh = xlabel('displacement y');
    xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2) * 1.08);
    
    hold on
    
    ax = gca;
    %ax.XDir = 'reverse';
    ax.XDir = 'reverse';
end