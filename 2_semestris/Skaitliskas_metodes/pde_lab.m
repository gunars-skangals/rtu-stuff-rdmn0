%% Example 1.
% Solution of the heat equation
clear all, clc
m=0; % the value of m is always zero if the solution is
% obtained in Cartesian coordinates
x=linspace(0,1,20); % We defined the mesh on the interval
% (0,1). The mesh consists of 20 points
% uniformly distributed on the interval (0,1)
% numerical solution using ode45
t=linspace(0,2,5); % We defined the mesh on the interval
% (0,2). The mesh consists of 5 points
% uniformly distributed on the interval (0,2)
u=pdepe(m,@pdeq,@pdeic,@pdebc,x,t);
% Solve the heat equation numerically
% using the given initial and boundary conditions 
%% Example 1 (continuation)
surf(x,t,u)
% plot the solution u(x,t) (the surface in the
% three-dimensional space)
title('Numerical solution')
% Create the title for the graph
xlabel('x')
ylabel('y')
% Create the labels for the axes

%% 
function [c,f,s] = pdeq(x,t,u,dudx)
    c=pi^2;
    f=dudx;
    s=0;
end;
function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t)
    pl=ul;
    ql=0;
    pr=pi*exp(-t);
    qr=1;
end;
function u0 = pdeic(x)
    u0=sin(pi*x);
end;