
%% EE361 HW#4

%% NAME: _Nail Tosun
%% STUDENT NUMBER: 2094563

function []=template_hw4()

%%
% change your .m file name and the function name above, to the following:
% name_surname_ID_hw4

% please add axis names, legends and titles
% in all your plots

% please be aware that there is MATLAB function at the end of the file
% which you will utilize in part (e)
% http://www.mathworks.com/help/matlab/ref/function.html

% examine the whole template before you start

%%
% PARAMETERS
%define the constant parameters
I = 20; % Amps
Nturn = 220; % turns
A = 8e-4; % m^2
lr = 1e-2; % m
m = 0.4; % kg
kspring = 20; % N/m
u0 = 4*pi*1e-7; % H/m

%define displacement vector in meters
%in the interval of [-0.15,0.15]
%with at least 2000 data points

%%
% part a

%use the following sytax in all your derivations and expressions
x = linspace(-0.15,0.15,5000);
y = linspace(-0.15,0.15,5000);
%%
% $z=sqrt(x^2+y^2)$

%%
% RELUCTANCE
%insert your code here (expression)
l=0.16;
R=(l-x)/(u0*A);
%%
% $R=l/(u_0*A)$
%...

%insert your code here (calculation)
R=(l-x)/(u0*A);
%%
% INDUCTANCE
%insert your code here (expression)

%insert your code here (calculation)
L = Nturn*Nturn./(R);

%%
% part b

%%
% STORED MAGNETIC ENERGY
%insert your code here (derivation and expression)

%insert your code here (calculation)
magneticEnergy = 0.5*I*I.*L;
%insert your code here (plot)
figure
plot(x,magneticEnergy); grid on
title('Stored magnetic energy as a function of x')
xlabel('x distance (m)')
ylabel('Stored Magnetic Energy (Joule)')
%your graph will be here


%%
% part c

%%
% part c-i

%%
% ELECTROMAGNETIC FORCE
%insert your code here (derivation and expression)
%% F = dW/dx Wstored = (0.5)*N^2*I^2*u0*A/((l-x) Fmagnetic = (0.5)*N^2*I^2*u0*A/((l-x)^2
%insert your code here (calculation)
forceMagnetic = ((10)*(20)*Nturn*Nturn*u0*A)./((l-x).^2);
plot(x,forceMagnetic); grid on
title('Mechanical spring force as a function of x');
xlabel('x distance (m)');
ylabel('Mechanical Spring Force (N)');
%%
% part c-ii

%%
% MECHANICAL SPRING FORCE
%insert your code here (expression)

%insert your code here (calculation)
forceSpring = -kspring.*x;
%%
% part c-iii

%%
% NET FORCE
%insert your code here (expression)

%insert your code here (calculation)
forceNet = forceSpring + forceMagnetic;
%%
% part c-iv

%insert your code here (plot)
figure; hold on; grid on
a1 = plot(x,forceSpring,'color','r'); M1 = "Spring Force ";
a2 = plot(x,forceMagnetic,'color','b'); M2 = "Electro-Magnetic Force";
a3 = plot(x,forceNet,'color','g'); M3 = "Net Force";
legend([a1,a2,a3], [M1, M2, M3]);
xlabel('x(m)')
ylabel('Force (N)');
title ('All forces in the system as a function of destination')
%your graph will be here

%%
% part c-v

%intervals of displacement at which the net force is in the direction of +x
for a = 1:5000
    if forceNet(a) < 0
        criticalPoint = (a-2500)*0.3/5000;
        display(criticalPoint);
        display(a);
        break
    end
end

%intervals of displacement at which the net force is in the direction of -x
Area = 0;
deltaX = 0.5/5000;
for i = 2500:5000
    Area = Area + 0.5*deltaX*(forceNet(i)); % Small error is okey due to the integration
    if  Area < 0
        interestPoint = i; %point that mass stop
        display(interestPoint);
        stoppingPoint = (interestPoint - 2500)*(0.3)/5000; % Normalizing 
        display(stoppingPoint);
        break
    end
end
%%
% COMMENT

%Insert
%your
%comment
%here
%and
%there
%...
%% Mass will be ossilate at always +x side which [0,0.0623m].
%%
% part d

%%
% ACCELERATION
%insert your code here (derivation and expression)
    %%F = m*a
%insert your code here (calculation)
massAcceleration = forceNet./m;
%insert your code here (plot)
figure
plot(x,massAcceleration);grid on 
xlabel('x(m)');
ylabel('Mass Accelaration (m/s^2)');
title('Mass Accelaration as a function of destination');
%your graph will be here


%%
% part e

%%
% part e-i
%%
% AGAINST TIME

% insert your code here
%define time vector in seconds
%in the interval of [0,5]
%time step should be 1 msec

% HINTS
%implement a subroutine as a MATLAB function
%which calculates the acceleration of the mass at a given position
%inside the "function" which is given at the end of this file

%implement another subroutine here, which calculates the new position
%of the mass using the acceleration, velocity and previous position data
%at that moment

%the function you have defined will be called by your subroutine
%at each time step

%one way to do this is using "for" loops

%displacement, velocity and acceleration should be treated as vectors (array)


% insert your code here (iterative calculation)
velocity(1) = 0;
timeVector(1) = 0;

for j = 2500:3538
    velocity(j-2498) = sqrt(abs((velocity(j-2499)^2)+deltaX*(massAcceleration(j)+massAcceleration(j+1))));
    timeVector(j-2498) =(velocity(j-2498)-velocity(j-2499))/(0.5*(massAcceleration(j)+massAcceleration(j+1))) + timeVector(j-2499);
    A(j-2498)= 0.5*(massAcceleration(j)+massAcceleration(j+1));
end

for j = 3538:-1:2500
    velocity(-j+4579) = -sqrt(abs((velocity(-j+4578)^2)-deltaX*(massAcceleration(j)+massAcceleration(j+1))));
    timeVector(-j+4579) =(velocity(-j+4579)-velocity(-j+4578))/(0.5*(massAcceleration(j)+massAcceleration(j+1))) + timeVector(-j+4578);
    A(-j+4579) = 0.5*(massAcceleration(j)+massAcceleration(j+1));
end   
for j = 2500:3538
    velocity(j-2498+2078) = sqrt(abs((velocity(j-2499+2078)^2)+deltaX*(massAcceleration(j)+massAcceleration(j+1))));
    timeVector(j-2498+2078) =(velocity(j-2498+2078)-velocity(j-2499+2078))/(0.5*(massAcceleration(j)+massAcceleration(j+1))) + timeVector(j-2499+2078);
    A(j-2498+2078) = 0.5*(massAcceleration(j)+massAcceleration(j+1));
end
for j = 3538:-1:2500
    velocity(-j+4579+2078) = -sqrt(abs((velocity(-j+4578+2078)^2)-deltaX*(massAcceleration(j)+massAcceleration(j+1))));
    timeVector(-j+4579+2078) =(velocity(-j+4579)-velocity(-j+4578+2078))/(0.5*(massAcceleration(j)+massAcceleration(j+1))) + timeVector(-j+4578+2078);
    A(-j+4579+2078) = 0.5*(massAcceleration(j)+massAcceleration(j+1));
end   
for j = 2500:3538
    velocity(j-2498+2078*2) = sqrt(abs((velocity(j-2499+2078)^2)+deltaX*(massAcceleration(j)+massAcceleration(j+1))));
    timeVector(j-2498+2078*2) =(velocity(j-2498+2078*2)-velocity(j-2499+2078*2))/(0.5*(massAcceleration(j)+massAcceleration(j+1))) + timeVector(j-2499+2078*2);
    A(j-2498+2078*2) = 0.5*(massAcceleration(j)+massAcceleration(j+1));
end
for j = 3538:-1:2500
    velocity(-j+4579+2078*2) = -sqrt(abs((velocity(-j+4578+2078*2)^2)-deltaX*(massAcceleration(j)+massAcceleration(j+1))));
    timeVector(-j+4579+2078*2) =(velocity(-j+4579+2078*2)-velocity(-j+4578+2078*2))/(0.5*(massAcceleration(j)+massAcceleration(j+1))) + timeVector(-j+4578+2078*2);
    A(-j+4579+2078*2) = 0.5*(massAcceleration(j)+massAcceleration(j+1));
end   

%%
destinationOfMass(1) = 0;
for m = 1:size(timeVector-1)
    destinationOfMass(1+m) = destinationOfMass(m) + 0.50*((velocity(m)+velocity(m+1))*(timeVector(m+1)-timeVector(m)));
end
destinationOfMass(6235) = destinationOfMass(6234); 

%%
% part e-ii

% insert your code here (plot)

 figure;
 subplot(3,1,1);
 plot(timeVector,destinationOfMass);
 grid on;
 set(gca,'FontSize',12);
 title('Destination of the mass over time');
 xlabel('time (secs)');
 ylabel('Destination (m)');
%
 subplot(3,1,2);
 plot(timeVector, velocity)
 grid on;
 set(gca,'FontSize',12);
 title('Velocity of the mass over time');
 xlabel('time (secs)');
 ylabel('Velocity (m/secs)');
%
 subplot(3,1,3);
 plot(timeVector,A)
 grid on;
 set(gca,'FontSize',12);
 title('Acceleration of the mass over time');
 xlabel('time (secs)');
 ylabel('Velocity (m/secs^2)');
% your graph will be here


%%
% part f

%%
% COMMENT

%Insert
%your
%comment
%here
%and
%there
%...


%%
% part g

%%
% COMMENT

%Insert
%your
%comment
%here
%and
%there
%...


%%
% The following is the MATLAB function which calculates
% the acceleration (output) of the mass at a given position (input)
% you may take advantage of your derivations in previous parts

    function [output] = calculate_acceleration(input)
        
        %%
        
        
    end

%% After you finished
% Run the following command from Matlab terminal (command window)
% generate a report of your .m file as pdf and
% ONLY upload the PDF file to ODTUClass.

%publish('name_surname_ID_hw4.m', 'pdf')

end