destinationOfMass(1) = 0;
for m = 1:size(timeVector-1)
    destinationOfMass(1+m) = destinationOfMass(m) + 0.50*((velocity(m)+velocity(m+1))*(timeVector(m+1)-timeVector(m)));
end
destinationOfMass(6235) = destinationOfMass(6234); 

 subplot(3,1,1);
 plot(timeVector,destinationOfMass);
 grid on;
 set(gca,'FontSize',12);
 title('Destination of the mass over time');
 xlabel('time (secs)');
 ylabel('Destination (m)');
 
  figure;
 subplot(3,1,1);
 plot(timeVector,destinationOfMass);
 grid on;
 set(gca,'FontSize',12);
 title('Destination of the mass over time');
 xlabel('time (secs)');
 ylabel('Destination (m)');
%