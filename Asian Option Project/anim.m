function anim(r,F,v)
N=length(F(:,1));
step=round(1+N*v/10);
figure
for i=1:step:N
  plot(r,F(i,:));
  %axis([0 1 0 1/2]);
  drawnow;
  pause(0.3);
end