t = linspace(0,pi)';
x = cos(t);
y = sin(t);
color = jet(50);
figure
hold on
for kk = 1:50
  plot((2+kk)*x,(2+kk)*y,'Color',color(kk,:),'LineWidth',20)
end
ylim([0 15])