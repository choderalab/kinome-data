% Plot data.
clf;

clear;

x = textread('tmp.txt');
ligand_concentration = 10e-6 * 2.^(-(0:7));

subplot(2,2,1);
plot(x(:,1:5), '.-');
legend('1 uM', '100 nM', '10 nM', '1 nM', '0 nM');
title('erlotinib');
xlabel('ligand 1:2 dilution');
ylabel('fluorescence 280/480 nm');
set(gca,'XTick', 1:8);
oldaxis = axis;
axis([0 9 0 oldaxis(4)]);
% Plot linear prediction for each increment.
hold on;
for i = 2:8
  prediction = 2*(x(i,1) - x(i,5)) + x(i-1,5);
  plot([i-1 i], [prediction x(i,1)], 'k-');
end

subplot(2,2,2);
plot(x(:,6:10), '.-');
legend('1 uM', '100 nM', '10 nM', '1 nM', '0 nM');
title('bosutinib');
xlabel('ligand 1:2 dilution');
ylabel('fluorescence 280/480 nm');
set(gca,'XTick', 1:8);
oldaxis = axis;
axis([0 9 0 oldaxis(4)]);
% Plot linear prediction for each increment.
hold on;
for i = 2:8
  prediction = 2*(x(i,6) - x(i,10)) + x(i-1,10);
  plot([i-1 i], [prediction x(i,6)], 'k-');
end

subplot(2,2,3);
plot(x(:,11), '.-');
title('2,4-diaminoquinazoline');
xlabel('1:2 dilution');
ylabel('fluorescence 280/480 nm');
set(gca,'XTick', 1:8);
oldaxis = axis;
axis([0 9 0 oldaxis(4)]);

subplot(2,2,4);
plot(x(:,12), '.');
title('2,4-diaminoquinazoline');
xlabel('protein (1-4) and buffer (5-8) control wells');
ylabel('fluorescence 280/480 nm');
set(gca,'XTick', 1:8);
oldaxis = axis;
axis([0 9 0 oldaxis(4)]);

print -depsc badass-data.eps
system('epstopdf badass-data.eps');
