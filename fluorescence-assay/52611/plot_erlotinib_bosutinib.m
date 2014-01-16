% Plot data.
clf;

clear;

x = textread('tmp.txt');
ligand_concentration = 10e-6 * 2.^(-(0:7));

subplot(2,1,1);
plot(x(:,1:5), '.-');
%legend('1 uM', '100 nM', '10 nM', '1 nM', '0 nM');
%legend('1 uM p38', '100 nM p38', '10 nM p38', '1 nM p38', '0 nM p38');
title('erlotinib');
xlabel('ligand 1:2 dilution (from 10 uM)');
ylabel('fluorescence 280-480 nm');
set(gca,'XTick', 1:8);
oldaxis = axis;
axis([0 9 0 oldaxis(4)]);
set(gca,'YTick',[]);

subplot(2,1,2);
plot(x(:,6:10), '.-');
%legend('1 uM', '100 nM', '10 nM', '1 nM', '0 nM');
legend('1 uM p38', '100 nM p38', '10 nM p38', '1 nM p38', '0 nM p38');
title('bosutinib');
xlabel('ligand 1:2 dilution (from 10 uM)');
ylabel('fluorescence 280-480 nm');
set(gca,'XTick', 1:8);
oldaxis = axis;
axis([0 9 0 oldaxis(4)]);
set(gca,'YTick',[]);

scale = 0.6;
exportfig(gcf, 'erlotinib-bosutinib.eps', 'width', scale*3, 'height', scale*6, 'color', 'cmyk');
system('epstopdf erlotinib-bosutinib.eps');
