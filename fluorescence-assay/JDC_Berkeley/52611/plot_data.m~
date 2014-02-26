% Plot data.
clf;

x = textread('tmp.txt');
ligand_concentration = 10e-6 * 2.^(-(0:7));

baseline = x(:,5);
signal = x(:,1) - x(:,5);

% Fit data.
complex_concentration = @(Pt,Lt,Kd) ((Pt + Lt + Kd) - sqrt((Pt + Lt + Kd).^2 - 4*Pt*Lt)) / 2;
model = @(Pt,Lt,a,b,Kd) a*complex_concentration(Pt,Lt,Kd) + b;
errfun = @(X) sum((model(X(1),ligand_concentration',X(2),X(3),X(4)) - signal).^2);
X0 = [1e-6 1e10 3000 250e-9];
X = fminsearch(errfun, X0);
Pt = X(1)
a = X(2)
b = X(3)
Kd = X(4)

Pt = 1.0e-6;
errfun = @(X) sum((model(Pt,ligand_concentration',X(1),X(2),X(3)) - signal).^2);
X0 = [1e10 3000 250e-9];
X = fminsearch(errfun, X0);
a = X(1)
b = X(2)
Kd = X(3)

Pt = 1.0e-6;
Kd = 250e-9;
errfun = @(X) sum((model(Pt,ligand_concentration',X(1),X(2),Kd) - signal).^2);
X0 = [(x(1,1)-x(8,1))/Pt x(8,1)];
X = fminsearch(errfun, X0);
a = X(1)
b = X(2)

disp('-------------------------------');
% a [PL] + b [L] + c [P] + d
Kd = 10000e-9;
Pt = 1.0e-6;
d = mean(x(5:8,12));
c = (x(1,12) - d)/Pt;
b = (x(1,5) - d)/ligand_concentration(1);
a = 1.3*(x(1,1)-b*(ligand_concentration(1)-Pt)-c*Pt-d)/Pt;
options = optimset;
options.TolX = 1e-50;
options.TolFun = 1e-20;
options.Display = 'iter';
X0 = [Kd Pt a b c d];
X = X0;
model = @(Kd,Pt,a,b,c,d) a*complex_concentration(Pt,ligand_concentration',Kd) + b*(ligand_concentration' - complex_concentration(Pt,ligand_concentration',Kd)) + c*(Pt - complex_concentration(Pt,ligand_concentration',Kd)) + d;
errfun = @(X) sum((X(3)*complex_concentration(X(2),ligand_concentration',X(1)) + X(4)*(ligand_concentration' - complex_concentration(X(2),ligand_concentration',X(1))) + X(5)*(X(2) - complex_concentration(X(2),ligand_concentration',X(1))) + X(6) - x(:,1)).^2) + sum((X(4)*ligand_concentration' + X(6) - x(:,5)).^2) + sum((X(5)*X(2) - x(1,12)).^2) + sum((X(6) - x(5:8,12)).^2);
[X, fval, exitflag] = fminsearch(errfun, X0, options)
Kd = X(1)
Pt = X(2)
a = X(3)
b = X(4)
c = X(5)
d = X(6)

subplot(2,2,1);
plot(ligand_concentration, x(:,1:5), '.');
hold on;
%plot(ligand_concentration, model(Pt,ligand_concentration,a,b,Kd),'b-');
plot(ligand_concentration, model(Kd,Pt,a,b,c,d),'b-');
legend('1 uM', '100 nM', '10 nM', '1 nM', '0 nM');
title('erlotinib');
xlabel('ligand concentration (uM)');
ylabel('fluorescence 280/480 nm');

subplot(2,2,3);
semilogy(ligand_concentration, x(:,1:5), '.');
hold on;
%semilogy(ligand_concentration, model(Pt,ligand_concentration,a,b,Kd),'b-');
semilogy(ligand_concentration, model(Kd,Pt,a,b,c,d),'b-');
legend('1 uM', '100 nM', '10 nM', '1 nM', '0 nM');
title('erlotinib');
xlabel('ligand concentration (uM)');
ylabel('fluorescence 280/480 nm');

subplot(2,2,2);
plot(x(:,1:5), '.');
hold on;
%plot(model(Pt,ligand_concentration,a,b,Kd),'b-');
plot(model(Kd,Pt,a,b,c,d),'b-');
legend('1 uM', '100 nM', '10 nM', '1 nM', '0 nM');
title('erlotinib');
xlabel('ligand 1:2 dilution');
ylabel('fluorescence 280/480 nm');

subplot(2,2,4);
semilogy(x(:,1:5), '.');
hold on;
%semilogy(model(Pt,ligand_concentration,a,b,Kd),'b-');
semilogy(model(Kd,Pt,a,b,c,d),'b-');
legend('1 uM', '100 nM', '10 nM', '1 nM', '0 nM');
title('erlotinib');
xlabel('ligand 1:2 dilution');
ylabel('fluorescence 280/480 nm');

