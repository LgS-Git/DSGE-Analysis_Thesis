%Params
tau_S  = 0.93;

lam  = 0.4;  
delt  = 1;
gam  = 1.67;
alph  = 0.25;
phi  = 1;
psiP  = 372.8;
thetaP = 9;
Tp  = 1/(thetaP - 1);

YB=  		 0.897735;
NB=  		 0.866025;
wB=  		 0.777461;
dB=  		 0.33665;

CkSB = @(tau) (wB*NB - Tp*YB + (1-tau).*delt.*dB)./YB;

%Optim weights
lamS = @(tau) lam.*CkSB(tau);

Wy = @(tau) (1+phi)./(1-alph) + (lam-lamS(tau))./(1-lamS(tau)).^2 .* ((1-alph).*(1-(1-tau).*delt)./CkSB(tau)).^2 .* ((1+phi)./(1-alph)).^2;
Wpi = @(tau) psiP * (1 - (lam-lamS(tau))./(1-lamS(tau)).*(1-(1-tau).*delt./CkSB(tau)));
Wdc = @(tau) lamS(tau).*(1-lamS(tau));

%Ricardian weights
Wy_ric = @(tau) ones(size(tau)).*(1+phi)./(1-alph);
Wpi_ric = @(tau) ones(size(tau)).*psiP;
Wdc_ric = @(tau) ones(size(tau)).*0;

%Keynesian weights
Wy_key = @(tau) (1+phi)/(1-alph) + 1./(1-CkSB(tau)) .* ((1-alph).*((1-(1-tau).*delt)./CkSB(tau))).^2 .* ((1+phi)/(1-alph)).^2;
%Wpi_key = @(tau) psiP * (1 - ((1-(1-tau).*delt)./CkSB(tau)));
Wpi_key = @(tau) psiP * (1 - (1-CkSB(tau))./(1-CkSB(tau)).*(1-(1-tau).*delt./CkSB(tau)));
Wdc_key = @(tau) CkSB(tau).*(1-CkSB(tau));

%percent diffs for standard calibration
WyDiff_ric = (Wy_ric(tau_S) - Wy(tau_S))/Wy(tau_S);
WpiDiff_ric = (Wpi_ric(tau_S) - Wpi(tau_S))/Wpi(tau_S);
WdcDiff_ric = (Wdc_ric(tau_S) - Wdc(tau_S))/Wdc(tau_S);

WyDiff_key = (Wy_key(tau_S) - Wy(tau_S))/Wy(tau_S);
WpiDiff_key = (Wpi_key(tau_S) - Wpi(tau_S))/Wpi(tau_S);
WdcDiff_key = (Wdc_key(tau_S) - Wdc(tau_S))/Wdc(tau_S);


%Plots
x0=10;
y0=10;
width = 400;
height = 400;

figure;
set(gcf,'position',[x0,y0,width,height])
grid on
hold on
fplot(Wpi);
fplot(Wpi_ric);
fplot(Wpi_key);
hold off
xlabel('\tau');
xlim([0 1]);
xline(0.93,"--");
legend(["optim","ric","key"],'FontSize',10);

figure;
set(gcf,'position',[x0,y0,width,height])
grid on
hold on
fplot(Wy);
fplot(Wy_ric);
fplot(Wy_key);
hold off
xlabel('\tau');
xlim([0 1]);
xline(0.93,"--");
legend(["optim","ric","key"],'FontSize',10);

figure;
set(gcf,'position',[x0,y0,width,height])
grid on
hold on
fplot(Wdc);
fplot(Wdc_ric);
fplot(Wdc_key);
hold off
xlabel('\tau');
xlim([0 1]);
xline(0.93,"--");
legend(["optim","ric","key"],'FontSize',10);

