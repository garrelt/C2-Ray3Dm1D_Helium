% This script is for plotting the ionization cross section fits according
% to Verner et al 1996

E=13.6:0.1:5440;
E2=E*eV;
Mb=1e-18*1e-4;  %Mb in m^2
nu=E2./(h_bar*2*pi);
nu2=nu;

E_th=[13.6 24.59 54.42];
E_max=[5e4 5e4 5e4];
E_0=[0.4298 13.61 1.720];
sig_0=[5.475e4 9.492e2 1.369e4]*Mb;
y_a=[3.288e1 1.469 3.288e1];
P=[2.963 3.188 2.963];
y_w=[0 2.039 0];
y_0=[0 4.434e-1 0];
y_1=[0 2.136 0];


for i=1:3
    xi=E./E_0(i)-y_0(i); % vector
    yi=sqrt(xi.^2+y_1(i)^2);  %vector
    Fi=((xi-1).^2+y_w(i)^2).*yi.^(0.5*P(i)-5.5).*(1+sqrt(yi/y_a(i))).^-P(i);
    sig(i,:)=sig_0(i)*Fi;
end

 [a,Hpos]=min(abs(E(1,:)-E_th(1)));
 [a,He0pos]=min(abs(E(1,:)-E_th(2)));
 [a,He1pos]=min(abs(E(1,:)-E_th(3)));
 
%figure;
%loglog(nu(Hpos:end),sig(1,Hpos:end))
%hold on;
%loglog(nu(He0pos:end),sig(2,He0pos:end),'r')
%loglog(nu(He1pos:end),sig(3,He1pos:end),'g')