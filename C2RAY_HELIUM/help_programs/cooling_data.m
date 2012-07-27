% This script loads the cooling data from the tables that are used in the
% code right now and plots them together with the cooling rates from
% Hui%Gnedin 1997 and Cen 1992
load('cooling_rates_our_tables.mat');
H1A=10.^H1A;
H1B=10.^H1B;
He1=10.^He1;
He2=10.^He2;
He0=10.^He0;
H0=10.^H0;


logT=1:0.01:9;
T=10.^logT;
% now, these ones are from Hui & Gnedin 1997 MNRAS 292---------------------
T_H(1)=157807; %H threshold in K
T_H(2)=285335; % He0 threshold in K
T_H(3)=631515; % He+ threshold in K
k_b=1.3806504e-16; % in erg/K

for i=1:3
lambda(i,:)=2*(T_H(i)./T);
end

% recombination (including dielectronic recombination) cooling from
% Hui& Gnedin
R_A(1,:)=1.778e-29.*T.*lambda(1,:).^1.965./(1+(lambda(1,:)/0.541).^0.502).^2.697;  % H case A ion
R_B(1,:)=3.435e-30.*T.*lambda(1,:).^1.970./(1+(lambda(1,:)/2.250).^0.376).^3.720; % H case B ion

R_A(2,:)=3.0e-14*lambda(2,:).^0.654*k_b.*T;           % He0 case A ion   (>5e3 K)
R_B(2,:)=1.26e-14*lambda(2,:).^0.750*k_b.*T;          % He0 case B ion   (>5e3 K)

D_2a=1.9e-3*T.^(-3/2).*exp(-0.75*lambda(3,:)/2).*(1+0.3*exp(-0.15*lambda(3,:)/2))*0.75*k_b*T_H(3);


R_A(3,:)=2*1.778e-29.*T.*lambda(3,:).^1.965./(1+(lambda(3,:)/0.541).^0.502).^2.697; % He1 case A ion 
R_B(3,:)=2*3.435e-30.*T.*lambda(3,:).^1.970./(1+(lambda(3,:)/2.250).^0.376).^3.720; % He1 case B ion
%        ^ in the paper, it says factor 8 here, correct however is a factor 2

% collisional ionization (cooling) from Hui & Gnedin
CI(1,:) =21.11*T.^(-3/2).*exp(-lambda(1,:)/2).*lambda(1,:).^-1.089./(1+(lambda(1,:)/0.354).^0.874).^1.101; % collisional H ionization
CR(1,:)=k_b*T_H(1)*CI(1,:);  % collisional H cooling
CI(2,:) =32.38*T.^(-3/2).*exp(-lambda(2,:)/2).*lambda(2,:).^-1.146./(1+(lambda(2,:)/0.416).^0.987).^1.056; % collisional He0 ionization
CR(2,:) =k_b*T_H(2)*CI(2,:); % collisional He0 cooling
CI(3,:) =19.95*T.^(-3/2).*exp(-lambda(3,:)/2).*lambda(3,:).^-1.089./(1+(lambda(3,:)/0.553).^0.735).^1.275; % collisional He1 ionization
CR(3,:) =k_b*T_H(3)*CI(3,:);% collisional He1 cooling

%--------------------------------------------------------------------------
% This part is not really about the cooling but about collisional
% ionization
% collisional ionization rates as included in the code:
colh0=5.835410275968903E-011  ;
colhe0=2.709585263842644E-011 ;
colhe1=5.707336249831606E-012 ;

    sqrtt0 =sqrt(T);
    acolh0 =colh0*sqrtt0.*exp(-T_H(1)./T);
    acolhe0=colhe0*sqrtt0.*exp(-T_H(2)./T);
    acolhe1=colhe1*sqrtt0.*exp(-T_H(3)./T);
    
% collisional ionization rates from Cen 1992:    
    Cen_H   = 5.85e-11*sqrt(T).*exp(-157809.1./T)./(1+sqrt(T/1e5));
    Cen_He0 = 2.38e-11*sqrt(T)./(1+sqrt(T/1e5)).*exp(-285335.4./T);
    Cen_He1 = 5.68e-12*sqrt(T)./(1+sqrt(T/1e5)).*exp(-631515./T);

    % Now, plotting also Janev Langer "Elemenrary Processes in
    % Hydrogen-Helium Plasmas"
    b=[-3.27139e1, 1.35365e1, ...
      -5.73932e0, 1.56315e0, ...
      -2.87705e-1, 3.48255e-2, ...
      -2.63197e-3, 1.11954e-4, ...
      -2.03914e-6];
  
   bHe=[-4.409864886e1 2.391596563e1   -1.07532302e1  3.05803875e0 ...
       -5.6851189e-1   6.79539123e-2   -5.0090561e-3   2.06723616e-4   -3.64916141e-6];
   
  
    
  units; % to get normal k_b back
  lnT=T*k_b/eV;lnT=log(lnT);
    Janev_H_ln=b(1)     +b(2)*lnT ...
              +b(3)*lnT.^2+b(4)*lnT.^3 ...
              +b(5)*lnT.^4+b(6)*lnT.^5 ...
              +b(7)*lnT.^6+b(8)*lnT.^7 ...
              +b(9)*lnT.^8;
           
     Janev_He_ln=bHe(1)     +bHe(2)*lnT ...
              +bHe(3)*lnT.^2+bHe(4)*lnT.^3 ...
              +bHe(5)*lnT.^4+bHe(6)*lnT.^5 ...
              +bHe(7)*lnT.^6+bHe(8)*lnT.^7 ...
              +bHe(9)*lnT.^8;         
     
    Janev_H=exp(Janev_H_ln);
    Janev_He=exp(Janev_He_ln);
    
    % From Dopita and Sutherland
    a_H_T=2.5e-10*(1+T/78945).^(-1).*sqrt(T).*exp(-T_H(1)./T);
    
    


%--------------------------------------------------------------------------
 test1=(CI(1,:)-Janev_H)./Janev_H;
% figure;semilogx(T,test1)

% collisional excitation cooling from Hui & Gnedin
EC(1,:)=7.5e-19*exp(-0.75*lambda(1,:)/2)./(1+sqrt(T./10^5));  %H0 excitation
EC(2,:)=9.1e-27*T.^(-0.1687).*exp(-13179./T);  %FROM BLACK, 1981

% Is this may be wrong? Do I need lambda(3,:) here? (was lambda(2,:))
EC(3,:)=5.54e-17*exp(-0.75*lambda(3,:)/2)./(1+sqrt(T./10^5))./(T.^0.397);  % He1 excitation


%to check, here some data from Cen, 1992:
H_coll_ion   = 1.27e-21*sqrt(T)./(1+sqrt(T/1e5)).*exp(-157809.1./T);
He0_coll_ion = 9.38e-22*sqrt(T)./(1+sqrt(T/1e5)).*exp(-285335.4./T);
He1_coll_ion = 4.95e-22*sqrt(T)./(1+sqrt(T/1e5)).*exp(-631515./T);

H_coll_ex    = 7.5e-19./(1+sqrt(T/1e5)).*exp(-118348./T);
He0_coll_ex  = 9.1e-27*T.^-0.1687./(1+sqrt(T/1e5)).*exp(-13179./T);
He1_coll_ex  = 5.54e-17*T.^-0.397./(1+sqrt(T/1e5)).*exp(-473638./T);



D_2b = 1.9e-3*T.^(-3/2).*exp(-4.7e5./T).*(1+0.3*exp(-9.4e4./T)); 
DR_2b= 0.75*k_b*T_H(3).*D_2b;


% The cooling data from Hummer 1994 is loaded here
recomrates_He1_Hummer; %--> beta_H_B, beta_H_1, beta_H_ff
                       %--> beta_He2_B, beta_He2_1, beta_He2_ff, TH, THe
recomrates_He0_Hummer_Storey %--> beta_He1_B, beta_He1_1,beta_He1_ff, T2

% make interpolation of the data
hummer_H_B=interp1(TH,beta_H_B,T);
hummer_H_A=interp1(TH,beta_H_1+beta_H_B,T);
hummer_H_ff=interp1(TH,beta_H_ff,T);
hummer_He2_B=interp1(THe,beta_He2_B,T);
hummer_He2_A=interp1(THe,beta_He2_B+beta_He2_1,T);
hummer_He2_ff=interp1(THe,beta_He2_ff,T);
hummer_He1_B=interp1(T2,beta_He1_B,T);
hummer_He1_A=interp1(T2,beta_He1_B+beta_He1_1,T);
hummer_He1_ff=interp1(T2,beta_He1_ff,T);



hummer_He1_A(isnan(hummer_He1_A))=0;
hummer_He1_B(isnan(hummer_He1_B))=0;
hummer_He1_ff(isnan(hummer_He1_ff))=0;

mycompilation=hummer_He1_B+hummer_He1_ff+DR_2b+CR(3,:)+EC(3,:);




% cooling: free-free      k--
%          recombination  A--/ B-
figure; title('H1')

loglog(H1A(:,1),H1A(:,2),'r--','Displayname','A, table')
title('H1')
hold on
loglog(H1B(:,1),H1B(:,2),'r-','Displayname','B, table')
loglog(TH,beta_H_ff+beta_H_B+beta_H_1,'g--','Displayname','A+ff, Hummer')
loglog(TH,beta_H_ff+beta_H_B,'g','Displayname','B+ff, Hummer')
loglog(T,R_A(1,:),'b--','LineWidth',1,'Displayname','A, Hui&Gn')
loglog(T,R_B(1,:),'b-','LineWidth',1,'Displayname','B, Hi&Gn')
loglog(TH,beta_H_B,'g-','LineWidth',1,'Displayname','B, Hummer')
loglog(TH,beta_H_B+beta_H_1,'g--','LineWidth',1,'Displayname','A, Hummer')
loglog(TH,beta_H_ff,'k--','LineWidth',1,'Displayname','ff, Hummer')
legend toggle


% cooling: collisional excitation   :
%          collisional ionization   -.
figure; title('H0')
loglog(H0(:,1),H0(:,2),'r','Displayname','table')
title('H0')
hold on
loglog(T,CR(1,:)+EC(1,:),'b','Displayname','CollExc+CollIon,Hui&Gn')
loglog(T,CR(1,:),'b-.','LineWidth',1,'Displayname','CollIon, Hui&Gn')
loglog(T,EC(1,:),'b:','LineWidth',1,'Displayname','CollExc, Hui&Gn')
legend toggle

% cooling: recombination   A--    B-
%          free-free       k--
figure; title('He2')
loglog(He2(:,1),He2(:,2),'r','Displayname','table')
title('He2')
hold on
loglog(THe,beta_He2_ff+beta_He2_B+beta_He2_1,'g--','Displayname','ff+A, Hummer')
loglog(THe,beta_He2_ff+beta_He2_B,'g','Displayname','ff+B, Hummer')
loglog(T,R_A(3,:),'b--','LineWidth',1,'Displayname','A, Hi&Gn'); 
plot(T,R_B(3,:),'b--','LineWidth',1,'Displayname','B, Hi&Gn')
loglog(THe,beta_He2_B,'g','LineWidth',1,'Displayname','B, Hummer')
loglog(THe,beta_He2_B+beta_He2_1,'g--','LineWidth',1,'Displayname','A, Hummer')
loglog(THe,beta_He2_ff,'k--','LineWidth',1,'Displayname','ff, Hummer')
legend toggle

% cooling: recombination
%          recombination dielectronic     k-
%          free-free                      k--
%          collisional ionization
%          collisional excitation
figure; title('He1')
loglog(He1(:,1),He1(:,2),'r','Displayname','table');
title('He1')
hold on
loglog(T,R_A(2,:)+CR(3,:)+DR_2b+EC(3,:),'b--','Displayname','A+dielectronic+CollIon+CollExc, Hui&Gn'); 
  plot(T,R_B(2,:)+CR(3,:)+DR_2b+EC(3,:),'b-','Displayname', 'B+dielectronic+CollIon+CollExc, Hui&Gn')
  plot(T,DR_2b,'k-','LineWidth',1,'Displayname','dielectronic1')
  plot(T,CR(3,:),'b-.','LineWidth',1,'Displayname','ColIon, Hui&Gn')
  plot(T,R_B(2,:),'-b','LineWidth',1,'Displayname','B, Hui&Gn')
  plot(T,R_A(2,:),'--b','LineWidth',1,'Displayname','A, Hui&Gn')
loglog(T2,beta_He1_B,'g','LineWidth',1,'Displayname','B, Hummer')
loglog(T2,beta_He1_B+beta_He1_1,'g--','LineWidth',1,'Displayname','A, Hummer')
loglog(T2,beta_He1_ff,'k--','LineWidth',1,'Displayname','ff, Hummer')
loglog(T,EC(3,:),'b:','LineWidth',1,'Displayname','ColExc, Hui&Gn')
newdata=EC(3,:)+DR_2b+CR(3,:);
legend toggle

part1=EC(3,:);%+CR(3,:)  ;
part2=beta_He1_B+beta_He1_ff;
part2=[part2,7e-25,8e-25,9e-25];
T2=[T2,4e4,6e4,9e4];
newtemp=10.^[1.00:0.01:9.0];
t_part1=interp1(T,part1,newtemp);
t_part2=interp1(T2,part2,newtemp);
for i=1:801;if isnan(t_part2(i)); t_part2(i)=0.0;end;end
total=t_part1+t_part2;
%figure;loglog(newtemp,total)

 
% cooling: collisional excitation (?)
%          collisional ionization
figure; title('He0');
loglog(He0(:,1),He0(:,2),'r','Displayname','table')
title('He0')
hold on
loglog(T,EC(2,:)+CR(2,:),'b','Displayname','ColExc+CollIon, Hui&Gn')
loglog(T,CR(2,:),'b-.','LineWidth',1,'Displayname','CollIon, Hui&Gn')
loglog(T,EC(2,:),'b:','LineWidth',1,'Displayname','ColExc, Black')

legend toggle

%--------------------------------------------------------------------------
