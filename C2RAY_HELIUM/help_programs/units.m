% alles in m, s, kg, J, K, C
c=2.99792458*10^8;                  % speed of light
cm= 0.01;                            
yr=3.156e7;                         % year
lyr=c*yr;                           % light year
AU=1.496e11;                        % astronomical unit
pc=3.262*lyr;                       % parsec
M_sol=1.989e30;                     % solar mass
L_sol=3.826e26;                     % solar luminosity
G=6.672e-11;                        % gravitational constant
h=0.7;                              % hubble parameter in 100 km /(s*Mpc)
rho_c=277.6*h^2*M_sol/(1000*pc)^3;  % critical density
k_b=1.38066e-23;                    % Boltzman constant
h_bar=6.626076e-34/2/pi;            % h bar             
m_p=1.672623e-27;                   % proton mass

e=1.602e-19;                        % Elementarladung
m_e=9.109e-31;                      % electron mass
eV=1.602e-19;                       % Elekronvolt
sigma_T=6.652e-25*cm^2;             % Thomson crossection
% To convert masses given in MeV to kg
MeV=m_e*c^2/0.511;                  % MeV in kg (m_e*c^2=0.511*MeV)
H_0=100*h*1000/(1e6*pc);             % Hubble constant

% to convert T_vir and v_circ,
% default values are for milkyway sized object
M=10^10*M_sol;
r_vir=200*1000*pc;
mu_m=2;


v_circ=sqrt(G*M/r_vir);
T_vir=mu_m*m_p*v_circ^2/(2*k_b);



% free fall time & cooling time
%n_g =                               % gas number density
%rho =                               % density
%Lambda=                             % ne t radiative cooling rate
%t_ff=3*pi/(32*G*rho)^0.5;
%t_cool=3*n_g*k_b*T/(2*Lambda);

