% Program to plot the frequency dependent cross sections of H He0 and He1
% in the three intervalls and compare them with the approximations made.

units %to have eV, h_bar ...
cross_sections; % to plot the cross sections according to Verner et al. 1996
plot_yes=1;

% begin of specification section ------------------------------------------

% freq steps
numsteps=150;
% how many subbins in the different frequency bins?

% numin2vec=1;                                       % 1
% numin2vec=[1,1.5];                                 % 2
% numin2vec=[1,1.3,1.7];                             % 3
%numin2vec=[1,1.15,1.3,1.5,1.7,1.9557];              % 6
%numin2vec=[1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9];  % 10
numin2vec=[1,1.02,1.05,1.07,1.1,1.15,1.2,1.25,1.3,1.35,1.4,...
    1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2,2.05,2.1,2.15];  % 26

%numin3vec=1;                                         % 1
%numin3vec=[1,2.0,4.0,10.0];                          % 4
% numin3vec=[1,1.5,2.0,3.0,4.0,7.0,10.0,20.0,50.0];   % 9
% numin3vec=[1,1.1,1.2,1.5,2,3,4,7,10,20,50];          % 11
% numin3vec=[1,1.05,1.1,1.2,1.4,1.7,2,3,5,7,10,15,20,30,50,70];          % 16
numin3vec=[1.0,1.05,1.1,...
           1.2,1.4,1.7,...
           2.0,2.5,3.0,...
           4.0,5.0,7.0,...
           10.0,15.0,20.0,...
           30.0,40.0,50.0,70.0,90.0];          % 20

% end of specification section--------------------------------------------- 


numin2=length(numin2vec);
numin3=length(numin3vec);




%-------------------------------\
nu_0_H=E_th(1)*eV/(h_bar*2*pi);    %\
nu_0_He0=E_th(2)*eV/(h_bar*2*pi);   %> thershold frequencies
nu_0_He1=E_th(3)*eV/(h_bar*2*pi);  %/
%-------------------------------/


% min freqs 
numin(1)=nu_0_H;
numin(2:numin2+1)=nu_0_He0*numin2vec;
numin(numin2+2:numin2+1+numin3)=nu_0_He1*numin3vec;

% Hpos, He0pos and He1pos are the positions of the threshold frequencies in
% the nu-vector.
a_H_th=sig(1,Hpos);
a_He0_th=sig(2,He0pos);
a_He1_th=sig(3,He1pos);



% max freqs
for i=1:numin2+numin3
numax(i)=numin(i+1);
end
numax(numin2+numin3+1)=100*nu_0_He1;


for i=1: numin2+numin3+1
    step(i)=(numax(i)-numin(i))/numsteps;
end

intm1=zeros(1,numin2+numin3+1);
intm2=intm1;intm3=intm2;pos=intm2;

for i=1:numin2+numin3+1
 [a,pos(i)]=min(abs(nu(1,:)-numin(i)));
 intm1(i)=sig(1,pos(i));
 intm2(i)=sig(2,pos(i));
 intm3(i)=sig(3,pos(i));
end
intm2(1)=0;
intm3(1:numin2+1)=0;

% set up frequency vectors
% the frequencies in the different frequency bins
nu=zeros(numin2+numin3+1,numsteps+1);
for i=1: numin2+numin3+1
nu(i,:)=numin(i):step(i):numax(i);
end



nulog=log10(nu);
siglog=log10(sig);
nu2log=log10(nu2);

% This is to calculate the functions f1 and f2 for ion and heat 

for i=1:length(numin)
    if (h_bar*2*pi*numin(i)/eV-E_th(1))>=11
        f2heat_h(i)=   (11/(h_bar*2*pi*numin(i)/eV-E_th(1))).^0.7;  
        f1heat_h(i)=   1;      
    else
        f2heat_h(i)=   0;  
        f1heat_h(i)=   0;
    end

    if (h_bar*2*pi*numin(i)/eV-E_th(2))>=11
        f2heat_he0(i)= (11/(h_bar*2*pi*numin(i)/eV-E_th(2))).^0.7; 
        f1heat_he0(i)=   1;    
    else    
        f2heat_he0(i)=   0;  
        f1heat_he0(i)=   0;
    end    

    if (h_bar*2*pi*numin(i)/eV-E_th(3))>=11
        f2heat_he1(i)= (11/(h_bar*2*pi*numin(i)/eV-E_th(3))).^0.7;
        f1heat_he1(i)=   1;          
    else    
        f2heat_he1(i)=   0;  
        f1heat_he1(i)=   0;
    end 
end



for i=1:length(numin)
    if (h_bar*2*pi*numin(i)/eV-E_th(1))>=28
        f2ion_h(i)=   (28/(h_bar*2*pi*numin(i)/eV-E_th(1))).^0.4;  
        f1ion_h(i)=   1;      
    else
        f2ion_h(i)=   0;  
        f1ion_h(i)=   0;
    end

    if (h_bar*2*pi*numin(i)/eV-E_th(2))>=28
        f2ion_he0(i)= (28/(h_bar*2*pi*numin(i)/eV-E_th(2))).^0.4; 
        f1ion_he0(i)=   1;    
    else    
        f2ion_he0(i)=   0;  
        f1ion_he0(i)=   0;
    end    

    if (h_bar*2*pi*numin(i)/eV-E_th(3))>=28
        f2ion_he1(i)= (28/(h_bar*2*pi*numin(i)/eV-E_th(3))).^0.4;
        f1ion_he1(i)=   1;          
    else    
        f2ion_he1(i)=   0;  
        f1ion_he1(i)=   0;
    end 
end

%for i=2:1+numin2; fprintf('%06.4f',f2heat_h(i)),fprintf('_dp, '); end; fprintf('\n')
%for i=2:1+numin2; fprintf('%06.4f',f2heat_he0(i)),fprintf('_dp, '); end; fprintf('\n')
%for i=2:1+numin2; fprintf('%06.4f',f2heat_he1(i)),fprintf('_dp, '); end; fprintf('\n')



% finding the power law fit coefficients in the sub-bins

for i=1:numin2+numin3;
trash=polyfit(nu2log(pos(i):pos(i+1)),siglog(1,(pos(i):pos(i+1))),1);
sh0(i)=abs(trash(1));
end
i=numin2+numin3+1;
trash=polyfit(nu2log(pos(i):end),siglog(1,(pos(i):end)),1);
sh0(i)=abs(trash(1));


for i=2:numin2+numin3;
trash=polyfit(nu2log(pos(i):pos(i+1)),siglog(2,(pos(i):pos(i+1))),1);
she0(i)=abs(trash(1));
end
i=numin2+numin3+1;
trash=polyfit(nu2log(pos(i):end),siglog(2,(pos(i):end)),1);
she0(i)=abs(trash(1));


for i=numin2+2:numin2+numin3;
trash=polyfit(nu2log(pos(i):pos(i+1)),siglog(3,(pos(i):pos(i+1))),1);
she1(i)=abs(trash(1));
end
i=numin2+numin3+1;
trash=polyfit(nu2log(pos(i):end),siglog(3,(pos(i):end)),1);
she1(i)=abs(trash(1));

              
 for i=1:numin2+numin3+1
 sig_H(i,:)   = intm1(i)*(nu(i,:)./numin(i)).^(-sh0(i));
 sig_He0(i,:) = intm2(i)*(nu(i,:)./numin(i)).^(-she0(i));
 sig_He1(i,:) = intm3(i)*(nu(i,:)./numin(i)).^(-she1(i));
 end
 
 
 
 for i=2:numin2+1
     sig_H_approx(i,:) = intm1(i)*(nu(i,:)./numin(i)).^(-she0(i));
 end
 
 for i=numin2+2:numin2+numin3+1
     sig_H_approx(i,:)=intm1(i)*(nu(i,:)./numin(i)).^(-she1(i));
     sig_He0_approx(i,:)=intm2(i)*(nu(i,:)./numin(i)).^(-she1(i));   
 end

sig_Hlog=log10(sig_H); 
sig_He0log=log10(sig_He0);
sig_He1log=log10(sig_He1);
sig_H_approxlog=log10(sig_H_approx);
sig_He0_approxlog=log10(sig_He0_approx);


figure1=figure;
subplot(2,2,1,'align')
i=1;
plot(nu2log(pos(i):pos(i+1)), siglog(1,pos(i):pos(i+1)),'b','LineWidth',1);hold on
plot(nulog(i,:),sig_Hlog(i,:),'k--','LineWidth',1); hold on
xlabel('\nu / Hz','FontSize', 18);
        ylabel('\sigma / m^2','FontSize', 18); title('bin 1')

        
      
subplot(2,2,2,'align')
for i=2:numin2+1
plot(nulog(i,:),sig_He0log(i,:),'k--','LineWidth',1); hold on
plot(nu2log(pos(i):pos(i+1)), siglog(2,pos(i):pos(i+1)),'r','LineWidth',1)

plot(nu2log(pos(i):pos(i+1)), siglog(1,pos(i):pos(i+1)),'b','LineWidth',1)
plot(nulog(i,:),sig_H_approxlog(i,:),'k--','LineWidth',1);

end
xlabel('\nu / Hz','FontSize', 18);
        ylabel('\sigma / m^2','FontSize', 18);title('bin 2')
               

subplot(2,2,3:4,'align')
for i=numin2+2:numin3+numin2
plot(nulog(i,:),sig_He1log(i,:),'k--','LineWidth',1); hold on
plot(nu2log(pos(i):pos(i+1)), siglog(3,pos(i):pos(i+1)),'g','LineWidth',1)

plot(nu2log(pos(i):pos(i+1)), siglog(2,pos(i):pos(i+1)),'r','LineWidth',1)
plot(nulog(i,:),sig_He0_approxlog(i,:),'k--','LineWidth',1);

plot(nu2log(pos(i):pos(i+1)), siglog(1,pos(i):pos(i+1)),'b','LineWidth',1)
plot(nulog(i,:),sig_H_approxlog(i,:),'k--','LineWidth',1);
end
i=numin2+numin3+1;
plot(nulog(i,:),sig_He1log(i,:),'k--','LineWidth',1);  hold on
plot(nulog(i,:),sig_He0_approxlog(i,:),'k--','LineWidth',1);
plot(nulog(i,:),sig_H_approxlog(i,:),'k--','LineWidth',1);

plot(nu2log(pos(i):end), siglog(3,pos(i):end),'g','LineWidth',1)
plot(nu2log(pos(i):end), siglog(2,pos(i):end),'r','LineWidth',1)
plot(nu2log(pos(i):end), siglog(1,pos(i):end),'b','LineWidth',1)
xlabel('\nu / Hz','FontSize', 18);
        ylabel('\sigma / m^2','FontSize', 18);title('bin 3')
        
 annotation(figure1,'textbox',...
    [0.136803095661 0.726882735473046 0.252482618624714 0.0809523809523829],...
    'String',{'blue=H','red= He0','green= He1','from Janev','black dashed: fits'},...
    'FontSize',11,...
    'FitBoxToText','off',...
    'LineStyle','none');       


        
  NumBndin2=numin2;      
        
fprintf('intm1(2:1+NumBndin2)=(/ & '); fprintf('\n');   for i=2:1+numin2; fprintf('%d',intm1(i)*1e4),fprintf('_dp, ');if i==5||i==10||i==15||i==20; fprintf('\n'); end;  end;fprintf('/) '); fprintf('\n'); fprintf('\n')
fprintf('intm2(2:1+NumBndin2)=(/ & '); fprintf('\n');   for i=2:1+numin2; fprintf('%d',intm2(i)*1e4),fprintf('_dp, ');if i==5||i==10||i==15||i==20; fprintf('\n'); end;  end;fprintf('/) '); fprintf('\n'); fprintf('\n')

fprintf('intm1(NumBndin2+2:NumBndin2+NumBndin3+1)=(/ & '); fprintf('\n');   for i=numin2+2:1+numin2+numin3; fprintf('%d',intm1(i)*1e4),fprintf('_dp, ');if i==NumBndin2+5||i==NumBndin2+10||i==NumBndin2+15||i==NumBndin2+20; fprintf('\n'); end;  end;fprintf('/) '); fprintf('\n'); fprintf('\n')
fprintf('intm2(NumBndin2+2:NumBndin2+NumBndin3+1)=(/ & '); fprintf('\n');   for i=numin2+2:1+numin2+numin3; fprintf('%d',intm2(i)*1e4),fprintf('_dp, ');if i==NumBndin2+5||i==NumBndin2+10||i==NumBndin2+15||i==NumBndin2+20; fprintf('\n'); end;  end;fprintf('/) '); fprintf('\n'); fprintf('\n')
fprintf('intm3(NumBndin2+2:NumBndin2+NumBndin3+1)=(/ & '); fprintf('\n');   for i=numin2+2:1+numin2+numin3; fprintf('%d',intm3(i)*1e4),fprintf('_dp, ');if i==NumBndin2+5||i==NumBndin2+10||i==NumBndin2+15||i==NumBndin2+20; fprintf('\n'); end;  end;fprintf('/) '); fprintf('\n'); fprintf('\n')


fprintf('sh0(2:1+NumBndin2)=(/ & '); fprintf('\n');   for i=2:1+numin2; fprintf('%06.4f',sh0(i)),fprintf('_dp, ');if i==5||i==10||i==15||i==20; fprintf('\n'); end; end; fprintf('\n') ; fprintf('\n')
fprintf('she0(2:1+NumBndin2)=(/ & '); fprintf('\n');  for i=2:1+numin2; fprintf('%06.4f',she0(i)),fprintf('_dp, ');if i==5||i==10||i==15||i==20; fprintf('\n'); end; end; fprintf('\n'); fprintf('\n')
fprintf('sh0(NumBndin2+2:NumBndin2+NumBndin3+1)=(/ & '); fprintf('\n');    for i=2+numin2:1+numin3+numin2; fprintf('%06.4f',sh0(i)),fprintf('_dp, ');if i==NumBndin2+5||i==NumBndin2+10||i==NumBndin2+15||i==NumBndin2+20; fprintf('\n'); end;  end;fprintf('/) '); fprintf('\n'); fprintf('\n')
fprintf('she0(NumBndin2+2:NumBndin2+NumBndin3+1)=(/ & '); fprintf('\n');   for i=2+numin2:1+numin3+numin2; fprintf('%06.4f',she0(i)),fprintf('_dp, ');if i==NumBndin2+5||i==NumBndin2+10||i==NumBndin2+15||i==NumBndin2+20; fprintf('\n'); end;  end;fprintf('/) '); fprintf('\n'); fprintf('\n')
fprintf('she1(NumBndin2+2:NumBndin2+NumBndin3+1)=(/ & '); fprintf('\n');   for i=2+numin2:1+numin3+numin2; fprintf('%06.4f',she1(i)),fprintf('_dp, ');if i==NumBndin2+5||i==NumBndin2+10||i==NumBndin2+15||i==NumBndin2+20; fprintf('\n'); end;  end;fprintf('/) '); fprintf('\n'); fprintf('\n')


fprintf('f2heat_h(NumBndin2+2:NumBndin2+NumBndin3+1)=(/ & '); fprintf('\n');   for i=2+numin2:1+numin3+numin2; fprintf('%06.4f',f2heat_h(i)),fprintf('_dp, ');if i==NumBndin2+5||i==NumBndin2+10||i==NumBndin2+15||i==NumBndin2+20; fprintf('\n'); end;end;fprintf('/) '); fprintf('\n'); fprintf('\n')
fprintf('f2heat_he0(NumBndin2+2:NumBndin2+NumBndin3+1)=(/ & '); fprintf('\n');  for i=2+numin2:1+numin3+numin2; fprintf('%06.4f',f2heat_he0(i)),fprintf('_dp, ');if i==NumBndin2+5||i==NumBndin2+10||i==NumBndin2+15||i==NumBndin2+20; fprintf('\n');end; end;fprintf('/) '); fprintf('\n'); fprintf('\n')
fprintf('f2heat_he1(NumBndin2+2:NumBndin2+NumBndin3+1)=(/ & '); fprintf('\n');  for i=2+numin2:1+numin3+numin2; fprintf('%06.4f',f2heat_he1(i)),fprintf('_dp, ');if i==NumBndin2+5||i==NumBndin2+10||i==NumBndin2+15||i==NumBndin2+20; fprintf('\n'); end;end;fprintf('/) '); fprintf('\n'); fprintf('\n')
   
fprintf('f2heat_h(2:1+NumBndin2)=(/ & '); fprintf('\n');   for i=2:1+numin2; fprintf('%06.4f',f2heat_h(i)),fprintf('_dp, ');if i==5||i==10||i==15||i==20; fprintf('\n'); end; end;fprintf('/) '); fprintf('\n'); fprintf('\n')
fprintf('f2heat_he0(2:1+NumBndin2)=(/ & '); fprintf('\n'); for i=2:1+numin2; fprintf('%06.4f',f2heat_he0(i)),fprintf('_dp, ');if i==5||i==10||i==15||i==20; fprintf('\n'); end; end;fprintf('/) '); fprintf('\n'); fprintf('\n')
fprintf('f2heat_he1(2:1+NumBndin2)=(/ & '); fprintf('\n'); for i=2:1+numin2; fprintf('%06.4f',f2heat_he1(i)),fprintf('_dp, ');if i==5||i==10||i==15||i==20; fprintf('\n'); end; end;fprintf('/) '); fprintf('\n'); fprintf('\n')



fprintf('f1heat_h(NumBndin2+2:NumBndin2+NumBndin3+1)=(/ & '); fprintf('\n');   for i=2+numin2:1+numin3+numin2; fprintf('%06.4f',f1heat_h(i)),fprintf('_dp, ');if i==NumBndin2+5||i==NumBndin2+10||i==NumBndin2+15||i==NumBndin2+20; fprintf('\n'); end;end;fprintf('/) '); fprintf('\n'); fprintf('\n')
fprintf('f1heat_he0(NumBndin2+2:NumBndin2+NumBndin3+1)=(/ & '); fprintf('\n');  for i=2+numin2:1+numin3+numin2; fprintf('%06.4f',f1heat_he0(i)),fprintf('_dp, ');if i==NumBndin2+5||i==NumBndin2+10||i==NumBndin2+15||i==NumBndin2+20; fprintf('\n');end; end;fprintf('/) '); fprintf('\n'); fprintf('\n')
fprintf('f1heat_he1(NumBndin2+2:NumBndin2+NumBndin3+1)=(/ & '); fprintf('\n');  for i=2+numin2:1+numin3+numin2; fprintf('%06.4f',f1heat_he1(i)),fprintf('_dp, ');if i==NumBndin2+5||i==NumBndin2+10||i==NumBndin2+15||i==NumBndin2+20; fprintf('\n'); end;end;fprintf('/) '); fprintf('\n'); fprintf('\n')
   
fprintf('f1heat_h(2:1+NumBndin2)=(/ & '); fprintf('\n');   for i=2:1+numin2; fprintf('%06.4f',f1heat_h(i)),fprintf('_dp, ');if i==5||i==10||i==15||i==20; fprintf('\n'); end; end;fprintf('/) '); fprintf('\n'); fprintf('\n')
fprintf('f1heat_he0(2:1+NumBndin2)=(/ & '); fprintf('\n'); for i=2:1+numin2; fprintf('%06.4f',f1heat_he0(i)),fprintf('_dp, ');if i==5||i==10||i==15||i==20; fprintf('\n'); end; end;fprintf('/) '); fprintf('\n'); fprintf('\n')
fprintf('f1heat_he1(2:1+NumBndin2)=(/ & '); fprintf('\n'); for i=2:1+numin2; fprintf('%06.4f',f1heat_he1(i)),fprintf('_dp, ');if i==5||i==10||i==15||i==20; fprintf('\n'); end; end;fprintf('/) '); fprintf('\n'); fprintf('\n')


fprintf('f2ion_h(NumBndin2+2:NumBndin2+NumBndin3+1)=(/ & '); fprintf('\n');   for i=2+numin2:1+numin3+numin2; fprintf('%06.4f',f2ion_h(i)),fprintf('_dp, ');if i==NumBndin2+5||i==NumBndin2+10||i==NumBndin2+15||i==NumBndin2+20; fprintf('\n'); end;end;fprintf('/) '); fprintf('\n'); fprintf('\n')
fprintf('f2ion_he0(NumBndin2+2:NumBndin2+NumBndin3+1)=(/ & '); fprintf('\n'); for i=2+numin2:1+numin3+numin2; fprintf('%06.4f',f2ion_he0(i)),fprintf('_dp, ');if i==NumBndin2+5||i==NumBndin2+10||i==NumBndin2+15||i==NumBndin2+20; fprintf('\n');end; end;fprintf('/) '); fprintf('\n'); fprintf('\n')
fprintf('f2ion_he1(NumBndin2+2:NumBndin2+NumBndin3+1)=(/ & '); fprintf('\n'); for i=2+numin2:1+numin3+numin2; fprintf('%06.4f',f2ion_he1(i)),fprintf('_dp, ');if i==NumBndin2+5||i==NumBndin2+10||i==NumBndin2+15||i==NumBndin2+20; fprintf('\n'); end;end;fprintf('/) '); fprintf('\n'); fprintf('\n')
   
fprintf('f2ion_h(2:1+NumBndin2)=(/ & '); fprintf('\n');   for i=2:1+numin2; fprintf('%06.4f',f2ion_h(i)),fprintf('_dp, ');if i==5||i==10||i==15||i==20; fprintf('\n'); end; end;fprintf('/) '); fprintf('\n'); fprintf('\n')
fprintf('f2ion_he0(2:1+NumBndin2)=(/ & '); fprintf('\n'); for i=2:1+numin2; fprintf('%06.4f',f2ion_he0(i)),fprintf('_dp, ');if i==5||i==10||i==15||i==20; fprintf('\n'); end; end;fprintf('/) '); fprintf('\n'); fprintf('\n')
fprintf('f2ion_he1(2:1+NumBndin2)=(/ & '); fprintf('\n'); for i=2:1+numin2; fprintf('%06.4f',f2ion_he1(i)),fprintf('_dp, ');if i==5||i==10||i==15||i==20; fprintf('\n'); end; end;fprintf('/) '); fprintf('\n'); fprintf('\n')



fprintf('f1ion_h(NumBndin2+2:NumBndin2+NumBndin3+1)=(/ & '); fprintf('\n');   for i=2+numin2:1+numin3+numin2; fprintf('%06.4f',f1ion_h(i)),fprintf('_dp, ');if i==NumBndin2+5||i==NumBndin2+10||i==NumBndin2+15||i==NumBndin2+20; fprintf('\n'); end;end;fprintf('/) '); fprintf('\n'); fprintf('\n')
fprintf('f1ion_he0(NumBndin2+2:NumBndin2+NumBndin3+1)=(/ & '); fprintf('\n'); for i=2+numin2:1+numin3+numin2; fprintf('%06.4f',f1ion_he0(i)),fprintf('_dp, ');if i==NumBndin2+5||i==NumBndin2+10||i==NumBndin2+15||i==NumBndin2+20; fprintf('\n');end; end;fprintf('/) '); fprintf('\n'); fprintf('\n')
fprintf('f1ion_he1(NumBndin2+2:NumBndin2+NumBndin3+1)=(/ & '); fprintf('\n'); for i=2+numin2:1+numin3+numin2; fprintf('%06.4f',f1ion_he1(i)),fprintf('_dp, ');if i==NumBndin2+5||i==NumBndin2+10||i==NumBndin2+15||i==NumBndin2+20; fprintf('\n'); end;end;fprintf('/) '); fprintf('\n'); fprintf('\n')
   
fprintf('f1ion_h(2:1+NumBndin2)=(/ & '); fprintf('\n');   for i=2:1+numin2; fprintf('%06.4f',f1ion_h(i)),fprintf('_dp, ');if i==5||i==10||i==15||i==20; fprintf('\n'); end; end;fprintf('/) '); fprintf('\n'); fprintf('\n')
fprintf('f1ion_he0(2:1+NumBndin2)=(/ & '); fprintf('\n'); for i=2:1+numin2; fprintf('%06.4f',f1ion_he0(i)),fprintf('_dp, ');if i==5||i==10||i==15||i==20; fprintf('\n'); end; end;fprintf('/) '); fprintf('\n'); fprintf('\n')
fprintf('f1ion_he1(2:1+NumBndin2)=(/ & '); fprintf('\n'); for i=2:1+numin2; fprintf('%06.4f',f1ion_he1(i)),fprintf('_dp, ');if i==5||i==10||i==15||i==20; fprintf('\n'); end; end;fprintf('/) '); fprintf('\n'); fprintf('\n')

;fprintf('/) ');



