% Chemotaxis model v 5.6.7.1 -- Dr + Flick + Mechanosensing + Chemochinesis
% implemented

% FM: to do:
% 1) Address \tau discrepancy with \nabla \neq 0
% 2) Address Pf shift with \nabla \neq 0
% 3) Implement dynamical system for chemokinesis

clc; clear; close all;
tic
% load sys
% load('C:\Users\Kwangmin\Dropbox\Mfiles\Project_Chemokinesis\Modeling\sys.mat');

CKonoff=1;   % 1 (ON) 0 (OFF) 2 (LTI approx of chemokinesis) - Chemokinesis switch
SDonoff=0;   % 1 (ON) 0 (OFF) 2 (Uniform concentration) - Steady concentration switch
MSonoff=1;   % 1 (ON) 0 (OFF) - Mechanosensing switch
FLonoff=1;   % 1 (Speed dependent flick) 2 (Always flick) 0 (Never flick) - Flick switch
DRonoff=1;   % 1 (ON) 0 (OFF) - Rotational diffusion switch

NumBin=5; % Number of speed bins for CMC plotting

showPlot=0;
directories = cd;
%   load([directories,'\DiffusionSolution_mid_14percent.mat']);  % microinjector
%   load([directories,'\DiffusionSolution_mid_33percent.mat']);  % 1mm transient channel
% load('C:\Users\Kwangmin\Dropbox\Mfiles\Project_Chemokinesis\DiffusionSolution_SteadyGaussian_time_100s.mat')
% load('C:\Users\Kwangmin\Dropbox\Mfiles\Project_Chemokinesis\DiffusionSolution_SteadyGaussian_time_100s_2D.mat')
load('C:\MIT\Data\DiffusionSolution_SteadyGaussian_time_10s_2D.mat')

% Chemokinesis switch
% vini=63.8284e-6; % From the experimental data, Consant swimming speed [m/s]   
% vmax=81.1722e-6;
C=0.25; %KS: used if this is a uniform concentration
cthresh=0.05; %KS: Threshold is set to 50nM. Here in the setup, all the cells will undergo chemokinesis (100nM is lowest) 
% v=vini;

% Constants
t0=1/2.0; % unbiased run length [s] %KS: this is used if mechanosensing is off (***)
tm=0.1; % Memory decay time scale [s]
dt=0.1; % interval between time steps [s]
alpha =30; % Constants the exponential bias in run length [s]
% nt=2250/dt;
nt=2000/dt;
t=dt:dt:(nt*dt);

% Assume the concentration (in arbitrary units) is given by the following
% equation C = C(x) = m + mkx. Here [x,y,z] at the origin (in meters),
% where the bacteria are released. K sets the length scale over which 
% concentration doubles;

L=3000e-6; %FM: channel length in m

% m=0.1; %KS: concentration in uM (100nM=0.1uM at the left edge)
% k=500;  % sets the gradient 
m=0.4; %KS: concentration in uM (100nM=0.1uM at the right edge)
k=-500;  % sets the gradient 
% xmaximum=.6e-3; %KS: this is in meters (=600um)
xmaximum=3e-3; %KS: this is in meters (=3mm)
ymaximum=3e-3; %KS: this is in meters (=3mm)
kd=10;   % micromolar (uM) %KS: this can change

Dr=0.0349; % rad^2/s (based on geometry from V.alginolyticus)
% dircos=cos(sqrt(2*Dr*dt));
RotDr=sqrt(4*Dr*dt); %KS: del_Theta

% Now run the model
% number of bacteria
nb=3000; 

% Initialize matrices
% Assume the bacteria starts off at the origin. It knows the local 
% concentration but dcdc=0;

x=zeros(nt,nb); y=zeros(nt,nb); %KS: Changed for 2D Gaussian
% x=zeros(nt,nb);
x(1,:)=rand(1,nb).*L; 
y(1,:)=rand(1,nb).*L; %KS: Changed for 2D Gaussian
% y=zeros(nt,nb); %FM: new initialization % First option matches better for 3-way channel

speeds=zeros(nt,nb);

avgspeed=zeros(1,nb); %average speed of each bacterium
avgx=zeros(nt,1); % average bacteria position at each time step
CMC=zeros(nt,1);

c=zeros(nt,nb);
tau=t0*ones(nt,nb);
dcdt=zeros(nt,nb);
reverse=zeros(nt,nb);
flick=zeros(nt,nb);
reorientation=zeros(nt,nb);
Nout=0; % for counting cells going out of the domain in x

for i=1:nb
    disp([num2str(i),'/',num2str(nb)]);
%     h=waitbar(i/nb);
%     vini=gamrnd(5.71149,4.94949)*1E-6; %KS: this can be modified for various speed at the population level
    vini=rand*60*1E-6; %KS: this can be modified for various speed at the population level
%     vini=60*1E-6; %KS: this can be modified for various speed at the population level
    v=vini;
%     c(1,i)=k*x(1,i)+m;
    if round(x(1,i)*1e6)==0 && round(y(1,i)*1e6)> 0 %KS: Changed for 2D Gaussian
     c(1,i)=solution(1,round(y(1,i)*1e6));  
    elseif round(y(1,i)*1e6)==0 && round(x(1,i)*1e6)> 0
     c(1,i)=solution(round(x(1,i)*1e6),1);       
    elseif round(x(1,i)*1e6)==0 && round(y(1,i)*1e6)==0
     c(1,i)=solution(1,1);
    else
     c(1,i)=solution(round(x(1,i)*1e6),round(y(1,i)*1e6));  
    end
    ven=0;
    speeds(1,i)=v;
    
    % choose random initial direction
     ax=(rand-0.5); ay=(rand-0.5); %KS: z-direction used here?
     dx=ax*v*dt/sqrt(ax^2+ay^2);
     dy=ay*v*dt/sqrt(ax^2+ay^2); 
     
    %KS: choose either the flick (fr=1) or reversal(fr=0) as a 1st reorientation
     if rand>0.5
      fr=1; %KS: odd number is flick cycle
     else
      fr=0; %KS: even number is reversal cycle
     end
     
    % Start at time step 2
     tmpx=x(1,i)+dx;
     tmpy=y(1,i)+dy;
     
     if showPlot==1
         figure;
         scatter(x(1,i)+dx,y(1,i)+dy);
         hold on
     end
    % axis([tmpx-tmpx*10 tmpx+tmpx*10 tmpy-tmpy*10 tmpy+tmpy*10])
    
    %KS: Start the time loop for this bacteria
%         h=waitbar(0,['Vibrio-model running: ',num2str(i),'/',num2str(nb)]);
       for j=2:nt;
%             waitbar(j/(nt-1),h);
            %KS: Chemokinesis switch (On/Off), For On : two options (stepwise or linear)
            if CKonoff==1 
                if c(j-1,i)>=cthresh  % chemokinesis                   
                    v=vini*1.3;  % step-wise 30% enhancement
%                   v=vini+(vmax-vini);  % step-wise
%                   v=vini+(vmax-vini)*(c(j,i)-cthresh);   % linear increase
                end    
            elseif CKonoff==0
              % just take the initial value as its speed      
            elseif CKonoff==2 %FM: added dynamic chemokinesis
                if c(j-1,i)>=cthresh
                    venTmp=lsim(sys,ones(1,100),linspace(0,dt,100),ven./0.04549900578);
                    ven=mean(venTmp);
                    venT(j,i)=ven;
                    v=vini*(1+ven);
                    ven=venTmp(end);
                    clear venTmp
                else
                    venTmp=lsim(sys,zeros(1,100),linspace(0,dt,100),ven./0.04549900578);
                    ven=venTmp(end);
                    venT(j,i)=ven;
                    v=vini*(1+ven);
                    clear venTmp
                end
            end 
            
          %KS: Cell is sense and moving  
          %KS: Either the cell is (1) going on the original direction or (2)New direction is chosen after reorientation
          %KS: New reorientational direction chosen at the last time step 
          %KS: Add 2-D rotational diffusion
          speeds(j,i)=v;
             AA1=ax/sqrt(ax^2+ay^2);  %cos(theta)
             AA2=ay/sqrt(ax^2+ay^2);  %sin(theta)
             OD=atan2(AA2,AA1); %KS: original direction  
              if rand>0.5
               ND=OD+RotDr; %KS: NewDirection = Theta + delTheta
              else
               ND=OD-RotDr; %KS: NewDirection = Theta - delTheta
              end
              
          %KS: Rotational diffusion switch (On/Off):
          if DRonoff==1
            dx=v*dt*cos(ND); %new direction
            dy=v*dt*sin(ND);  
          elseif DRonoff==0  
            dx=v*dt*cos(OD); %old direction
            dy=v*dt*sin(OD); 
          end
          x(j,i)=x(j-1,i)+dx;
          y(j,i)=y(j-1,i)+dy;
          clear AA1 AA2 OD ND
             
            %KS: Steady-concentration switch (On/Off)
            if SDonoff==1      %Linear steady concentration
                  c(j,i)=k*x(j,i)+m;     
            elseif SDonoff==2  %Uniform concenration
                  c(j,i)=C;      
            elseif SDonoff==0  %Transient concentration
                 a1=round(j*dt);   % time
                 a2=round(x(j,i)*1e6);  
                 a3=round(y(j,i)*1e6);   
                 if a1==0; a1=a1+1; end
                 if a2<=0; a2=1; end
                 if a2>3000; a2=3000; end
                 if a3<=0; a3=1; end %KS: Changed for 2D Gaussian
                 if a3>3000; a3=3000; end %KS: Changed for 2D Gaussian
%                  c(j,i)=solution(a1,a2);  % grab data from the matrix, every um, every second
%                  c(j,i)=solution(1,a2);  %KS: Steady 1D Gaussian
                 c(j,i)=solution(a2,a3);   %KS: Changed for 2D Gaussian
                 clear('a1'); clear('a2'); clear('a3');
            end
            dcdt(j,i)=(c(j,i)-c(j-1,i))/dt;
             
           %KS: Setting new Run time integration (chemotaxis + mechanosensing)
            term1=kd*((kd+c(1:j,i)).^(-2));
            term2=(dcdt(1:j,i));
            term3=exp((t(1:j)-t(j))/tm)';
            dpbdt=sum((term1.*term2.*term3)*dt/tm);
            
             %KS: Mechanosensing switch (On/Off): mechanosensing modulates unbiased run time w.r.t. swimming speed
             if MSonoff==1 
%              %FM: mechanosensing
%               tau(j,i)=-0.3942*1/(1+exp(-0.2019*((v/1E-6)-18.88)))+0.8452;
              % FM: integration chemotaxis + mechanosensing
              t0=-0.3942*1/(1+exp(-0.2019*((v/1E-6)-18.88)))+0.8452;
             elseif MSonoff==0 
               % KS:just take the initial t0 value as its run time  
             end    
           tau(j,i)=t0*exp(alpha*dpbdt);
            
           %KS: Flick switch (On/Off/Always): 3 options
           if FLonoff==1 
            Pf=0.7069*1/(1+exp(-0.2642*((v/1E-6)-36.13)))+0.06261;  
           elseif FLonoff==2
            Pf=1; %Always flick without any speed dependence, run-reverse-flicker
           elseif FLonoff==0 
            Pf=0; %Never flick, This is just run-reverser  
           end    
             
             if ((dt/tau(j,i))> rand) && mod(fr,2)==1 %KS: odd number --> flick cycle 
                 reorientation(j,i)=1; %KS: for counts
                 flickrnd=rand;
                 if Pf>flickrnd %KS:it's flicking then choose (+/-) pi/2 ?
                  flick(j,i)=1; %KS: for counts
                     if rand>0.5
                         axtmp=ax;
                         aytmp=ay;
                         ax=aytmp;
                         ay=-axtmp;
                         fr=fr+1; %update this number for the next cycle
%                          dxtmp=dx;
%                          dytmp=dy;
%                          dx=dytmp;
%                          dy=-dxtmp;
                     else
                         axtmp=ax;
                         aytmp=ay;
                         ax=-aytmp;
                         ay=axtmp;
                         fr=fr+1; %update this number for the next cycle
%                          dxtmp=dx;
%                          dytmp=dy;
%                          dx=-dytmp;
%                          dy=dxtmp;
                     end
                 else %it's not flicking, just reversing
                    reverse(j,i)=1; %KS: for counts
                    ax=-ax;
                    ay=-ay;  
                    fr=fr+1; %update this number for the next cycle
%                      dxtmp=dx;
%                      dytmp=dy;
%                      dx=-dxtmp;
%                      dy=-dytmp;
                 end%KS:end of flick cycle
             elseif ((dt/tau(j,i))> rand) && mod(fr,2)==0 %KS: even number --> reversal cycle    
                    reorientation(j,i)=1; %KS: for counts
                    reverse(j,i)=1; %KS: for counts
                    ax=-ax;
                    ay=-ay;  
                    fr=fr+1; %update this number for the next cycle
             end%KS:end of reorientation cycle
                 
                    %Tumble -- new direction
                         %ax=(rand-0.5); ay=(rand-0.5); az=(rand-0.5);
                         %dx=ax*v*dt/sqrt(ax^2+ay^2+az^2);
                         %dy=ay*v*dt/sqrt(ax^2+ay^2+az^2);
                         %dz=az*v*dt/sqrt(ax^2+ay^2+az^2);
                         
%           %KS: Either the cell is (1) going on the original direction or (2)New direction is chosen after reorientation (add brownian diffusion to run)
%           %KS: Move at this time-step even after reorientation to make reorientation instantaneous 
%           %KS: Add 2-D rotational diffusion
%           speeds(j,i)=v;
%              AA1=ax/sqrt(ax^2+ay^2);  %cos(theta)
%              AA2=ay/sqrt(ax^2+ay^2);  %sin(theta)
%              OD=atan2(AA2,AA1); %KS: original direction  
%               if rand>0.5
%                ND=OD+RotDr; %KS: NewDirection = Theta + delTheta
%               else
%                ND=OD-RotDr; %KS: NewDirection = Theta - delTheta
%               end
%               
%           %KS: Rotational diffusion switch (On/Off):
%           if DRonoff==1
%             dx=v*dt*cos(ND); %new direction
%             dy=v*dt*sin(ND);  
%           elseif DRonoff==0  
%             dx=v*dt*cos(OD); %old direction
%             dy=v*dt*sin(OD); 
%           end
%           x(j,i)=x(j-1,i)+dx;
%           y(j,i)=y(j-1,i)+dy;
%           clear AA1 AA2 OD ND
          
        % Retake the step if swims out the channel 
         while x(j,i)<0 || x(j,i)>xmaximum || y(j,i)<0 || y(j,i)>ymaximum
           %Tumble -- new direction
            ax=(rand-0.5); 
            ay=(rand-0.5);         
            dx=ax*v*dt/sqrt(ax^2+ay^2);
            dy=ay*v*dt/sqrt(ax^2+ay^2);
            x(j,i)=x(j-1,i)+dx;
            y(j,i)=y(j-1,i)+dy;
            Nout=Nout+1;
         end
           
         if showPlot==1
             plot(x(j-1:j),y(j-1:j),'-o')
             if j>5
                 scatter(x(j-4:j),y(j-4:j),50,[0 0 1;0.25 0 0.75;0.5 0 0.5;0.75 0 0.25;1 0 0],'filled')
                 title(strcat('v= ',num2str(speeds(j,i)),', P_f=',num2str(Pf)));
                 axis([x(j,i)-100e-6 x(j,i)+100e-6 y(j,i)-100e-6 y(j,i)+100e-6])
             end
             pause(0.02)
         end            
       end%j(end of time-loop)
  
%   if mod(i,1000)==0  % to save files during computation    
%   save(['C:\Data\Kwangmin\Model_CaseB_Tm01_a30_CKoff_01111_',num2str(i),'.mat']);     
%   end     
% close(h);        
end  % i
% close(h);

%%
if showPlot==1
    hold off
end

colors{1}= rgb('red');
colors{2}= rgb('green');
colors{3}= rgb('blue');
colors{4}= rgb('magenta');
colors{5}= rgb('cyan');
colors{6}= rgb('gold');
colors{7}= rgb('lime');
colors{8}= rgb('dodgerblue');
colors{9}= rgb('pink');
colors{10}= rgb('black');

% save('D:\Kwangmin\Data\130917_Chemokinesis_analysis\Modeling\Model_CaseB_Tm01_a30_CKoff_01111.mat');

figure(1);
plot(x(:,1),y(:,1),'ro-','MarkerFaceColor','r','MarkerSize',3); hold on;
plot(x(1,1),y(1,1),'ks-','MarkerFaceColor','k','MarkerSize',3);
plot(x(:,2),y(:,2),'go-','MarkerFaceColor','g','MarkerSize',3);
plot(x(1,2),y(1,2),'ks-','MarkerFaceColor','k','MarkerSize',3);
plot(x(:,3),y(:,3),'bo-','MarkerFaceColor','b','MarkerSize',3);
plot(x(1,3),y(1,3),'ks-','MarkerFaceColor','k','MarkerSize',3);
plot(x(:,4),y(:,4),'mo-','MarkerFaceColor','m','MarkerSize',3);
plot(x(1,4),y(1,4),'ks-','MarkerFaceColor','k','MarkerSize',3);
plot(x(:,5),y(:,5),'co-','MarkerFaceColor','c','MarkerSize',3);
plot(x(1,5),y(1,5),'ks-','MarkerFaceColor','k','MarkerSize',3);
xlim([0 30e-4]);
ylim([0 30e-4]);
daspect([1 1 1]);

% Summarizing result
avgx=mean(x,2); % Average bacteria location at each time step (um/s)
CMC=-(avgx-(xmaximum/2))./(xmaximum/2); % 300um is the half width of the channel
avgx=avgx*(10^6); % Average bacteria location at each time step (um/s)

fmodel=sum(reorientation,1)./(nt*dt); %KS: reorientation frequency f for each cell (nt*dt=time in sec)
fFmodel=sum(flick,1)./(nt*dt); %KS: flicking frequency f for each cell (nt*dt=time in sec)
pFmodel=(sum(flick,1)./sum(reorientation,1))*2; %KS: assign speed for each bacterium

avgspeed=mean(speeds,1)*(10^6); %average speed of each bacterium

% PrTheo=(-0.3942*1./(1.+exp(-0.2019.*(speeds./1E-6-18.88)))+0.8452).^(-1);
% PfTheo=(0.7069*1./(1+exp(-0.2642.*((speeds/1E-6)-36.13)))+0.06261);

figure(2); % CMC plot
plot(t',CMC,'b.-');
ylim([-0.1 1]);
xlabel('Time (s)');
ylabel('CMC');

figure(3); % speed distribution (um/s)
hist(avgspeed,20); xlim([0 90]);
xlabel('Swimming speed');
ylabel('Count');

figure(4); % f: reorientation frequency
for i=1:nb
 plot(avgspeed(1,i),fmodel(1,i),'bo-'); hold on;
end
xlim([0 80]); %ylim([0 5.5]); 
xlabel('Swimming speed');
ylabel('reorientation frequency');

figure(5); % pF: probability of flick
for i=1:nb
 plot(avgspeed(1,i),pFmodel(1,i),'bo-');  hold on;
end
xlim([0 80]); ylim([0 1]);
xlabel('Swimming speed');
ylabel('Probability of flicking');

toc

% Bin equal number of cells in each speed bin
% avgspeed --> already in um/s
NumBin=5;
SpeedList=sort(avgspeed,'ascend');
Div=round(length(SpeedList)/NumBin);

idx=zeros(NumBin+1,1);
SpdDiv=zeros(NumBin+1,1);
idxspeed=zeros(NumBin,Div);
meanspdbin=zeros(NumBin,1);
for i=1:NumBin
idx(i,1)=(Div*(i-1))+1;
SpdDiv(i,1)=SpeedList(1,idx(i,1));
end
idx(NumBin+1,1)=length(SpeedList);
SpdDiv(end,1)=SpeedList(1,idx(NumBin+1,1));

for i=1:NumBin
 if i<NumBin
   idxspeed(i,:)=find(avgspeed(1,:)>=SpdDiv(i,1) & avgspeed(1,:)<SpdDiv(i+1,1)); %index for speed binning
 else 
   idxspeed(i,:)=find(avgspeed(1,:)>=SpdDiv(i,1) & avgspeed(1,:)<=SpdDiv(i+1,1)); %index for speed binning
 end
end

figure(6); % CMC plot
for i=1:NumBin
 meanspdbin(i,1)=mean(avgspeed(:,idxspeed(i,:)));
 avgxspd=mean(x(:,idxspeed(i,:)),2); % Average bacteria location at each time step (um/s)
 CMCspd=-(avgxspd-(xmaximum/2))./(xmaximum/2); % 300um is the half width of the channel
 plot(t',CMCspd,'-','color',colors{i}); hold on;
 clear avgxspd CMCspd
end
% xlim([0 500]);
ylim([-0.1 1.0]);
xlabel('Time (s)');
ylabel('CMC');

figure(7); 
for i=1:NumBin
 plot(mean(avgspeed(:,idxspeed(i,:))),mean(mean(x(end-10:end,idxspeed(i,:)))),'bo','color',colors{i}); hold on;
end
ylim([0 xmaximum]);
xlabel('Speed');
ylabel('Mean X location');

figure(8); % Speed vs CMC plot
tinterval=10; % final time points for averaging in seconds
for i=1:NumBin
 meanspdbin(i,1)=mean(avgspeed(:,idxspeed(i,:)));
 avgxspd=mean(x(:,idxspeed(i,:)),2); % Average bacteria location at each time step (um/s)
 CMCspd=-(avgxspd-(xmaximum/2))./(xmaximum/2); % 300um is the half width of the channel
 plot(mean(avgspeed(:,idxspeed(i,:))),mean(CMCspd(end-(tinterval/dt):end,1)),'o','color',colors{i},'MarkerFaceColor',colors{i},'MarkerSize',10); hold on;
%  plot(t',CMCspd,'-','color',colors{i}); hold on;
 clear avgxspd CMCspd
end
xlim([0 65]);
ylim([0 1.0]);
xlabel('Speed');
ylabel('Steadystate CMC');

% for i=1:NumBin
% %  figure,
% %  hist(x(end,idxspeed(i,:)),15);
%  figure,
%  plot(x(end-20:end,idxspeed(i,:)),y(end-20:end,idxspeed(i,:)),'bo');
%  xlim([0 xmaximum]);
%  xlabel('X location');
% end
% 
for i=1:NumBin
 figure,
 hist(mean(x(end-20:end,idxspeed(i,:))),15); hold on;
 xlim([0 xmaximum]);
 xlabel('X location');
end


