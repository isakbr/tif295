%% Same sized particles
clear all

clc
clf
load 120nm_ag_100_meas_19.txt
load lamp_spectrum.txt
load 120nm_pd_100_meas.txt
load 120nm_au_100_meas.txt
load start_100_19particles_nm_ag_100_meas.txt
load 120nm_pd_120_extrasession.txt

%Au
[lambda_au,particles,background,lamp]= getData(X120nm_au_100_meas,lamp_spectrum,19);
norm_particles_au=(particles-background)./lamp;

%Pd
[lambda_pd,particles,background,lamp]= getData(X120nm_pd_100_meas,lamp_spectrum,19);
norm_particles_pd=(particles-background)./lamp;

%Ag
[lambda_ag,particles,background,lamp]= getData(X120nm_ag_100_meas_19,lamp_spectrum,19);
norm_particles_ag=(particles-background)./lamp;

%Pd from extra session
[lambda_pd_ex,particles,background,lamp]= getData(X120nm_pd_120_extrasession,lamp_spectrum,19);
norm_particles_pd_ex=(particles-background)./lamp;


y_ax=ones(length(norm_particles_ag),1);

%Silver
figure(2)
for i=1:length(norm_particles_ag(1,:))
     plot3(lambda_ag,y_ax*i,norm_particles_ag(:,i),'color',[0 0.4470 0.7410],'Linewidth',2);
     hold on
     grid on
end
grid on
set(gca,'FontSize',28)
yl=ylabel('Particle','Interpreter','latex');
zl=zlabel('$I_{\textrm{norm}}$','Interpreter','latex');
xl=xlabel('$\lambda$ [nm]','Interpreter','latex');
xl.FontSize=34;
yl.FontSize=34;
zl.FontSize=34;


%Gold
figure(3)
for i=1:length(norm_particles_au(1,:))
     plot3(lambda_au,y_ax*i,norm_particles_au(:,i),'color',[0.8500 0.3250 0.0980],'Linewidth',2);
     hold on
     grid on
end
grid on
set(gca,'FontSize',28)
yl=ylabel('Particle','Interpreter','latex');
zl=zlabel('$I_{\textrm{norm}}$','Interpreter','latex');
xl=xlabel('$\lambda$ [nm]','Interpreter','latex');
xl.FontSize=34;
yl.FontSize=34;
zl.FontSize=34;


%Palladium
figure(4)
for i=1:length(norm_particles_pd(1,:))
     plot3(lambda_pd,y_ax*i,norm_particles_pd(:,i),'color',[0.4660 0.6740 0.1880],'Linewidth',2);
     hold on
     grid on
end
grid on
set(gca,'FontSize',28)
yl=ylabel('Particle','Interpreter','latex');
zl=zlabel('$I_{\textrm{norm}}$','Interpreter','latex');
xl=xlabel('$\lambda$ [nm]','Interpreter','latex');
xl.FontSize=34;
yl.FontSize=34;
zl.FontSize=34;

% 
% y_ax_ex=ones(length(norm_particles_ag),1);
% %Palladium extra session
% figure(5)
% for i=1:10:length(norm_particles_pd_ex(1,:))
%      plot3(lambda_pd_ex,y_ax_ex*i,norm_particles_pd_ex(:,i),'color',[0.4660 0.6740 0.1880],'Linewidth',2);
%      legend('Pd')
%      hold on
%      grid on
% end

%%
clc
clf
part_nr=19;
meas_nr=1;

peak_lam_pd=zeros(meas_nr,part_nr);
peak_lam_au=zeros(meas_nr,part_nr);
peak_lam_ag=zeros(meas_nr,part_nr);
C_pd=zeros(meas_nr,part_nr);
C_au=zeros(meas_nr,part_nr);
C_ag=zeros(meas_nr,part_nr);
ex_pd=zeros(meas_nr,part_nr);
ex_au=zeros(meas_nr,part_nr);
ex_ag=zeros(meas_nr,part_nr);
FWHM_pd=zeros(meas_nr,part_nr);
FWHM_au=zeros(meas_nr,part_nr);
FWHM_ag=zeros(meas_nr,part_nr);

 for i=1:part_nr
    [peak,pos_max]=max(norm_particles_pd(:,i));
    lam_max=lambda_pd(pos_max);
 [peak_lam_pd(:,i),C_pd(:,i),ex_pd(:,i),FWHM_pd(:,i)]=PeakFitLorentzCentroid(norm_particles_pd(:,i),lambda_pd,lam_max,50,1,0);
 end

 for i=1:part_nr
    [peak,pos_max]=max(norm_particles_au(:,i));
    lam_max=lambda_au(pos_max);
 [peak_lam_au(:,i),C_au(:,i),ex_au(:,i),FWHM_au(:,i),conf_au(:,i)]=PeakFitLorentzCentroid(norm_particles_au(:,i),lambda_au,lam_max,50,1,0);
 end

 for i=1:part_nr
    [peak,pos_max]=max(norm_particles_ag(:,i));
    lam_max=lambda_ag(pos_max);
 [peak_lam_ag(:,i),C_ag(:,i),ex_ag(:,i),FWHM_ag(:,i)]=PeakFitLorentzCentroid(norm_particles_ag(:,i),lambda_ag,lam_max,80,1,0);
 end


%% Round for table
clc
peak_lam_au=round(peak_lam_au);
FWHM_au=round(FWHM_au);
ex_au=round(ex_au,2);
peak_lam_pd=round(peak_lam_pd);
FWHM_pd=round(FWHM_pd);
ex_pd=round(ex_pd,2);
peak_lam_ag=round(peak_lam_ag);
FWHM_ag=round(FWHM_ag);
ex_ag=round(ex_ag,2);
%% Disk to rod with square wave pulses
clear all
clc
clf
load lamp_spectrum.txt
load disk_to_rod_Pd_h_square_wave_meas.txt
load disk_to_rod_Au_h_square_wave_meas.txt
load disc_to_rod_pd_square_wave_10percentageH2.txt
load isak_au_10per_600s.txt
load disk_to_rod_Au_h_square_wave_meas.txt
% % Au particles square pulses
[lambda_pd_ex,particles,background,lamp]= getData(disk_to_rod_Au_h_square_wave_meas,lamp_spectrum,14);
norm_particles_pd_ex=(particles-background)./lamp;

%Pd particles disk to rod square pulses
% [lambda_pd_ex,particles,background,lamp]= getData(disc_to_rod_pd_square_wave_10percentageH2,lamp_spectrum,14);
% norm_particles_pd_ex=(particles-background)./lamp;

%Number of wavelengths
int_nr=length(norm_particles_pd_ex(:,1,1));
%Number of particles
part_nr=length(norm_particles_pd_ex(1,1,:));
%Number of measurements
meas_nr=length(norm_particles_pd_ex(1,:,1));
%% Plot one particle square pulses
y_ax=ones(int_nr,1);
clf

figure(1)
for i=1:10:meas_nr
     plot3(lambda_pd_ex,y_ax*i,norm_particles_pd_ex(:,i,3),'color',[0 0.4470 0.7410],'Linewidth',2);
     hold on
     plot3(lambda_pd_ex,y_ax*i,norm_particles_pd_ex(:,i,12),'color',[0.8500 0.3250 0.0980],'Linewidth',2);
     legend('Pd: Disk- 1st particle','Pd: Disk- 10th particle')
     hold on
     grid on
end

figure(2)
for i=1:7:length(norm_particles_pd_ex(1,:,1))
     plot3(lambda_pd_ex,y_ax*i,norm_particles_pd_ex(:,i,12),'color',[0 0.4470 0.7410],'Linewidth',2);
     hold on
     grid on
end
grid on
set(gca,'FontSize',28)
yl=ylabel('Measurement','Interpreter','latex');
zl=zlabel('$I_{\textrm{norm}}$','Interpreter','latex');
xl=xlabel('$\lambda$ [nm]','Interpreter','latex');
xl.FontSize=34;
yl.FontSize=34;
zl.FontSize=34;

%% Square pulses of hydrogen
clc
clf

col = [0, 0, 1;
       0.6350 0.0780 0.1840;
       0.3010 0.7450 0.9330;
       0, 0, 0;
       1, 0, 1;
       0, 1, 1;
       0.6350 0.0780 0.184;
       0 0.4470 0.7410;
       0.8500 0.3250 0.0980;
       0.4660 0.6740 0.1880;
       0.4940 0.1840 0.5560;
       1, 0, 0;
       0, 1, 0;
       0.9290 0.6940 0.1250;
   ] ;        %Colors

%Get peaks
[lam,intensity]=get_peaks(norm_particles_pd_ex);
%Particle number
part_nr=length(norm_particles_pd_ex(1,1,:));
intensities=reshape(intensity(1,:,:),[meas_nr,part_nr]);



str=string(missing);
figure(1)
for i=1:part_nr
plot(intensities(2:5:end,i),'Linewidth',1.5,'color',col(i,:))

str(i)="Particle"+string(i);
hold on
end
[legh,objh]=legend(str,'Location', 'eastoutside','Interpreter','latex','LineWidth', 2,'FontSize',22);
lineh = findobj(objh,'type','line');
set(lineh,'Linewidth',2)

grid on
set(gca,'FontSize',22)
yl=ylabel('$I_{\textrm{p,norm}}$','Interpreter','latex');
xl=xlabel('Measurements','Interpreter','latex');
xl.FontSize=32;
yl.FontSize=32;

%% different size
clc
clf
clear all
load lamp_spectrum.txt
load start_100_19particles_nm_ag_100_meas.txt
[lambda,particles,background,lamp]= getData(start_100_19particles_nm_ag_100_meas,lamp_spectrum,18);
norm_particles=(particles-background)./lamp;
y_ax=ones(length(norm_particles),1);

figure(1)
for i=1:length(norm_particles(1,:))
     plot3(lambda,y_ax*i, norm_particles(:,i),'color',[0 0.4470 0.7410],'Linewidth',2);
     hold on
     grid on
end

%% Disk to rod with/without hydrogen
clc
clf
load start_100_19particles_nm_ag_100_meas.txt
load lamp_spectrum.txt
load disk_to_rod_16_particles_Ar90_H_10.txt
load disk_to_rod_16_particles_Ar.txt
load disk_to_rod_Au_14_particles_Ar90.txt
load disk_to_rod_Au_14_particles_Ar90_H10.txt

%Au particles
%Disk to rod 16 particles. These are reversed for some reason, i.e.
%background is at last position
% [lambda_90_10,particles_90_10,background_90_10,lamp]= getData(disk_to_rod_Au_14_particles_Ar90_H10,lamp_spectrum,14);
% norm_particles_90_10=(particles_90_10-background_90_10)./lamp;
%  
% [lambda_100,particles_100,background_100,lamp]= getData(disk_to_rod_Au_14_particles_Ar90,lamp_spectrum,14);
% norm_particles_100=(particles_100-background_100)./lamp;


%Disk to rod 16 particles. These are reversed for some reason, i.e.
%background is at last position
[lambda_90_10,particles_90_10,background_90_10,lamp]= getDataRev(disk_to_rod_16_particles_Ar90_H_10,lamp_spectrum,16);
norm_particles_90_10=(particles_90_10-background_90_10)./lamp;
 
[lambda_100,particles_100,background_100,lamp]= getDataRev(disk_to_rod_16_particles_Ar,lamp_spectrum,16);
norm_particles_100=(particles_100-background_100)./lamp;

clc
clf
y_ax=ones(length(norm_particles_100),1);

figure(1)
clf

figure(1)
for i=2:14
     plot3(lambda_100,y_ax*(i-1), norm_particles_100(:,i),'color',[0 0.4470 0.7410],'Linewidth',2);
      hold on
      grid on
end
grid on
set(gca,'FontSize',28)
legend('100\% Ar','Location', 'northeast','Interpreter','latex','LineWidth', 2,'FontSize',28);
yl=ylabel('Particle','Interpreter','latex');
zl=zlabel('$I_{\textrm{norm}}$','Interpreter','latex');
xl=xlabel('$\lambda$ [nm]','Interpreter','latex');
xl.FontSize=34;
yl.FontSize=34;
zl.FontSize=34;
% zlim([0 0.24])

figure(2)
for i=2:14
     plot3(lambda_90_10,y_ax*(i-1), norm_particles_90_10(:,i),'color',[0.8500 0.3250 0.0980],'Linewidth',2);
     hold on
     grid on
end
grid on
set(gca,'FontSize',28)
legend('90\% Ar and 10\% H$_2$','Location', 'northeast','Interpreter','latex','LineWidth', 2,'FontSize',28);
yl=ylabel('Particle','Interpreter','latex');
zl=zlabel('$I_{\textrm{norm}}$','Interpreter','latex');
xl=xlabel('$\lambda$ [nm]','Interpreter','latex');
xl.FontSize=34;
yl.FontSize=34;
zl.FontSize=34;
% zlim([0 0.24])

%% Isoterm: From first session
clear all
clc
clf
load lamp_spectrum.txt
load isoterm_180nm_pd_30deg.txt
load isoterm_first.txt

[lambda,particles,background,lamp]= getData(isoterm_180nm_pd_30deg,lamp_spectrum,19);
norm_particles=(particles-background)./lamp;


part_nr=length(norm_particles(1,1,:));
meas_nr=length(norm_particles(1,:,1));


figure(1)
data_top=zeros(1035,1);
%  for i=1:length(norm_particles(1,:))
k=1;
 for i=1:length(norm_particles(1,:,1))
    data_top(k)=max(norm_particles(:,i,1));
    k=k+1;
 end
 plot(data_top)
 
y_ax=ones(length(norm_particles(:,1,1)),1);
figure(2)
for i=1:40:length(norm_particles(1,:,1))
     plot3(lambda,y_ax*i,norm_particles(:,i,3),'color',[0 0.4470 0.7410],'Linewidth',2);
     hold on
     grid on
end

%% Isoterm: Partial pressure
clc

t=floor((0:10.207:10554.038))+1;
lt_ar=4*ones(1,length(t));
lt_h2=5*ones(1,length(t));

ix_ar = sub2ind(size(isoterm_first),t,lt_ar); 
ix_h2 = sub2ind(size(isoterm_first),t,lt_h2); 
 
% t=floor((0:10.207:10564.245));
 %Total flow T=A_max*p_ar+H2_max*p_h2 where A_max=200 ml/min, H2_max=70
 %ml/min and p_ar,p_h2 given in percent
 tot_flow=2*isoterm_first(ix_ar)+0.1*isoterm_first(ix_h2);
 
 %Get partial pressure assuming 1 bar in room, convert to mbar
 ppress=1e3*(0.1*isoterm_first(ix_h2)./tot_flow);
 g=get_peaks(norm_particles);
%% Fit peaks with Lorentzian
clc
 for i=1:19
    [peak,pos_max]=max(norm_particles(:,i,1));
    lam_max=lambda(pos_max);
 PeakFitGaussCentroid(norm_particles(:,i,1),lambda,lam_max,50,1,1)
 end
 

%% Plotting
clc
clf
[peak_0,intensity]=get_peaks(norm_particles);
 pks_median=peak_0(1,:,:);


particle_nr=10;
spec_to_plot=(pks_median(1,:,particle_nr)-pks_median(1,30,particle_nr))/pks_median(1,30,particle_nr);
plot(-spec_to_plot(30:970),ppress(30:970),'x')

ppress_r=round(ppress,1);
ppress_non_dub=unique(ppress_r);
%%
clc
plot(is_term,pressures,'x')
xlim([0 0.35])
%% Isoterm extra session + Partial pressure
clear all
clc
clf

load lamp_spectrum.txt
load isoterm_23T_Pd_120nm_isak_viktor_extrasession.txt
load Isoterm_isak_viktor_T30_extrasession.txt

[lambda,particles,background,lamp]= getData(isoterm_23T_Pd_120nm_isak_viktor_extrasession,lamp_spectrum,19);
norm_particles=(particles-background)./lamp;


t=floor((0:10.207:7747.113))+1;
lt_ar=4*ones(1,length(t));
lt_h2=5*ones(1,length(t));

ix_ar = sub2ind(size(Isoterm_isak_viktor_T30_extrasession),t,lt_ar); 
ix_h2 = sub2ind(size(Isoterm_isak_viktor_T30_extrasession),t,lt_h2); 
 
% t=floor((0:10.207:10564.245));
 %Total flow T=A_max*p_ar+H2_max*p_h2 where A_max=200 ml/min, H2_max=70
 %ml/min and p_ar,p_h2 given in percent
tot_flow=2*Isoterm_isak_viktor_T30_extrasession(ix_ar)+0.7*Isoterm_isak_viktor_T30_extrasession(ix_h2);
 
%Get partial pressure assuming 1 bar in room, convert to mbar
ppress=1e3*(0.7*Isoterm_isak_viktor_T30_extrasession(ix_h2)./tot_flow);
g=get_peaks(norm_particles);

part_nr=length(norm_particles(1,1,:));
meas_nr=length(norm_particles(1,:,1));
y_ax=ones(length(norm_particles(:,1,1)),1);

peak_int=zeros(meas_nr,1);
%%
clc
clf
figure(1)
for i=120:30:meas_nr
     plot3(lambda,y_ax*i,norm_particles(:,i,3),'Linewidth',2);
     plot3(y_ax*660.9020,lambda,y_ax*max(norm_particles(:,i,3)),'Linewidth',4,'color','k')
     hold on
     legend('Pd: Disk- 1st particle')
     hold on
     grid on
     pause(1)
     xlim([600 720])
end
xlim([600 720])
%%
[peak,pos_max]=max(norm_particles(:,420,1));
    lam_max=lambda(pos_max);
PeakFitLorentzCentroid(norm_particles(:,420,1),lambda,lam_max,95,1,1)


%% Fit data with Lorentzian
clc
norm_particles_to_plot=norm_particles(:,120:724,:);

peaks=zeros(length(norm_particles_to_plot(1,:,1)),part_nr);
C=zeros(length(norm_particles_to_plot(1,:,1)),part_nr);
ex=zeros(length(norm_particles_to_plot(1,:,1)),part_nr);
FWHM=zeros(length(norm_particles_to_plot(1,:,1)),part_nr);

for i=1:part_nr
    [peak,pos_max]=max(norm_particles_to_plot(:,1,i));
    lam_max=lambda(pos_max);
[peaks(:,i),C(:,i),ex(:,i),FWHM(:,i)]= PeakFitLorentzCentroid(norm_particles_to_plot(:,:,i),lambda,lam_max,120,1,1);
end
%%
clc
%% Save data
writematrix(peaks,'isoterm_plot_extra_delta_lam120pt2')
writematrix(C,'isoterm_plot_extra_delta_centroid120')
writematrix(ex,'isoterm_plot_extra_delta_ex120')
writematrix(FWHM,'isoterm_plot_extra_delta_FWHM120')

%% Data fitted with Lorentzian to plot
load isoterm_plot_extra.txt
load isoterm_plot_extra_delta_lam120.txt
load isoterm_plot_extra_delta_lam120pt2.txt
load isoterm_plot_extra_delta_centroid120.txt
load isoterm_plot_extra_delta_ex120.txt
load isoterm_plot_extra_delta_FWHM120.txt

peaks=isoterm_plot_extra_delta_lam120;
peaks_pt2=isoterm_plot_extra_delta_ex120;
FWHM=isoterm_plot_extra_delta_FWHM120;
clc
clf
particle=3;

figure(1)
% plot((peaks(1:16:end,particle)-peaks(1,particle)),ppress(120:16:724),'-x','Markersize',8,'Linewidth',2,'color',[0.8500 0.3250 0.0980])
plot((peaks(1:16:end,particle)-peaks(1,particle)),ppress(120:16:724),'-x','Markersize',8,'Linewidth',2)

ylim([0,72])
xlim([0 2.3])
grid on
set(gca,'FontSize',28)
yl=ylabel('$H_2 \textrm{ Partial pressure [mbar]}$','Interpreter','latex');
xl=xlabel('$\Delta\lambda_{max}$ [nm] ','Interpreter','latex');
xl.FontSize=34;
yl.FontSize=34;


figure(2)
plot(-(peaks_pt2(1:16:end,particle)-peaks_pt2(1,particle))/peaks_pt2(1,particle),ppress(120:16:724),'-x','Markersize',8,'Linewidth',2,'color',[0.8500 0.3250 0.0980])
% plot(-(peaks_pt2(1:16:end,particle)-peaks_pt2(1,particle))/peaks_pt2(1,particle),ppress(120:16:724),'-x','Markersize',8,'Linewidth',2)

ylim([0,72])
grid on
ax=set(gca,'FontSize',28);
yl=ylabel('$H_2 \textrm{ Partial pressure [mbar]}$','Interpreter','latex');
xl=xlabel('$\Delta \textrm{PI}$ (a.u.) ','Interpreter','latex');
xl.FontSize=32;
yl.FontSize=32;
 xticklabels({'0','-0.05','-0.1','-0.15','-0.2','-0.25','-0.3'})
%  xticklabels({'0','-0.1','-0.2','-0.3'})
xlim([0 0.32])

figure(3)
plot(FWHM(1:16:end,particle),ppress(120:16:724),'-x','Markersize',8,'Linewidth',2)
% plot(FWHM(1:16:end,particle),ppress(120:16:724),'-x','Markersize',8,'Linewidth',2,'color',[0.8500 0.3250 0.0980])
ylim([0,72])
grid on
set(gca,'FontSize',28)
yl=ylabel('$H_2 \textrm{ Partial pressure [mbar]}$','Interpreter','latex');
xl=xlabel('$ \textrm{FWHM}$ [nm] ','Interpreter','latex');
xl.FontSize=32;
yl.FontSize=32;

%%
clc

peak_0=zeros(part_nr,meas_nr);
for j=1:part_nr
        peak_0(j,:)=PeakFitLorentzCentroid(norm_particles(:,:,j),lambda,685,120,1,0);
end
%%
clc
clf
load Isoterm_delta_lambda.txt
peak0=Isoterm_delta_lambda;
% [peak_0,intensity]=get_peaks(norm_particles);
% pks_median=peak_0(1,:,:);
spec_to_plot=peak0;

particle_nr=1;
% spec_to_plot=(pks_median(1,:,particle_nr)-pks_median(1,140,particle_nr))/pks_median(1,140,particle_nr);

plot(spec_to_plot(140:724),ppress(140:724),'-x','Linewidth',1)
ylim([0 72])
grid on
set(gca,'FontSize',24)
yl=ylabel('$H_2 \textrm{ Partial pressure [mbar]}$','Interpreter','latex');
xl=xlabel('-$\Delta\lambda_{\textrm{norm}}$ ','Interpreter','latex');
xl.FontSize=28;
yl.FontSize=28;

% hold on
% plot(peak0(end-12,:),ppress,'x')

%%

lambda=linspace(300,700);
MLWA=getQ_ext(4,1,lambda);
plot(lambda,MLWA)
