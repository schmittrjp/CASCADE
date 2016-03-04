function [hydraulicData] = hydraulicCalc(AggData,plot)
%HYDRAULICSOLVER Summary of this function goes here
%   Detailed explanation goes here

%% hydraulic calculations 
disp('Running hydraulic calculations')

% general variables 
global ID_arcid ID_FromN ID_ToN ID_ElUs ID_ElUsRaw ID_ElDs ID_ElDsRaw ID_Slp ID_SlpRaw ID_ElDiff ID_Length ID_StrO ID_MicroWSAre ID_FldPlnWdth ID_Ad ID_FX ID_FY ID_TX ID_TY ID_Wac ID_Q15% Clear temporary variables
global kk wac q15 slp v kst_analytic rho s taucrit fminout d90
%some hydraulic constants constants for the hydraulic solver 
global  wac q15 slp  kst rho s fminout kk taucrit v kst_analytic d90
rho=1000; % density, water 
s=2600; % relative density, sediment 
kst=35; % Strickler coefficient

%% define options for hydraulic solving procedure
options = optimset;
 options.Algorithm = 'interior-point';
%  options.Algorithm = 'sqp';

options.GradObj = 'off';
options.Hessian = 'lbfgs';
options.Display = 'off';
options.MaxIter = 50000;
options.MaxFunEvals = 50000;
options.TolFun = 0.00000000001;
options.TolX = 0.00000000001;
%options.UseParallel='always';
option.OutputFcn=@outfun;
Hyd_out_store=zeros(length(AggData),3);
%% Loop through records for all reaches
reach_with_high_v=0;

for ii=1:length(AggData)
waitbar(ii/length(AggData));
wac=AggData(ii,ID_Wac); % multiply with 2, because the active channel measurments are only 1-sided
q15=AggData(ii,ID_Q15);
slp=AggData(ii,ID_Slp);

ii;
%define initial conditions
% h_init=wac*0.1;
h_init=q15.^0.299*slp.^-0.206; % from: 1. Huang HQ, Nanson GC (2002) A stability criterion inherent in laws governing alluvial channel flow. Earth Surface Processes and Landforms 27(9):929–944.: p 939, Tbale II

%solve for h
kk=1;
taucrit=0.047;
% [h,fval,exitflag,output]=fmincon(@hydraulicSolver2,h_init,[],[],[],[],0.1,20,[],options);
[h,fval,exitflag,output]=fmincon(@hydraulicSolver,h_init,[],[],[],[],0.1,20,[],options);


if taucrit>1
asd=1;
end

if h>20
asd=0; 
end

  
%calculate v 
v=q15/(h*wac);
if v>20
asd=0; 
end


% % estimate d90 using the common critical shear stress (0.047)
% d90=(rho*h*slp)/((s-rho)*taucrit);
% 
% % check Suska and Brownlie conditions..... 
% 
% if h/d90<10  % for low water levels in comparison to the grain size, use the formula proposed by Suszka L (1991) Modification of transport rate formula for steep channels. Fluvial Hydraulics of Mountain Regions, Lecture Notes in Earth Sciences., eds Armanini PA, Silvio PGD (Springer Berlin Heidelberg), pp 59–70. Available at: http://link.springer.com/chapter/10.1007/BFb0011182 [Accessed March 2, 2015].
%     taucrit=0.0851*(h/d90)^(-0.0261);
% else  
%     Re_d50=(1.6*9.81*d90)^0.5*d90/10^(-6); % particle Reynolds number
%     taucrit=(0.22*Re_d50^(-0.6)+0.06*10^(-7.7*Re_d50^(-0.6)));
% end
% 
% % .... and estimate d90 using the common critical shear stress (0.047)
% d90=(rho*h*slp)/((s-rho)*taucrit);

%kst=21.2/d90.^(1/6)

% if v>10; h=q15^0.299*slp^-0.206; v=q15/(h*wac); reach_with_high_v=reach_with_high_v+1;  end  % for some reaches there is no credible solution (very high velocities)

kst_analytical=21.2/(d90^(1/6));

Fr=v/sqrt(9.81*h);
if h>20
asd=0; 
end
Hyd_out_store(ii,1)=h; % Output storage for hydrological variables 
%if d90 > 2E-2; Hyd_out_store(ii,2)=d90/2.1; else Hyd_out_store(ii,2)=d90;  end ; %for gravel/cobble: d50 is ca. half of the d90. 
Hyd_out_store(ii,2)=d90/2.1;
Hyd_out_store(ii,3)=kst_analytical;
Hyd_out_store(ii,4)=v;
Hyd_out_store(ii,5)=Fr;
Hyd_out_store(ii,6)=taucrit;
% Hyd_out_store(ii,6)=exitflag;
end 
hydraulicData=Hyd_out_store;

%% plotting of results 

if plot==1
figure('Name','Hydraulic characteristics') 
xlabels_hyd={'Flow stage [m]', 'D_{50} [m]', 'K_{Strickler}', 'Velocity [m s^{-1}]','Froude','\tau_{* crit}'}
    for ff=1:6
        subplot(1,6,ff)
        boxplot(Hyd_out_store(:,ff))
        xlabel(xlabels_hyd{ff})
        if ff==2; 
           % find the outliers in grain size 
            h = findobj(gcf,'tag','Outliers');
            xdata = get(h,'XData');
            ydata = get(h,'YData');
            outlier_d50=find(ismember(Hyd_out_store(:,2),ydata{1})); % find the reaches which grain size is classified as outliers 
           % find the whiskers
            h = findobj(gcf,'tag','Upper Whisker');
            ydata = get(h,'YData');
            set(gca,'Ylim',[0 ydata{1}(2)]);
        end
      
    end

figure('Name','Hydraulic correlation') 
%  [h,~, AXpm]=gplotmatrix( [Hyd_out_store(:,2) raw_data_sort(:,ID_meanWac) raw_data_sort(:,ID_Q15) raw_data_sort(:,ID_Slp)]);    
% xlabel(AXpm(:,1),'D90', 'W_{AC}', 'Q_{1.5}', 'Slp')
% xlabels_hyd{[1 2 4 5 6]}
data_plot_mat=Hyd_out_store;
data_plot_mat(outlier_d50,:)=[];
[h,~, AXpm]=gplotmatrix( data_plot_mat(:,[1 2 5 6]),[],[],[],[],[],[],[],xlabels_hyd([1 2 5 6]));    
clear data_plot_mat h
% Channel_hyd_ANN_inputs=[raw_data(:,ID_meanWac) raw_data(:,ID_Q15) raw_data(:,ID_Slp) Hyd_out_store];
% D50=Hyd_out_store(:,2);
% 
% figure
% scatter(Hyd_out_store(:,2)./Hyd_out_store(:,1),raw_data_sort(:,ID_Slp))
% figure
% scatter(raw_data_sort(:,ID_Slp),Hyd_out_store(:,2))
% figure
% scatter(Hyd_out_store(:,1),Hyd_out_store(:,6))

set(gca,'Yscale','log','Xscale','log')

grainsizeClasses=fliplr([100 56 25 10 5 2 0.5 0.1])/1000;
% networkPlotter(AggData,hydraulicData(:,2),grainsizeClasses, AggData(:,ID_StrO),'Grain size distribution','Grain size [m]',0)

%% store data
% save('HYMO_and_con_data')
end 
end

