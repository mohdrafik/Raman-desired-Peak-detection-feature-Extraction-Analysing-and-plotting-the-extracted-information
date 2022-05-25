
clc;
close all;
clear;

read_data1 = dlmread("20220523_009_Mod2_2_PC2_3mm_hole_focus_graphene_532nm_ND0.3_grating600_20x_arbpol_mono582_pin100_60s_+0.35V.txt");

% read_data1 = readmatrix("C:\Users\mbari\Downloads\Compressed\raman_data\20220523_009_Mod2_2_PC2_3mm_hole_focus_graphene_532nm_ND0.3_grating600_20x_arbpol_mono582_pin100_60s_+0.35V");
figure
plot(read_data1(:,1),read_data1(:,2));
[row,col] = size(read_data1);
% x_2D= read_data1(2640:2740,);
for i=1:row
    if (read_data1(i,1)>=2640)
            index_2Dmin = i ;            
            break;
    end
end
for i=1:row
    if (read_data1(i,1)>=2740)
            index_2Dmax = i ;
            break;
    end
end

fprintf("index_2Dmin = %d    %d \n", index_2Dmin,read_data1(index_2Dmin));
fprintf("index_2Dmax = %d    %d  \n", index_2Dmax ,read_data1(index_2Dmax));
Data_2DPeak = read_data1(index_2Dmin:index_2Dmax,1:2); % Data for 2D peak 

for i=1:row
    if (read_data1(i,1)>=1550)
            index_Gpeakmin = i ;
            break;
    end
end
for i=1:row
    if (read_data1(i,1)>=1650)
            index_GpeakMax = i ;
            break;
    end
end

fprintf("index_GpeakMin = %d   \t  %d  \n", index_Gpeakmin, read_data1(index_GpeakMax));
fprintf("index_GpeakMax = %d   \t  %d  \n", index_GpeakMax, read_data1(index_Gpeakmin));

Data_Gpeak = read_data1(index_Gpeakmin:index_GpeakMax,1:2);     % data for the Gpeak

figure 
subplot(2,1,1)
plot(read_data1(:,1),read_data1(:,2));
axis tight;
subplot(2,2,3)
plot(Data_2DPeak(:,1),Data_2DPeak(:,2));
 axis tight;
% ylim([min(read_data1(:,2)), max(read_data1(:,2))]);
subplot(2,2,4)
plot(Data_Gpeak(:,1),Data_Gpeak(:,2));
axis tight;

[row2d,col2d] =size(Data_2DPeak);
xydata = zeros(size(Data_2DPeak));
for i=1:row2d
   for j=1:col2d
 xydata(i,j) =  Data_2DPeak(i,j);
   end
end

 y1=smooth(Data_2DPeak(:,2)); 
 scatter(Data_2DPeak(:,1),Data_2DPeak(:,2)); 
 hold on ;
 plot(Data_2DPeak(:,1),y1);

 %% this is for lorentz fit use this inbuilt function // Moh Rafik
 
 [yfit, PARAMS, RESNORM ,RESIDUAL, JACOBIAN] =  lorentzfit(Data_2DPeak(:,1), Data_2DPeak(:,2));


%% these parameter need to be written in excel
x_peakL =PARAMS(2);  % x value of the peak ntb wrtn  which is PARAM(2);        % PARAMS =[p1 p2 p3 c];
c = PARAMS(4);   % constant value need to be wriitten.
Ypeak =  (max(yfit)-c); % I max or peak value need to be written.
FWHM =     2*sqrt(PARAMS(3));    % FWHM   is need to be written.

%% Area need to be written  and need to be calculated
% FWHM = 2*sqrt(PARAMS(3));


Y_FWHM = ((max(yfit)-c)./2)*ones(length(Data_2DPeak),1) + c ;
x_peakL_vect = x_peakL*ones(length(yfit(:,1)),1);


y_peakL = max(yfit) - PARAMS(4);

for i =1:length(yfit(:,1))
    if yfit(i,:) == PARAMS(2)
     xpeakLMan = i; 
     break;
    end 
end

% fprintf(" MANUALLY xpeak value of x and y  = %f  %d\n", xpeakLMan, (yfit(xpeakLMan)-PARAMS(3)) );
fprintf(" Automatically  xpeak value of x and y = %4d  %4d\n", x_peakL, y_peakL);


figure
plot(xydata(:,1), yfit); 
title(" yfit plot seperate")
hold on
plot(Data_2DPeak(:,1),y1,'ro');
% axis tight;
hold on 
plot(Data_2DPeak(:,1),Y_FWHM,'-mo');
hold on
plot(x_peakL_vect, yfit,'-g*');
axis tight;