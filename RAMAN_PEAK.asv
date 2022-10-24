
clc;
close all;
clear;

fileid = fopen("RamanData_file.txt",'w+');
% % fprintf(fileid, "voltage, HWFM, x_PEAK,y_PEAK,CONSTANT\n");
fprintf(fileid, "DATA FOR THE RAMAN 2D_peak \t\t\t\t\t\t\t\t\t\t\t");
fprintf(fileid, "DATA FOR THE RAMAN G_peak \n");
fprintf(fileid, "FWHM    x_peakL     Ypeak       c\t\t\t\t\t\t\t\t\t");
fprintf(fileid, "FWHM    x_peakL     Ypeak       c\n");
fclose(fileid);


%% here start your for loops

%try thsis way 
% % 
%  C={'a.dat', 'b.dat', 'c.dat',...};
%  for i = 1:n
%      fid(i) = fopen(C{i},'rt')
%  end
%

%  header = {'FWHM','x_peakL','Ypeak','constant'};
% fileid1 = xlswrite('RamanData_file1.xlsx',header, A1:D1);


read_data1 = dlmread("20220523_009_Mod2_2_PC2_3mm_hole_focus_graphene_532nm_ND0.3_grating600_20x_arbpol_mono582_pin100_60s_+0.35V.txt");
% read_data1 = readmatrix("C:\Users\mbari\Downloads\Compressed\raman_data\20220523_009_Mod2_2_PC2_3mm_hole_focus_graphene_532nm_ND0.3_grating600_20x_arbpol_mono582_pin100_60s_+0.35V");

figure
plot(read_data1(:,1),read_data1(:,2));
plot_datafile = i;
axis tight;
title("Raman Raw Data Plot ", plot_datafile);
[row,col] = size(read_data1);
% x_2D= read_data1(2640:2740,);



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

fprintf("index_GpeakMin = %.2f   \t  %.2f  \n", index_Gpeakmin, read_data1(index_GpeakMax));
fprintf("index_GpeakMax = %.2f   \t  %.2f  \n", index_GpeakMax, read_data1(index_Gpeakmin));

Data_Gpeak = read_data1(index_Gpeakmin:index_GpeakMax,1:2);     % data for the Gpeak


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

figure
subplot(2,2,[1,2])
plot(read_data1(:,1),read_data1(:,2));
title("Raman Raw Data")
subplot(2,2,3)
plot(Data_Gpeak(:,1),Data_Gpeak(:,2),'-r.');
title("G peak")
subplot(2,2,4)
plot(Data_2DPeak(:,1),Data_2DPeak(:,2),'-m.');
title("2D peak")
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
for i=1:2
    if i ==1
        fprintf("working for 2d peak i =1\n"); fprintf("i = %.1f\n",i);
        [yfit, PARAMS, RESNORM ,RESIDUAL, JACOBIAN] =  lorentzfit(Data_2DPeak(:,1), Data_2DPeak(:,2));  % this is for 2D peak selection data process.
        [row2d,col2d] =size(Data_2DPeak);
        xydata = zeros(size(Data_2DPeak));
        for ii=1:row2d
            for j=1:col2d
                xydata(ii,j) =  Data_2DPeak(ii,j);
            end
        end
    end

    if i ==2
        fprintf("working for G peak i =2\n"); fprintf("i = %.1f\n",i);
        [yfit, PARAMS, RESNORM ,RESIDUAL, JACOBIAN] =  lorentzfit(Data_Gpeak(:,1), Data_Gpeak(:,2));  % this is for G peak selection data process.
        [rowG,colG] =size(Data_Gpeak);
        xydata = zeros(size(Data_Gpeak));
        for ii=1:rowG
            for jj=1:colG
                xydata(ii,jj) =  Data_Gpeak(ii,jj);
            end
        end

    end

    % these parameter need to be written in excel
    x_peakL =PARAMS(2);  % x value of the peak ntb wrtn  which is PARAM(2);        % PARAMS =[p1 p2 p3 c];
    c = PARAMS(4);                 % constant value need to be wriitten.
    Ypeak =  (max(yfit)-c);        % I max or peak value need to be written.
    FWHM =   2*sqrt(PARAMS(3));    % FWHM   is need to be written.

    % Area need to be written  and need to be calculated
    % FWHM = 2*sqrt(PARAMS(3));

    P_mat_param = [x_peakL,Ypeak,c,FWHM];
    fprintf("parameter Display = \t");
    disp(P_mat_param);

    %     Y_FWHM = ((max(yfit)-c)./2)*ones(length(Data_2DPeak),1) + c ;
    Y_FWHM = ((max(yfit)-c)./2)*ones(length(yfit(:,1)),1) + c ;
    x_peakL_vect = x_peakL*ones(length(yfit(:,1)),1);


    y_peakL = max(yfit) - PARAMS(4);

    for ip =1:length(yfit(:,1))
        if yfit(ip,1) == PARAMS(2)
            xpeakLMan = ip;
            break;
        end
    end

    %     fprintf(" MANUALLY xpeak value of x and y  = %.2f  %.2f\n", xpeakLMan, (yfit(xpeakLMan)-PARAMS(3)));
    fprintf(" Automatically  xpeak value of x and y = %.2f  %.f\n", x_peakL, y_peakL);

    %% data is written in the file

    fileid = fopen("RamanData_file.txt",'a+');
    fprintf(fileid, "%.2f    %.2f    %.2f    %.2f\t\t\t\t\t\t\t\t",FWHM , x_peakL, Ypeak , c);
    % fprintf("\n");

    %% plotting the Peak and FWHM OF EACH 2D AND G PEAK
    figure
    plot(xydata(:,1), xydata(:,2),'bo');
    hold on
    plot(xydata(:,1), yfit(:,1));
    % axis tight;
    hold on
    plot(xydata(:,1),Y_FWHM,'-m.');
    hold on
    plot(x_peakL_vect, yfit,'-g.');
    axis tight;
    switch(i)
        case 1
            title(" Lorentz Fit Data for 2D PEAK", plot_datafile);
        case 2
            title(" Lorentz Fit Data for G PEAK " , plot_datafile);

    end
end

fprintf("\n")




%% here end for loop

% after for loop completed
% fclose(fileid);