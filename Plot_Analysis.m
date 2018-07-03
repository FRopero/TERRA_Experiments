%Script to show three plots results of the algorithm

%Input
dir = 'D:\OneDrive\Universidad\Doctorado\Matlab\Results\'; %Directory where the data results are saved
it = 250; %Number of executions of the algorithm

%Plot 2D Delta=2,4,8,16,32,64 ; R=2 ; N = 64
files = {'DFF_1';'DFF_2';'DFF_3';'DFF_4';'DFF_5';'DFF_6'};
Plot2D(dir,files,0);

%ScatterGOA Delta=1 ; R=2,4,8,16 ; N = 6
files = {'DRR_1';'DRR_2';'DRR_3';'DRR_4'};
ScatterGOA(dir,files,it,10);

%Distance Scatter R =1,3,9 ; Delta = 2,4,8
files = {'SPP_9';'SPP_8';'SPP_7';'SPP_6';'SPP_5';'SPP_4';'SPP_3';'SPP_2';'SPP_1'};
Scatter(dir,files,it,'d',1,10,false);

%Test Anova Delta = 8
files = {'SPP_9';'SPP_6';'SPP_3'};
Anova(dir,files,it);

%Test Anova Delta = 4
files = {'SPP_8';'SPP_5';'SPP_2'};
Anova(dir,files,it);

%Test Anova Delta = 2
files = {'SPP_7';'SPP_4';'SPP_1'};
Anova(dir,files,it);

function [] = Scatter(dir,files,it,key,alg,scale,violin)

%INPUTS
% dir = directory where the result files are stored
% files = result files
% it = nº of instances for each file
% key = parameter to plot (d=distance;t=time;s=stop)
% alg = algoritm type to show its results (0=all;1=No_GOA;2=GOA(GC);3=GOA(HC);4=GOA(MC))
% scale = scale for the ratius=fugv/fuav
% violin = to show as a violin diagram (true or false)

    %Distance Scatter Plot 
    f = figure;
    p = uipanel('Parent',f,'BorderType','none'); 
    if (key=='d')
        sct = 'Distance';
    elseif (key=='t')
        sct = 'Time';
    elseif (key=='s')
        sct =  'Stops';
    end
    if (alg==0)
        g = 'by algorithm';
    elseif (alg==1)
        g = 'without GOA';
    elseif (alg==2)
        g = 'with GOA(GC)';
    elseif (alg==3)
        g = 'with GOA(HC)';
    elseif (alg==4)
        g = 'with GOA(MC)';
    end
        
    p.Title = [sct ' Scatter Plot ' g]; 
    p.TitlePosition = 'centertop'; 
    p.FontSize = 12;
    p.FontWeight = 'bold';
    
    distance_ugv = [];
    distance_uav = [];
    total_distance = [];
    single_total_distance = [];
    time_ugv = [];
    time_uav = [];
    total_time = [];
    single_total_time = [];
    stops = [];
    single_total_stops = [];
    max = 0;
    min = +Inf;
    
    for k=1:length(files)
        %Reset  variables
        distance_ugv = [];
        distance_uav = [];
        total_distance = [];
        single_total_distance = [];
        time_ugv = [];
        time_uav = [];
        total_time = [];
        single_total_time = [];
        stops = [];
        single_total_stops = [];
        
        %Load file
        f = [dir strjoin(files(k,1)) '\data.mat'];
        load(f);
        
        for i=1:it
            if (key=='d')
                distance_ugv = [distance_ugv ; data_results(i).It.f_ugv_d];
                distance_uav = [distance_uav ; data_results(i).It.f_uav_d];
                total_distance = [total_distance ; data_results(i).It.ftotal_d];
                if (alg>0)
                    single_total_distance(i,1) = distance_ugv(i,alg);
                    single_total_distance(i,2) = distance_uav(i,alg);
                    single_total_distance(i,3) = (distance_ugv(i,alg)/distance_uav(i,alg))*scale;
                    single_total_distance(i,4) = total_distance(i,alg);
                end
                if (alg>0)

                    if (single_total_distance(i,4) > max)
                        max = single_total_distance(i,4);
                    end
                    if (single_total_distance(i,2) < single_total_distance(i,1))
                        if (single_total_distance(i,2) < min)
                            min = single_total_distance(i,2);
                        end
                    else
                        if (single_total_distance(i,1) < min)
                            min = single_total_distance(i,1);
                        end
                    end
                else
                    td = [];
                    td = [td data_results(i).It.ftotal_d];
                    for j=1:4
                        if (td(j) > max)
                            max = td(j);
                        end
                        if (td(j) < min)
                            min = td(j);
                        end
                    end
                end
            elseif (key=='t')
                time_ugv = [time_ugv ; data_results(i).It.f_ugv_t];
                time_uav = [time_uav ; data_results(i).It.f_uav_t];
                total_time = [total_time ; data_results(i).It.ftotal_t];
                if (alg>0)
                    single_total_time(i,1) = time_ugv(i,alg);
                    single_total_time(i,2) = time_uav(i,alg);
                    single_total_time(i,3) = (time_ugv(i,alg)/time_uav(i,alg))*scale;
                    single_total_time(i,4) = total_time(i,alg);
                end
                if (alg>0)
                    if (single_total_time(i,4) > max)
                        max = single_total_time(i,4);
                    end
                    if (single_total_time(i,2) < single_total_time(i,1))
                        if (single_total_time(i,2) < min)
                            min = single_total_time(i,2);
                        end
                    else
                        if (single_total_time(i,1) < min)
                            min = single_total_time(i,1);
                        end
                    end
                else
                    td = [];
                    td = [td data_results(i).It.ftotal_t];
                    for j=1:4
                        if (td(j) > max)
                            max = td(j);
                        end
                        if (td(j) < min)
                            min = td(j);
                        end
                    end
                end
            elseif (key=='s')
                stops = [stops ; data_results(i).It.stops];
                if (alg>0)
                    single_total_stops(i,1) = stops(i,alg);
                    if (single_total_stops(i,1)>max)
                        max = single_total_stops(i,1);
                    end
                    if (single_total_stops(i,1)<min)
                        min = single_total_stops(i,1);
                    end
                else
                    td = [];
                    td = [td data_results(i).It.stops];
                    for j=1:4
                        if (td(j) > max)
                            max = td(j);
                        end
                        if (td(j) < min)
                            min = td(j);
                        end
                    end
                end
            end     
        end
        
        %pos = mod(2+k,9)+1;
        pos = k;
        if (pos>0)
            ah(k) = subplot(3,3,pos,'Parent',p);
            if (alg>0)
                if (key=='d')
                    boxplot(single_total_distance,'labels',{'f_ugv','f_uav','Ratio','f_total'})
                    title(sprintf('Boxplot %s',string(files(k,1))));
                    if (violin)
                        distributionPlot(single_total_distance); 
                    end
                elseif (key=='t')
                    boxplot(single_total_time,'labels',{'f_ugv','f_uav','Ratio','f_total'})
                    if (violin)
                        distributionPlot(single_total_time); 
                    end
                elseif (key=='s')
                    boxplot(single_total_stops,'labels',{'stops'})
                    if (violin)
                        distributionPlot(single_total_stops); 
                    end
                end
            else
                if (key=='d')
                    boxplot(total_distance,'labels',{'GA','GOA(GC)','GOA(HC)','GOA(MC)'});
                    if (violin)
                        distributionPlot(total_distance); 
                    end
                elseif (key=='t')
                    boxplot(total_time,'labels',{'GA','GOA(GC)','GOA(HC)','GOA(MC)'});
                    if (violin)
                        distributionPlot(total_time); 
                    end
                elseif (key=='s')
                    boxplot(stops,'labels',{'GA','GOA(GC)','GOA(HC)','GOA(MC)'});
                    if (violin)
                        distributionPlot(stops); 
                    end
                end
            end
        end
    end
    
    %Set interval
    set(ah,'ylim',[min max]) 
    
end

function [p_valor] = Anova(dir,files,it)

    y_anova = [];
    
    for k=1:length(files)
        total_distance = [];
        %Load file
        f = [dir strjoin(files(k,1)) '\data.mat'];
        load(f);
        
        for i=1:it
           total_distance = [total_distance ; data_results(i).It(1).ftotal_d];
        end
        y_anova = [y_anova total_distance];
    end
    
    p_valor = anova1(y_anova);
    fprintf('p_valor = %1.6f\n',p_valor);
end

function [] = Plot2D(dir,files,alg)
    
    f = figure;
    p = uipanel('Parent',f,'BorderType','none'); 
    if (alg==0)
        g = 'by algorithm';
    elseif (alg==1)
        g = 'without GOA';
    elseif (alg==2)
        g = 'with GOA(GC)';
    elseif (alg==3)
        g = 'with GOA(HC)';
    elseif (alg==4)
        g = 'with GOA(MC)';
    end
        
    %p.Title = [' Fuav Plot ' g]; 
    p.TitlePosition = 'centertop'; 
    p.FontSize = 12;
    p.FontWeight = 'bold';

    y_uav = [];
    y_ugv = [];
    y_total = [];
    stops = [];
    y_uav1 = [];
    y_ugv1 = [];
    y_total1 = [];
    stops1 = [];
    y_uav2 = [];
    y_ugv2 = [];
    y_total2 = [];
    stops2 = [];
    y_uav3 = [];
    y_ugv3 = [];
    y_total3 = [];
    stops3 = [];
    y_uav4 = [];
    y_ugv4 = [];
    y_total4 = [];
    stops4 = [];
    
    x = [2,4,8,16,32,64];
    for k=1:length(files)
        %Load file
        f = [dir strjoin(files(k,1)) '\data.mat'];
        load(f);
        distance_ugv = [];
        distance_uav = [];
        distance_total = [];
        st = [];
        for i=1:it_completed
            distance_uav = [distance_uav ; data_results(i).It.f_uav_d];
            distance_ugv = [distance_ugv ; data_results(i).It.f_ugv_d];
            distance_total = [distance_total; data_results(i).It.ftotal_d];
            st = [st ; data_results(i).It.stops];
        end
        y_uav1 = [y_uav1 median(distance_uav(:,1))];
        y_ugv1 = [y_ugv1 median(distance_ugv(:,1))];
        y_total1 = [y_total1 median(distance_total(:,1))];
        stops1 = [stops1 median(st(:,1))];
        y_uav2 = [y_uav2 median(distance_uav(:,2))];
        y_ugv2 = [y_ugv2 median(distance_ugv(:,2))];
        y_total2 = [y_total2 median(distance_total(:,2))];
        stops2 = [stops2 median(st(:,2))];
        y_uav3 = [y_uav3 median(distance_uav(:,3))];
        y_ugv3 = [y_ugv3 median(distance_ugv(:,3))];
        y_total3 = [y_total3 median(distance_total(:,3))];
        stops3 = [stops3 median(st(:,3))];
        y_uav4 = [y_uav4 median(distance_uav(:,4))];
        y_ugv4 = [y_ugv4 median(distance_ugv(:,4))];
        y_total4 = [y_total4 median(distance_total(:,4))];
        stops4 = [stops4 median(st(:,4))];
    end
  
    if (alg==1)
        g = 'without GOA';
    elseif (alg==2)
        g = 'with GOA(GC)';
    elseif (alg==3)
        g = 'with GOA(HC)';
    elseif (alg==4)
        g = 'with GOA(MC)';
    end
    
    y_ugv = [y_ugv1;y_ugv2;y_ugv3;y_ugv4];
    subplot(1,4,1,'parent',p)
    plot(x,y_ugv,'-+');
    xlabel('Delta');
    ylabel('Fugv');
    
    y_uav = [y_uav1;y_uav2;y_uav3;y_uav4];
    subplot(1,4,2,'parent',p)
    plot(x,y_uav,'-+');
    xlabel('Delta');
    ylabel('Fuav');
    
    y_total = [y_total1;y_total2;y_total3;y_total4];
    subplot(1,4,3,'parent',p)
    plot(x,y_total,'-+');
    xlabel('Delta');
    ylabel('Ftotal');
    
    stops = [stops1;stops2;stops3;stops4];
    subplot(1,4,4,'parent',p)
    plot(x,stops,'-*');
    xlabel('Delta');
    ylabel('Stops');
    axis([0 70 15 25]);
    
end

function [] = Scatter3D(dir,files,it,alg)
% Plot 3d Ratio = Fugv/Fuav to see the best solution based on the Fugv and
% Fuav.
% Eje X = Fugv
% Eje Y = Fuav
% Eje Z = Ftotal
x = [];
y = [];
z = [];

for k=1:length(files)
    %Load file
    f = [dir strjoin(files(k,1)) '\data.mat'];
    load(f);
    distance_uav = [];
    distance_ugv = [];
    total_distance = [];
    for i=1:it
        distance_uav = [distance_uav ; data_results(i).It.f_uav_d];
        distance_ugv = [distance_ugv ; data_results(i).It.f_ugv_d];
        total_distance = [total_distance ; data_results(i).It.ftotal_d];
    end
    x = [x distance_ugv(:,alg)'];
    y = [y distance_uav(:,alg)'];
    %z = [z total_distance(:,alg)'];
    tmp = [];
    if (k==1 || k==10)
        for t=1:it
            tmp = [tmp 2];
        end
    elseif (k==2 || k==12)
        for t=1:it
            tmp = [tmp 4];
        end
    elseif (k==3 || k==14)
        for t=1:it
            tmp = [tmp 8];
        end
    elseif (k==4 || k==11)
        for t=1:it
            tmp = [tmp 2/3];
        end
    elseif (k==5 || k==13)
        for t=1:it
            tmp = [tmp 4/3];
        end
    elseif (k==6 || k==15)
        for t=1:it
            tmp = [tmp 8/3];
        end            
    elseif (k==7)
        for t=1:it
           tmp = [tmp 2/4];
        end
    elseif (k==8)
        for t=1:it
           tmp = [tmp 1];
        end  
    elseif (k==9)
        for t=1:it
           tmp = [tmp 2];
        end
    end
    z = [z tmp];
end

figure;
scatter3(x,y,z,'.');
xlabel('Fugv');
ylabel('Fuav');
zlabel('Ftotal');
title('TERRA GOA(MC)');
end

function [] = ScatterGOA(dir,files,it,scale)
    %f = figure;
    %p = uipanel('Parent',f,'BorderType','none'); 
    %p.Title = 'TERRA Scatter Plot'; 
    %p.TitlePosition = 'centertop'; 
    %p.FontSize = 12;
    %p.FontWeight = 'bold';
    %xlabel('Radius');
    %ylabel('TERRA');
    minn = [];
    maxx = [];
    fugv = [];
    fugv1 = [];
    fuav = [];
    fuav1 = [];
    ftotal = [];
    ftotal1 = [];
    t = 1;
    
    for k=1:length(files)
        %Reset  variables
        distance_ugv = [];
        distance_uav = [];
        total_distance = [];
        total_distance_a1 = [];
        total_distance_a2 = [];
        total_distance_a3 = [];
        total_distance_a4 = [];
        
        %Load file
        f = [dir strjoin(files(k,1)) '\data.mat'];
        load(f);
        
        %Load Data File
        for i=1:it
            ma = [];
            mi = [];
            
            distance_ugv = [distance_ugv ; data_results(i).It.f_ugv_d];
            distance_uav = [distance_uav ; data_results(i).It.f_uav_d];
            total_distance = [total_distance ; data_results(i).It.ftotal_d];

            total_distance_a1(i,1) = distance_ugv(i,1);
            total_distance_a1(i,2) = distance_uav(i,1);
            total_distance_a1(i,3) = (distance_ugv(i,1)/distance_uav(i,1))*scale;
            total_distance_a1(i,4) = total_distance(i,1);
            
            mi = [mi min(total_distance_a1(i,:))];
            ma = [ma max(total_distance_a1(i,:))];
            
            total_distance_a2(i,1) = distance_ugv(i,2);
            total_distance_a2(i,2) = distance_uav(i,2);
            total_distance_a2(i,3) = (distance_ugv(i,2)/distance_uav(i,2))*scale;
            total_distance_a2(i,4) = total_distance(i,2);
            
            mi = [mi min(total_distance_a2(i,:))];
            ma = [ma max(total_distance_a2(i,:))];
            
            total_distance_a3(i,1) = distance_ugv(i,3);
            total_distance_a3(i,2) = distance_uav(i,3);
            total_distance_a3(i,3) = (distance_ugv(i,3)/distance_uav(i,3))*scale;
            total_distance_a3(i,4) = total_distance(i,3);
            
            mi = [mi min(total_distance_a3(i,:))];
            ma = [ma max(total_distance_a3(i,:))];
            
            total_distance_a4(i,1) = distance_ugv(i,4);
            total_distance_a4(i,2) = distance_uav(i,4);
            total_distance_a4(i,3) = (distance_ugv(i,4)/distance_uav(i,4))*scale;
            total_distance_a4(i,4) = total_distance(i,4);
               
            mi = [mi min(total_distance_a4(i,:))];
            ma = [ma max(total_distance_a4(i,:))];

            minn = [minn min(mi)];
            maxx = [maxx max(ma)];
                                    
        end
        
        %Show Data File
        %ah(t) = subplot(4,4,k,'parent',p);
        %boxplot(total_distance_a4,'labels',{'f_ugv','f_uav','Ratio','f_total'});
        fugv(k,3) = median(total_distance_a4(:,1));
        fuav(k,3) = median(total_distance_a4(:,2));
        ftotal(k,3) = median(total_distance_a4(:,4));
        t = t + 1;
        %ah(t) = subplot(4,4,k+4,'parent',p);
        %boxplot(total_distance_a3,'labels',{'f_ugv','f_uav','Ratio','f_total'});
        t = t + 1;
        fugv(k,2) = median(total_distance_a3(:,1));
        fuav(k,2) = median(total_distance_a3(:,2));
        ftotal(k,2) = median(total_distance_a3(:,4));
        %ah(t) = subplot(4,4,k+8,'parent',p);
        %boxplot(total_distance_a2,'labels',{'f_ugv','f_uav','Ratio','f_total'});
        t = t + 1;
        fugv(k,1) = median(total_distance_a2(:,1));
        fuav(k,1) = median(total_distance_a2(:,2));
        ftotal(k,1) = median(total_distance_a2(:,4));
        %ah(t) = subplot(4,4,k+12,'parent',p);
        %boxplot(total_distance_a1,'labels',{'f_ugv','f_uav','Ratio','f_total'});
        t = t + 1;
        fugv1(1,k) = median(total_distance_a1(:,1));
        fuav1(1,k) = median(total_distance_a1(:,2));
        ftotal1(1,k) = median(total_distance_a1(:,4));
    end
    
    %Set interval
    %set(ah,'ylim',[min(minn) max(maxx)]) 
    
    f1 = figure;
    p1 = uipanel('Parent',f1,'BorderType','none'); 
    p1.Title = 'GOA Plot'; 
    p1.TitlePosition = 'centertop'; 
    p1.FontSize = 12;
    p1.FontWeight = 'bold';
    
    %%% Fugv %%%
%     for h=1:4
%         fugv(h,4) = mean(fugv(h,1:3));
%     end
    
    y =  fugv1-fugv';
    subplot(2,3,1,'parent',p1);
    plot([2,4,8,16],y);
    title('Y = Fugv(No GOA)-Fugv(GOA)')
    xlabel('R (km)')
    ylabel('Y (km)')
    
    %%% Fuav %%%
%     for h=1:4
%         fugv(h,4) = mean(fugv(h,1:3));
%     end
    
    y =  fuav1-fuav';
    subplot(2,3,2,'parent',p1);
    plot([2,4,8,16],y);
    title('Y = Fuav(No GOA)-Fuav(GOA)')
    xlabel('R (km)')
    ylabel('Y (km)')
    
    %%% Ftotal %%%
%     for h=1:4
%         fugv(h,4) = mean(fugv(h,1:3));
%     end
    
    y =  ftotal1-ftotal';
    subplot(2,3,3,'parent',p1);
    plot([2,4,8,16],y);
    title('Y = Ftotal(No GOA)-Ftotal(GOA)')
    xlabel('R (km)')
    ylabel('Y (km)')
    
    Vugv = 0.13;
    Vuav = 30;
    
    tugv = fugv ./Vugv;
    tugv1 = fugv1 ./Vugv;
    
    y = tugv1-tugv';
    
    subplot(2,3,4,'parent',p1);
    plot([2,4,8,16],y);
    title('Y = Tugv(No GOA)-Tugv(GOA)')
    xlabel('R (km)')
    ylabel('Y (h)')
    
    tuav = fuav ./Vuav;
    tuav1 = fuav1 ./Vuav;
    y = tuav1-tuav';
    subplot(2,3,5,'parent',p1);
    plot([2,4,8,16],y);
    title('Y = Tuav(No GOA)-Tuav(GOA)')
    xlabel('R (km)')
    ylabel('Y (h)')
    
    Ttotal = tugv+tuav;
    Ttotal1 = tugv1 + tuav1;
    y = Ttotal1-Ttotal';
    subplot(2,3,6,'parent',p1);
    plot([2,4,8,16],y);
    title('Y = Ttotal(No GOA)-Ttotal(GOA)')
    xlabel('R (km)')
    ylabel('Y (h)')
    
    
end