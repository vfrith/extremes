% This script carries out the analysis on the regional daily data from ERA5
% or ERA-Heat produced on BluePebble using reg_mean_calc_ERA5.m and
% reg_mean_calc_ERAHeat.m.


%% Load data if required
% Load the data produced on BluePebble and scped locally
if ~exist('reg_daily_05','var')
    load('reg_daily_05.mat')
    load('reg_daily_2.mat')
    load('reg_daily_WWA.mat')
    load('reg_daily.mat')
    load('mon_val.mat')
end

%% Create a useable date array
years = nan(length(mon_val),1);
year = 1979;
years(1) = year;

for i = 2:length(mon_val)
    if mon_val(i) == 1 && mon_val(i-1) == 12
        year = year + 1;
    end
    years(i) = year;
end


%% Extremes analysis -> tasmax
% Find months of maximum temperatures
[~,max_ids] = max(reg_daily_05,[],2);

% Take just the baseline period (1979-2020, 42 years)
reg_daily_baseline = reg_daily_05(:,1:15341);
mon_baseline = mon_val(1:15341);

% Create empty arrays for output
reg_05_sd = nan(237,1);
reg_05_mean = nan(237,1);
reg_05_exm = nan(237,1);
reg_05_exm_daily = nan(237,length(mon_val));

% Find long-term SD and mean for calendar month of max for each region
for r = 1:237
    
    % For each region, find which month had max temp
    max_mon = mon_val(max_ids(r));
    
    % Select only that month
    mon_ids = mon_baseline == max_mon;
    
    % Take SD and mean for that month
    reg_05_sd(r) = std(reg_daily_baseline(r,mon_ids));
    reg_05_mean(r) = mean(reg_daily_baseline(r,mon_ids));
    
    % Calculate extreme magnitude (following Vikki's definition)
    reg_05_exm(r) = (max(reg_daily_05(r,:)) - reg_05_mean(r)) ./ reg_05_sd(r);
    
    % Go through every day to calculate the metric
    reg_05_exm_daily(r,:) = (reg_daily_05(r,:) - reg_05_mean(r)) ./ reg_05_sd(r);
end
    
% % Save the files
% save('reg_sd.mat','reg_sd')
% save('reg_mean.mat','reg_mean')
% save('reg_exm.mat','reg_exm')
% save('reg_exm_daily.mat','reg_exm_daily')


%% As above but with a moving baseline (10 year trailing)
% Find months of maximum temperatures
[~,max_ids] = max(reg_daily_05,[],2);

% Create empty arrays for output
reg_05_sd_running = nan(237,33);
reg_05_mean_running = nan(237,33);
reg_05_exm_running = nan(237,33);
reg_05_exm_daily_running = nan(237,length(mon_val));


% Find long-term SD and mean for calendar month of max for each region
for r = 1:237
    
    % For each region, find which month had max temp
    max_mon = mon_val(max_ids(r));
        
    for y = 1988:2020 % Do a 10-year moving baseline
        
        y0 = y-9;
        tstart = min(find(years == y0));
        tend = max(find(years == y));
        tend2 = max(find(years == y+1));
        
        
        % Take just the baseline period (1979-2020, 42 years)
        reg_daily_baseline = reg_daily_05(:,tstart:tend);
        mon_baseline = mon_val(tstart:tend);
 
        % Select only maximum month
        mon_ids = mon_baseline == max_mon;
        
        % Take SD and mean for that month
        reg_05_sd_running(r,y-1987) = std(reg_daily_baseline(r,mon_ids));
        reg_05_mean_running(r,y-1987) = mean(reg_daily_baseline(r,mon_ids));
        
        % Calculate extreme magnitude (following Vikki's definition)
        reg_05_exm_running(r,y-1987) = (max(reg_daily_05(r,tend+1:tend2)) - reg_05_mean_running(r,y-1987)) ./ reg_05_sd_running(r,y-1987);
        
        % Go through every day to calculate the metric
        reg_05_exm_daily_running(r,tend+1:tend2) = (reg_daily_05(r,tend+1:tend2) - reg_05_mean_running(r,y-1987)) ./ reg_05_sd_running(r,y-1987);
        
        if r == 9
            disp([num2str(years(tstart)),' to ',num2str(years(tend)),': Extreme = ',num2str(reg_05_exm_running(r,y-1987))])
        end
        
    end
end

running_exm = max(reg_05_exm_running,[],2);
[running_exm2,running_max_id] = nanmax(reg_05_exm_daily_running,[],2);


% %% As above but with a moving baseline (11 year centred)
% 
% % Find months of maximum temperatures
% [~,max_ids] = max(reg_daily_05,[],2);
% 
% % Create empty arrays for output
% reg_05_sd_running2 = nan(237,32);
% reg_05_mean_running2 = nan(237,32);
% reg_05_exm_running2 = nan(237,32);
% reg_05_exm_daily_running2 = nan(237,length(mon_val));
% 
% 
% % Find long-term SD and mean for calendar month of max for each region
% for r = 1:237
%     
%     % For each region, find which month had max temp
%     max_mon = mon_val(max_ids(r));
%         
%     for y = 1984:2015 % Do a 10-year moving baseline
%         
%         y0 = y-5;
%         y1 = y+5;
%         tstart = min(find(years == y0));
%         tend = max(find(years == y1));
%         
%         
%         % Take just the baseline period (1979-2020, 42 years)
%         reg_daily_baseline = reg_daily_05(:,tstart:tend);
%         mon_baseline = mon_val(tstart:tend);
%  
%         % Select only maximum month
%         mon_ids = mon_baseline == max_mon;
%         
%         % Take SD and mean for that month
%         reg_05_sd_running2(r,y-1983) = std(reg_daily_baseline(r,mon_ids));
%         reg_05_mean_running2(r,y-1983) = mean(reg_daily_baseline(r,mon_ids));
%         
% %         % Calculate extreme magnitude (following Vikki's definition)
% %         reg_05_exm_running(r,y-1983) = (max(reg_daily_05(r,tend+1:tend2)) - reg_05_mean_running(r,y)) ./ reg_05_sd_running(r,y);
% %         
% %         % Go through every day to calculate the metric
% %         reg_05_exm_daily_running(r,tend+1:tend2) = (reg_daily_05(r,tend+1:tend2) - reg_05_mean_running(r,y)) ./ reg_05_sd_running(r,y);
% %         
% %         if r == 9
% %             disp([num2str(years(tstart)),' to ',num2str(years(tend)),': Extreme = ',num2str(reg_05_exm_running(r,y))])
% %         end
%         
%     end
% end
% 
% % running_exm = max(reg_05_exm_running,[],2);
% % [running_exm2,running_max_id] = nanmax(reg_05_exm_daily_running,[],2);


%% Plot
% Load regions to make choropleth map
regs_new_05 = ncread('region_05_regrid.nc','region');
reg_exm_plot = regs_new_05;
reg_exm_dates = regs_new_05;
reg_exm_mons = regs_new_05;
reg_exm_plot_running = regs_new_05;
reg_exm_dates_running = regs_new_05;
reg_exm_mons_running = regs_new_05;

% Assign value to each region
for r = 1:237
    reg_exm_plot(regs_new_05 == r-1) = reg_05_exm(r);
    reg_exm_dates(regs_new_05 == r-1) = years(max_ids(r));
    reg_exm_mons(regs_new_05 == r-1) = mon_val(max_ids(r));
    reg_exm_plot_running(regs_new_05 == r-1) = running_exm(r);
    reg_exm_dates_running(regs_new_05 == r-1) = years(running_max_id(r));
    reg_exm_mons_running(regs_new_05 == r-1) = mon_val(running_max_id(r));

end

% Load lat-long and coast data for plotting
lon_new = double(-179.875:0.25:179.875);
lat_new = double(-90:0.25:90);
[lat,lon]=meshgrid(lat_new,lon_new);

% Coastline
S = shaperead('landareas','UseGeoCoords',true);
coast = load('coast.mat');

% Plot of magnitude of maximum extreme
figure 
colormap(lajolla)

% Set up axes
axesm('MapProjection','Robinson', ...
    'MLineLocation', 30,...
    'PlineLocation', 30, 'MLabelParallel', 'south')

% Plot the data
pcolorm(lat,lon,reg_exm_plot_running)

% Adjust the plot
framem('FEdgeColor', 'black', 'FLineWidth', 1)
gridm('Gcolor',[0.3 0.3 0.3])
tightmap
box off
axis off

% Specify colour limit if necessary
if exist('collim','var')
    caxis(collim)
end

% Add coastline
hold on
geoshow([S.Lat], [S.Lon],'Color','black');

set(gcf, 'color', 'w');
set(gca,'Fontsize',14)

cbar = colorbar;
ylabel(cbar, 'Extreme magnitude (SDs)', 'fontsize', 14)

title('Greatest historic extreme (ERA-5 running mean)')


% Plot of dates
figure 

% Set up axes
axesm('MapProjection','Robinson', ...
    'MLineLocation', 30,...
    'PlineLocation', 30, 'MLabelParallel', 'south')

% Plot the data
pcolorm(lat,lon,reg_exm_dates_running)

% Adjust the plot
framem('FEdgeColor', 'black', 'FLineWidth', 1)
gridm('Gcolor',[0.3 0.3 0.3])
tightmap
box off
axis off

% Specify colour limit if necessary
if exist('collim','var')
    caxis(collim)
end

% Add coastline
hold on
geoshow([S.Lat], [S.Lon],'Color','black');

set(gcf, 'color', 'w');
set(gca,'Fontsize',14)

cbar = colorbar;
ylabel(cbar, 'Extreme magnitude (SDs)', 'fontsize', 14)

title('Year of greatest historic extreme (ERA-5 running mean)')


% %% Extremes analysis -> WWA region
% % Find months of maximum temperatures
% [~,max_ids] = max(reg_daily_WWA);
% 
% % Take just the baseline period (1979-2020, 42 years)
% reg_daily_baseline = reg_daily_WWA(1:15340);
% mon_baseline = mon_val(1:15340);
% 
% % Find which month had max temp
% max_mon = mon_val(max_ids);
% 
% % Select only that month
% mon_ids = mon_baseline == max_mon;
% 
% % Take SD and mean for that month
% reg_WWA_sd = std(reg_daily_baseline(mon_ids));
% reg_WWA_mean = mean(reg_daily_baseline(mon_ids));
% 
% % Calculate extreme magnitude (following Vikki's definition)
% reg_WWA_exm = (max(reg_daily_WWA) - reg_WWA_mean) ./ reg_WWA_sd;
% 
% % Go through every day to calculate the metric
% reg_WWA_exm_daily = (reg_daily_WWA - reg_WWA_mean) ./ reg_WWA_sd;
%     
% % % Save the files to evaluate them offline
% % save('reg_sd.mat','reg_sd')
% % save('reg_mean.mat','reg_mean')
% % save('reg_exm.mat','reg_exm')
% % save('reg_exm_daily.mat','reg_exm_daily')





% %% Count number of events > threshs
% num_events_05 = nan(237,43,5);
% for r = 1:237
%     for n = 1:5
%         starts = 1;
%         ends = 365;
%         for i = 1979:2020
%             if rem(i,4) == 0
%                 yrlen = 366;
%                 
%             else
%                 yrlen = 365;
%             end
%             
%             
%             num_events_05(r,i-1978,n) = sum(reg_05_exm_daily(r,starts:ends)>n);
%             
%             starts = starts+yrlen;
%             ends = ends+yrlen;
%         end
%         
%         num_events_05(r,43,n) = sum(reg_05_exm_daily(r,15342:15539)>n);
% 
%     end
% end
% 
% events05 = squeeze(sum(num_events_05>0,1));
% 
% figure
% 
% x = 1979:2021;
% Y = events05/237 * 100;
% 
% plot(1979:2020,events05(1:42,:)/237 * 100,'linewidth',2)
% set(gca,'Fontsize',14)
% ylabel('Percent of regions experiencing events')
% 
% set(gcf, 'color', 'w');
% set(gca,'Fontsize',14)
% 
% 
% for i = 1:5
%         % Do the regression
%         X = [ones(length(x),1) x'];
%         ensfit1 = X\Y(:,i);
%         
%         % Find the R-squared
%         mdl1 = fitlm(X,Y(:,i));
%         disp([num2str(i),' SD event mean gradient: ',num2str(ensfit1(2)),', R-squared: ',num2str(mdl1.Rsquared.Ordinary)])
% end
% 
% smooth_sum = nan(33,5);
% for i = 1:33
%     vals = i:i+9;
%     smooth_sum(i,:) = mean(events05(vals,:),1);
% end
% 
% hold on
% plot([0 1],[0 1],'-')
% plot([0 1],[0 1],'-')
% plot(1984:2016,smooth_sum/237 * 100,'--','linewidth',2)
% legend('1 SD','2 SD','3 SD','4 SD','5 SD','location','eastoutside')
% xlim([1978 2022])
% 
% hold on
% plot([0 1],[0 1],'-')
% plot([0 1],[0 1],'-')
% plot(2021,events05(43,:)/237 * 100,'o','linewidth',2)
% legend('1 SD','2 SD','3 SD','4 SD','5 SD','location','eastoutside')
% xlim([1978 2022])
% 


%% Count number of events > threshs for running baseline
num_events_05 = nan(237,33,5);
for r = 1:237
    for n = 1:5
        starts = 3653;
        ends = 3652+365;
        for i = 1989:2020
            if rem(i,4) == 0
                yrlen = 366;
                
            else
                yrlen = 365;
            end
            
            
            num_events_05(r,i-1988,n) = sum(reg_05_exm_daily_running(r,starts:ends)>n);
            
            starts = starts+yrlen;
            ends = ends+yrlen;
        end
        
        num_events_05(r,33,n) = sum(reg_05_exm_daily_running(r,15342:15539)>n);

    end
end


events05 = squeeze(sum(num_events_05>0,1));

figure

x = 1989:2021;
Y = events05/237 * 100;


plot(1989:2020,events05(1:32,:)/237 * 100,'linewidth',2)
set(gca,'Fontsize',14)
ylabel('Percent of regions experiencing events')

set(gcf, 'color', 'w');
set(gca,'Fontsize',14)


for i = 1:5
        % Do the regression
        X = [ones(length(x),1) x'];
        ensfit1 = X\Y(:,i);
        
        % Find the R-squared
        mdl1 = fitlm(X,Y(:,i));
        disp([num2str(i),' SD event mean gradient: ',num2str(ensfit1(2)),', R-squared: ',num2str(mdl1.Rsquared.Ordinary)])
end

smooth_sum = nan(23,5);
for i = 1:23
    vals = i:i+9;
    smooth_sum(i,:) = mean(events05(vals,:),1);
end

hold on
plot([0 1],[0 1],'-')
plot([0 1],[0 1],'-')
plot(1994:2016,smooth_sum/237 * 100,'--','linewidth',2)
legend('1 SD','2 SD','3 SD','4 SD','5 SD','location','eastoutside')
xlim([1978 2022])

hold on
plot([0 1],[0 1],'-')
plot([0 1],[0 1],'-')
plot(2021,events05(33,:)/237 * 100,'o','linewidth',2)
legend('1 SD','2 SD','3 SD','4 SD','5 SD','location','eastoutside')
xlim([1988 2022])

