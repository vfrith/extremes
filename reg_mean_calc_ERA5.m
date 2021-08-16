% Compute regional daily means of ERA5 for the NW Pacific heatwave paper.

%% Load region data
frac_2 = ncread('/work/ak0920/region_regrid.nc','area_frac');
frac_05 = ncread('/work/ak0920/region_05_regrid.nc','area_frac');
load('/work/ak0920/frac_WWA.mat');

%% Load ERA5 data
% Choose which variable to use:
% froot = '/bp1store/geog-tropical/data/ERA-5/day/ta/tas/'; % Take the file name...
froot = '/bp1store/geog-tropical/data/ERA-5/day/tasmax/'; % Take the file name...
files = dir([froot '*.nc']); % Then check if any files exist with this root

% If not, then display error
if isempty(files)
    disp('Cannot find any ERA5 data.')
else
    % Define new variables to save output to (region x timestep)
    reg_daily_2 = nan(68,15539); % Time series is at least this long
    reg_daily_05 = nan(237,15539); % Time series is at least this long
    reg_daily_WWA = nan(15539,1); % Time series is at least this long
    mon_val = nan(15539,1); % This is to record the month number
    
    % Set some values before starting loop
    count = 1;
    year = 1979;
    
    % Go through each ERA5 file
    for f = 1:length(files)
        file = [files(f).folder,'/',files(f).name];
                
        % Load each month of ERA5
        tas = ncread(file,'t2m'); % This may need updated if using other heat stress vars
        
        % Find the month number
        mon_num = file(end-4:end-3);
                
        % Go through each day of current month
        for d = 1:length(tas(1,1,:))
            
            % Take just the current timestep
            tas_day = tas(:,:,d);
                        
            % Reshape the ERA5 data to be consistent with regions
            tas_day = flipud(rot90(tas_day,2));
            tas_day2 = tas_day;
            tas_day(1:720,:) = tas_day2(721:1440,:);
            tas_day(721:1440,:) = tas_day2(1:720,:);
                        
            % Go through each 2 Mm2 region to find average
            for r = 1:68
                reg_daily_2(r,count) = nansum(nansum(frac_2(:,:,r).*tas_day)); % Find daily mean
            end
                        
            % Go through each 0.5 Mm2 region to find average
            for r = 1:237
                reg_daily_05(r,count) = nansum(nansum(frac_05(:,:,r).*tas_day)); % Find daily mean
            end
            
            % Finally, find average for WWA region
            reg_daily_WWA(count) = nansum(nansum(frac_WWA.*tas_day)); % Find daily mean
            
            % Record month number
            mon_val(count) = str2double(mon_num);
            
            % Log the year for sanity check
            if count > 1
                if mon_val(count) == 1 && mon_val(count-1) == 12
                    year = year + 1;
                    disp(num2str(year))
                end
            end
            
            count = count + 1; % Step through for next daily timestep
            
        end
    end
    
    % Save the output for quicker loading next time
    save('reg_daily_2.mat','reg_daily_2')
    save('reg_daily_05.mat','reg_daily_05')
    save('reg_daily_WWA.mat','reg_daily_WWA')
    save('mon_val.mat','mon_val')
end
