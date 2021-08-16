% Compute regional daily means of ERA-Heat for the NW Pacific heatwave paper.
%
% For reference, this takes approx. 5.5 hours to run and requires 4 or more GB of memory.
%

%% Load region data
frac_05 = ncread('/work/ak0920/region_05_regrid.nc','area_frac');
% Subset the mask to exclude Antarctica
frac_05 = frac_05(:,121:721,:);

%% Load ERA5 data
% Choose which variable to use:
froot = '/bp1store/geog-tropical/data/ERA-Heat/utci/'; % Take the file name...
files = dir([froot '*.nc']); % Then check if any files exist with this root

% If not, then display error
if isempty(files)
    disp('Cannot find any ERA-Heat data.')
else
    % Define new variables to save output to (region x timestep)
    reg_utci_05 = nan(237,length(files));
    mon_val_heat = nan(length(files),1); % This is to record the month number
    day_val_heat = nan(length(files),1); % This is to record the day number
    year_val_heat = nan(length(files),1); % This is to record the year
    
    % Set some values before starting loop
    year = 1979;
    
    % Go through each ERA5 file
    for f = 1:length(files)
        file = [files(f).folder,'/',files(f).name];
        
        % disp(file)
        
        % Load each month of ERA5
        utci = ncread(file,'utci'); % This may need updated if using other heat stress vars
        
        % Find the month number
        mon_num = file(end-15:end-14);
        day_num = file(end-13:end-12);
        year_num = file(end-19:end-16);
        
        % disp('Loaded date info okay')
        
        % Take just the current timestep
        utci_max = nanmax(utci,[],3);
        
        % disp('Calculate daily max okay')
        
        % Rotate the ERA-Heat data to be consistent with regions
        utci_max = flipud(rot90(utci_max,2));
        
        % disp('Corrected orientation okay')

        % Go through each 0.5 Mm2 region to find average
        for r = 1:237
            reg_utci_05(r,f) = nansum(nansum(frac_05(:,:,r).*utci_max)); % Find daily mean
        end
        
        % disp('High res regions okay')
        
        % Record month number
        mon_val_heat(f) = str2double(mon_num);
        day_val_heat(f) = str2double(day_num);
        year_val_heat(f) = str2double(year_num);
        
        % Log the year for sanity check
        if f > 2
            if mon_val_heat(f) == 1 && mon_val_heat(f-1) == 12
                year = year + 1;
                disp(num2str(year))
            end
        end  
    end
    
    % Save the output for quicker loading next time
    save('reg_utci_05.mat','reg_utci_05')
    save('mon_val_heat.mat','mon_val_heat')
    save('day_val_heat.mat','day_val_heat')
    save('year_val_heat.mat','year_val_heat')
end
