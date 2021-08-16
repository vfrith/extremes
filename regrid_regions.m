% Regrid the Daithi Stone regions for use with global ERA-5 data.
% 
% Paper:
% https://escholarship.org/content/qt6v8524c9/qt6v8524c9_noSplash_56a62f3a33ebbd80bf1bc138a6ae30e2.pdf
% Region data as netCDFs:
% https://portal.nersc.gov/c20c/data/C20C/WRAF/All-Hist/est1/v4-1/fx/


% %% Load the raw data for the 2Mm2 regions
% if ~exist('regs_raw','var')
%     regs_raw = double(ncread('region_fx-WRAF2-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.nc','region'));
% end
% 
% 
% %% Regrid
% if ~exist('regs_new','var')
%     
%     % Create old and new grids
%     lon_raw = double(ncread('region_fx-WRAF2-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.nc','lon'));
%     lat_raw = double(ncread('region_fx-WRAF2-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.nc','lat'));
%     [lats_raw,lons_raw]=meshgrid(lat_raw,lon_raw);
%     
%     lon_new = double(-179.875:0.25:179.875);
%     lat_new = double(-90:0.25:90);
%     [lats_new,lons_new]=meshgrid(lat_new,lon_new);
%     
%     regs_new = griddata(lons_raw,lats_raw,regs_raw,lons_new,lats_new,'nearest');
% 
% end
% 
% 
% %% Find area of each region for averaging
% if ~exist('areas_glob','var')
%     areas_glob = calc_latlon_area(lats_new,lons_new,'areaquad');
% end
% 
% % Check if this is consistent with the original dataset:
% areas = double(ncread('region_fx-WRAF2-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.nc','area'));
% 
% areas_new = zeros(68,1);
% for i = 1:68
%     areas_new(i) = nansum(nansum(areas_glob(regs_new == i-1)));
% end
% 
% err = (areas_new-areas)./areas;
% 
% regs_err = regs_new;
% for i = 1:68
%     regs_err(regs_new == i-1) = err(i);
% end
% 
% % Calculate fractional areas for averaging
% regs_frac = repmat(areas_glob,1,1,68);
% 
% for i = 1:68
%     mask = (regs_new == i-1)*1;
%     regs_frac(:,:,i) = (regs_frac(:,:,i).*mask) ./ nansum(nansum(regs_frac(:,:,i).*mask));
% end
% 
% 
% %% Start to save
% 
% fname_long = 'region_regrid.nc';
% Variable = 'region';
% 
% layer = double(ncread('region_fx-WRAF2-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.nc','layer'));
% 
% 
% 
% % Create netCDF and derived variable
% nccreate(fname_long,Variable,'Dimensions',{'lon',length(lon_new),'lat',length(lat_new)},'Datatype','double','Format','netcdf4_classic')
% ncwrite(fname_long,Variable,regs_new);
% ncwriteatt(fname_long,Variable,'standard_name','region_index');
% ncwriteatt(fname_long,Variable,'long_name','Region_Index');
% ncwriteatt(fname_long,Variable,'units','index');
% ncwriteatt(fname_long,Variable,'grid_mapping','rotated_latitude_longitude');
% 
% % Add lat and long data
% nccreate(fname_long,'lat','Dimensions',{'lat',length(lat_new)},'Datatype','double','Format','netcdf4_classic')
% ncwrite(fname_long,'lat',lat_new);
% ncwriteatt(fname_long,'lat','standard_name','latitude');
% ncwriteatt(fname_long,'lat','units','degrees north');
% ncwriteatt(fname_long,'lat','axis','Y');
% 
% nccreate(fname_long,'lon','Dimensions',{'lon',length(lon_new)},'Datatype','double','Format','netcdf4_classic')
% ncwrite(fname_long,'lon',lon_new);
% ncwriteatt(fname_long,'lon','standard_name','longitude');
% ncwriteatt(fname_long,'lon','units','degrees east');
% ncwriteatt(fname_long,'lon','axis','X');
% 
% nccreate(fname_long,'layer','Dimensions',{'layer',68},'Datatype','double','Format','netcdf4_classic')
% ncwrite(fname_long,'layer',layer);
% ncwriteatt(fname_long,'layer','standard_name','layer');
% ncwriteatt(fname_long,'layer','units','integer');
% 
% nccreate(fname_long,'area','Dimensions',{'layer',68},'Datatype','double','Format','netcdf4_classic')
% ncwrite(fname_long,'area',areas);
% ncwriteatt(fname_long,'area','standard_name','region_area');
% ncwriteatt(fname_long,'area','units','km2');
% 
% nccreate(fname_long,'area_frac','Dimensions',{'lon',length(lon_new),'lat',length(lat_new),'layer',68},'Datatype','double','Format','netcdf4_classic')
% ncwrite(fname_long,'area_frac',regs_frac);
% ncwriteatt(fname_long,'area_frac','standard_name','region_area_fraction');
% ncwriteatt(fname_long,'area_frac','units','fraction');
% 
% 
% % Write some general attributes
% %     ncwriteatt(fname_long,'/','collection','HEAT derived variable')
% ncwriteatt(fname_long,'/','creation_date',datestr(now))
% %     ncwriteatt(fname_long,'/','domain',fname(11:12))
% ncwriteatt(fname_long,'/','title','WRAF2 v4.1 regridded')
% 
%   


%% Load the raw data for the 0.5 Mm2 regiond
if ~exist('regs_raw_05','var')
    regs_raw_05 = double(ncread('region_fx-WRAF05-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.nc','region'));
end


%% Regrid
if ~exist('regs_new_05','var')
    
    % Create old and new grids
    lon_raw = double(ncread('region_fx-WRAF05-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.nc','lon'));
    lat_raw = double(ncread('region_fx-WRAF05-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.nc','lat'));
    [lats_raw,lons_raw]=meshgrid(lat_raw,lon_raw);
    
    lon_new = double(-179.875:0.25:179.875);
    lat_new = double(-90:0.25:90);
    [lats_new,lons_new]=meshgrid(lat_new,lon_new);
    
    regs_new_05 = griddata(lons_raw,lats_raw,regs_raw,lons_new,lats_new,'nearest');

end


%% Find area of each region for averaging
if ~exist('areas_glob','var')
    areas_glob = calc_latlon_area(lats_new,lons_new,'areaquad');
end

% Check if this is consistent with the original dataset:
areas = double(ncread('region_fx-WRAF05-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.nc','area'));

areas_new_05 = zeros(237,1);
for i = 1:237
    areas_new_05(i) = nansum(nansum(areas_glob(regs_new_05 == i-1)));
end

err = (areas_new_05-areas)./areas;

regs_err_05 = regs_new_05;
for i = 1:237
    regs_err_05(regs_new_05 == i-1) = err(i);
end

% Calculate fractional areas for averaging
regs_frac_05 = repmat(areas_glob,1,1,237);

for i = 1:237
    mask = (regs_new_05 == i-1)*1;
    regs_frac_05(:,:,i) = (regs_frac_05(:,:,i).*mask) ./ nansum(nansum(regs_frac_05(:,:,i).*mask));
end


%% Start to save

fname_long = 'region_05_regrid.nc';
Variable = 'region';

layer = double(ncread('region_fx-WRAF05-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.nc','layer'));

% Create netCDF and derived variable
nccreate(fname_long,Variable,'Dimensions',{'lon',length(lon_new),'lat',length(lat_new)},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
ncwrite(fname_long,Variable,regs_new_05);
ncwriteatt(fname_long,Variable,'standard_name','region_index');
ncwriteatt(fname_long,Variable,'long_name','Region_Index');
ncwriteatt(fname_long,Variable,'units','index');
ncwriteatt(fname_long,Variable,'grid_mapping','rotated_latitude_longitude');

% Add lat and long data
nccreate(fname_long,'lat','Dimensions',{'lat',length(lat_new)},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
ncwrite(fname_long,'lat',lat_new);
ncwriteatt(fname_long,'lat','standard_name','latitude');
ncwriteatt(fname_long,'lat','units','degrees north');
ncwriteatt(fname_long,'lat','axis','Y');

nccreate(fname_long,'lon','Dimensions',{'lon',length(lon_new)},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
ncwrite(fname_long,'lon',lon_new);
ncwriteatt(fname_long,'lon','standard_name','longitude');
ncwriteatt(fname_long,'lon','units','degrees east');
ncwriteatt(fname_long,'lon','axis','X');

nccreate(fname_long,'layer','Dimensions',{'layer',237},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
ncwrite(fname_long,'layer',layer);
ncwriteatt(fname_long,'layer','standard_name','layer');
ncwriteatt(fname_long,'layer','units','integer');

nccreate(fname_long,'area','Dimensions',{'layer',237},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
ncwrite(fname_long,'area',areas);
ncwriteatt(fname_long,'area','standard_name','region_area');
ncwriteatt(fname_long,'area','units','km2');

nccreate(fname_long,'area_frac','Dimensions',{'lon',length(lon_new),'lat',length(lat_new),'layer',237},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
ncwrite(fname_long,'area_frac',regs_frac_05);
ncwriteatt(fname_long,'area_frac','standard_name','region_area_fraction');
ncwriteatt(fname_long,'area_frac','units','fraction');


% Write some general attributes
ncwriteatt(fname_long,'/','creation_date',datestr(now))
ncwriteatt(fname_long,'/','title','WRAF0.5 v4.1 regridded')

    
%% Area of box from WWA study
% This is only one box, so no need to save as a netCDF
area_WWA = nan(1440,721);
area_WWA(229:244,541:569) = areas_glob(229:244,541:569);

frac_WWA = area_WWA./ nansum(nansum(area_WWA(229:244,541:569)));

save('frac_WWA.mat', 'frac_WWA')


% %% Regrid for CanESM and JRA-55
% 
% % Load the raw data for the 0.5 Mm2 regions
% if ~exist('regs_raw_05','var')
%     regs_raw_05 = double(ncread('region_fx-WRAF05-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.nc','region'));
% end
% 
% Can_lat = ncread('tasmax_day_CanESM5_ssp585_r1i1p1f1_gn_20150101-21001231.nc','lat');
% Can_lon = ncread('tasmax_day_CanESM5_ssp585_r1i1p1f1_gn_20150101-21001231.nc','lon');
% JRA_lat = flipud(ncread('jra-55_6hourly_tas_2021070100_2021073118.nc','lat'));
% JRA_lon = ncread('jra-55_6hourly_tas_2021070100_2021073118.nc','lon');
% 
% 
% 
% %% Regrid
% % if ~exist('regs_Can','var')
%     
%     % Create old and new grids
%     lon_raw = double(ncread('region_fx-WRAF05-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.nc','lon'));
%     lat_raw = double(ncread('region_fx-WRAF05-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.nc','lat'));
%     [lats_raw,lons_raw]=meshgrid(lat_raw,lon_raw);
%     
% 
%     [lats_JRA,lons_JRA]=meshgrid(JRA_lat,JRA_lon);
%     [lats_Can,lons_Can]=meshgrid(Can_lat,Can_lon);
% 
%     lats_JRA2 = lats_JRA;
%     lons_JRA2 = lons_JRA;
%     lats_Can2 = lats_Can;
%     lons_Can2 = lons_Can;
% 
%     lats_JRA2(321:640,:) = lats_JRA(1:320,:);lats_JRA2(1:320,:) = lats_JRA(321:640,:);
%     lons_JRA2(321:640,:) = lons_JRA(1:320,:);lons_JRA2(1:320,:) = lons_JRA(321:640,:)-360;
%     lats_Can2(65:128,:) = lats_Can(1:64,:);lats_Can2(1:64,:) = lats_Can(65:128,:);
%     lons_Can2(65:128,:) = lons_Can(1:64,:);lons_Can2(1:64,:) = lons_Can(65:128,:)-360;
%     
%     
%     regs_JRA = griddata(lons_raw,lats_raw,regs_raw_05,lons_JRA2,lats_JRA2,'nearest');
%     regs_Can = griddata(lons_raw,lats_raw,regs_raw_05,lons_Can2,lats_Can2,'nearest');
% 
% % end
% 
% 
% %% Find area of each region for averaging
% if ~exist('areas_glob','var')
%     areas_glob_Can = calc_latlon_area(lats_Can,lons_Can,'areaquad');
%     areas_glob_JRA = calc_latlon_area(lats_JRA,lons_JRA,'areaquad');
% end
% 
% % Check if this is consistent with the original dataset:
% areas = double(ncread('region_fx-WRAF05-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.nc','area'));
% 
% areas_Can = zeros(237,1);
% areas_JRA = zeros(237,1);
% for i = 1:237
%     areas_Can(i) = nansum(nansum(areas_glob_Can(regs_Can == i-1)));
%     areas_JRA(i) = nansum(nansum(areas_glob_JRA(regs_JRA == i-1)));
% end
% 
% errCan = (areas_Can-areas)./areas;
% errJRA = (areas_JRA-areas)./areas;
% 
% regs_err_Can = regs_Can;
% regs_err_JRA = regs_JRA;
% for i = 1:237
%     regs_err_Can(regs_Can == i-1) = errCan(i);
%     regs_err_JRA(regs_JRA == i-1) = errJRA(i);
% end
% 
% % Calculate fractional areas for averaging
% regs_Can_frac = repmat(areas_glob_Can,1,1,237);
% regs_JRA_frac = repmat(areas_glob_JRA,1,1,237);
% 
% for i = 1:237
%     mask = (regs_Can == i-1)*1;
%     regs_Can_frac(:,:,i) = (regs_Can_frac(:,:,i).*mask) ./ nansum(nansum(regs_Can_frac(:,:,i).*mask));
%     mask = (regs_JRA == i-1)*1;
%     regs_JRA_frac(:,:,i) = (regs_JRA_frac(:,:,i).*mask) ./ nansum(nansum(regs_JRA_frac(:,:,i).*mask));
% end
% 
% 
% %% Start to save CanESM5 grid
% fname_long = 'region_Can_regrid.nc';
% Variable = 'region';
% 
% layer = double(ncread('region_fx-WRAF05-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.nc','layer'));
% 
% 
% 
% % Create netCDF and derived variable
% nccreate(fname_long,Variable,'Dimensions',{'lon',length(Can_lon),'lat',length(Can_lat)},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
% ncwrite(fname_long,Variable,regs_Can);
% ncwriteatt(fname_long,Variable,'standard_name','region_index');
% ncwriteatt(fname_long,Variable,'long_name','Region_Index');
% ncwriteatt(fname_long,Variable,'units','index');
% ncwriteatt(fname_long,Variable,'grid_mapping','rotated_latitude_longitude');
% 
% % Add lat and long data
% nccreate(fname_long,'lat','Dimensions',{'lat',length(Can_lat)},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
% ncwrite(fname_long,'lat',Can_lat);
% ncwriteatt(fname_long,'lat','standard_name','latitude');
% ncwriteatt(fname_long,'lat','units','degrees north');
% ncwriteatt(fname_long,'lat','axis','Y');
% 
% nccreate(fname_long,'lon','Dimensions',{'lon',length(Can_lon)},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
% ncwrite(fname_long,'lon',Can_lon);
% ncwriteatt(fname_long,'lon','standard_name','longitude');
% ncwriteatt(fname_long,'lon','units','degrees east');
% ncwriteatt(fname_long,'lon','axis','X');
% 
% nccreate(fname_long,'layer','Dimensions',{'layer',237},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
% ncwrite(fname_long,'layer',layer);
% ncwriteatt(fname_long,'layer','standard_name','layer');
% ncwriteatt(fname_long,'layer','units','integer');
% 
% nccreate(fname_long,'area','Dimensions',{'layer',237},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
% ncwrite(fname_long,'area',areas_Can);
% ncwriteatt(fname_long,'area','standard_name','region_area');
% ncwriteatt(fname_long,'area','units','km2');
% 
% nccreate(fname_long,'area_frac','Dimensions',{'lon',length(Can_lon),'lat',length(Can_lon),'layer',237},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
% ncwrite(fname_long,'area_frac',regs_Can_frac);
% ncwriteatt(fname_long,'area_frac','standard_name','region_area_fraction');
% ncwriteatt(fname_long,'area_frac','units','fraction');
% 
% 
% % Write some general attributes
% %     ncwriteatt(fname_long,'/','collection','HEAT derived variable')
% ncwriteatt(fname_long,'/','creation_date',datestr(now))
% %     ncwriteatt(fname_long,'/','domain',fname(11:12))
% ncwriteatt(fname_long,'/','title','WRAF05 v4.1 regridded to CanESM5 grid')
% 
% 
% %% Start to save JRA-55 grid
% % 
% % 
% % 
% 
% 
%  fname_long = 'region_JRA_regrid.nc';
% Variable = 'region';
% 
% layer = double(ncread('region_fx-WRAF05-v4-1_WRAF_All-Hist_est1_v4-1_4-1-0_000000-000000.nc','layer'));
% 
% 
% 
% % Create netCDF and derived variable
% nccreate(fname_long,Variable,'Dimensions',{'lon',length(JRA_lon),'lat',length(JRA_lat)},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
% ncwrite(fname_long,Variable,regs_JRA);
% ncwriteatt(fname_long,Variable,'standard_name','region_index');
% ncwriteatt(fname_long,Variable,'long_name','Region_Index');
% ncwriteatt(fname_long,Variable,'units','index');
% ncwriteatt(fname_long,Variable,'grid_mapping','rotated_latitude_longitude');
% 
% % Add lat and long data
% nccreate(fname_long,'lat','Dimensions',{'lat',length(JRA_lat)},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
% ncwrite(fname_long,'lat',JRA_lat);
% ncwriteatt(fname_long,'lat','standard_name','latitude');
% ncwriteatt(fname_long,'lat','units','degrees north');
% ncwriteatt(fname_long,'lat','axis','Y');
% 
% nccreate(fname_long,'lon','Dimensions',{'lon',length(JRA_lon)},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
% ncwrite(fname_long,'lon',JRA_lon);
% ncwriteatt(fname_long,'lon','standard_name','longitude');
% ncwriteatt(fname_long,'lon','units','degrees east');
% ncwriteatt(fname_long,'lon','axis','X');
% 
% nccreate(fname_long,'layer','Dimensions',{'layer',237},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
% ncwrite(fname_long,'layer',layer);
% ncwriteatt(fname_long,'layer','standard_name','layer');
% ncwriteatt(fname_long,'layer','units','integer');
% 
% nccreate(fname_long,'area','Dimensions',{'layer',237},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
% ncwrite(fname_long,'area',areas_JRA);
% ncwriteatt(fname_long,'area','standard_name','region_area');
% ncwriteatt(fname_long,'area','units','km2');
% 
% nccreate(fname_long,'area_frac','Dimensions',{'lon',length(JRA_lon),'lat',length(JRA_lon),'layer',237},'Datatype','double','Format','netcdf4_classic','DeflateLevel',2)
% ncwrite(fname_long,'area_frac',regs_JRA_frac);
% ncwriteatt(fname_long,'area_frac','standard_name','region_area_fraction');
% ncwriteatt(fname_long,'area_frac','units','fraction');
% 
% 
% % Write some general attributes
% %     ncwriteatt(fname_long,'/','collection','HEAT derived variable')
% ncwriteatt(fname_long,'/','creation_date',datestr(now))
% %     ncwriteatt(fname_long,'/','domain',fname(11:12))
% ncwriteatt(fname_long,'/','title','WRAF05 v4.1 regridded to JRA-55 grid')

