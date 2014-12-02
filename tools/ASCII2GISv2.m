function [VelOut,goodrows] = ASCII2GISv2(drange,ref,tav)
% WinRiver ASCII to GIS Import Format 

% This program reads in an ASCII file or files generated from WinRiver 
% Classic ASCII output and outputs the Georeferenced mean velocity,
% temperature, depth, and backscatter data to a file for input into GIS.

% drange = [depth1 depth2] %range of depths over which to average the data
% ('dfs' option)
% (full range of data if blank)  %Added 12-20-10

% drange = single value for 'hab' option (height above bottom in m)

% ref = 'dfs' or 'hab';  %'dsf' = depth from surface; 'hab' = height above
% bottom

% tav = averaging time in seconds (leave empty for no averaging)

%Updated directional averaging PRJ 2/8/11
%Updated save path PRJ 3/10/11
%Added *.anv file export PRJ 5-11-11
%Added averaging capability PRJ 3-20-12
%Added ability to plot profiles along path (v2) 9-10-14
%Made conversion to work with VMT package [FLE 11-26-2014]

% P.R. Jackson 6-25-10


%% USer inputs
append_data = 1;
comp_us = 1;        % Compute shear velocity
plot_profiles = 0;  %Plots median profiles for each averaging timestep
plot_projected = 0; % Turn on projected velocity for plotting profiles
proj_dir = 144.4;   % projection direction for velocity projection (profiles)
plot_english = 0;
nma = 1;            % half window size for moving average of vertical profiles

if isempty(tav)
    avg_data = 0;
else
    avg_data = 1;
end

%% Check the inputs

if nargin == 0
    drange = [];
    ref = 'dfs';
elseif nargin < 2
    ref = 'dfs';
end

%% Read and Convert the Data

% Determine Files to Process
% Ask the user to select files:
% -----------------------------
guiprefs = getpref('VMT');
current_file = fullfile(guiprefs.ascii.path,guiprefs.ascii.file{1});
[zFileName,zPathName] = uigetfile({'*_ASC.TXT','ASCII (*_ASC.TXT)'}, ...
    'Select the ASCII Output Files', ...
    current_file, ...
    'MultiSelect','on');

if ischar(zPathName) % The user did not hit "Cancel"
% Determine number of files to be processed
if  isa(zFileName,'cell')
    z=size(zFileName,2);
    zFileName=sort(zFileName);       
else
    z=1;
    zFileName={zFileName};
end
%msgbox('Loading Data...Please Be Patient','Conversion Status','help','replace');
% Read in Selected Files
% Initialize row counter

% Query for an output file name and location
    [ofile,opath] = uiputfile('*.csv','Save file name',zPathName);
    outfile = [opath ofile];

% Begin master loop

VelOut = [];  %Matrix for output of velocity data

wbh = waitbar(0,'Converting Data Files...Please Be Patient');

for zi=1:z
    %Clear variables
    clear DAVeast DAVnorth DAVmag DAVdir DAVvert ustar zo cod i j
    
    % Open txt data file
    if  isa(zFileName,'cell')
        fileName=fullfile(zPathName,zFileName{zi});
        fileName=char(fileName);
    else
        fileName=strcat(zPathName,zFileName);
    end

    % screenData 0 leaves bad data as -32768, 1 converts to NaN
    screenData=1;

    % tfile reads the data from an RDI ASCII output file and puts the
    % data in a Matlab data structure with major groups of:
    % Sup - supporing data
    % Wat - water data
    % Nav - navigation data including GPS.
    % Sensor - Sensor data
    % Q - discharge related data
	ignoreBS = 0;
    [A]=tfile(fileName,screenData,ignoreBS);
    %Extract only Lat lon data
    latlon(:,1)=A.Nav.lat_deg(:,:);
    latlon(:,2)=A.Nav.long_deg(:,:);
    
    %Rescreen data to remove missing data (30000 value)
    indx1 = find(abs(latlon(:,1)) > 90);
    indx2 = find(abs(latlon(:,2)) > 180);
    latlon(indx1,1)=NaN;
    latlon(indx2,2)=NaN;
    
    indx3 = find(~isnan(latlon(:,1)) & ~isnan(latlon(:,2)));
    latlon = latlon(indx3,:);
    
    
    %Extract the Depths
    BeamDepths  = A.Nav.depth(indx3,:);
    Depth = nanmean(A.Nav.depth(indx3,:),2);
    
    %Filter Backscatter 
    A = VMT_FilterBS(1,A);
    
    
    %Extract the averaged velocities and backscatter (layer average)
    if isempty(drange)  
        disp(['Extracting DFS Range = Full'])
        DAVeast  = VMT_LayerAveMean(A.Wat.binDepth(:,indx3),A.Wat.vEast(:,indx3));
        DAVnorth = VMT_LayerAveMean(A.Wat.binDepth(:,indx3),A.Wat.vNorth(:,indx3));
        DAVvert  = VMT_LayerAveMean(A.Wat.binDepth(:,indx3),A.Wat.vVert(:,indx3));
        DABack   = VMT_LayerAveMean(A.Wat.binDepth(:,indx3),A.Clean.bsf(:,indx3));
        %DAVeast  = nanmean(A.Wat.vEast(:,indx3),1)';
        %DAVnorth = nanmean(A.Wat.vNorth(:,indx3),1)';
        %DAVvert  = nanmean(A.Wat.vVert(:,indx3),1)';
        %DABack   = nanmean(A.Clean.bsf(:,indx3),1)';
        DAVeast  = DAVeast';
        DAVnorth = DAVnorth';
        DAVvert  = DAVvert';
        DABack   = DABack';
    elseif strcmp(ref,'dfs')
        disp(['Extracting DFS Range = ' num2str(drange(1)) ' to ' num2str(drange(2)) ' m'])
        indxr = find(A.Wat.binDepth(:,1) >= drange(1) & A.Wat.binDepth(:,1) <= drange(2));
        DAVeast  = VMT_LayerAveMean(A.Wat.binDepth(indxr,indx3),A.Wat.vEast(indxr,indx3));
        DAVnorth = VMT_LayerAveMean(A.Wat.binDepth(indxr,indx3),A.Wat.vNorth(indxr,indx3));
        DAVvert  = VMT_LayerAveMean(A.Wat.binDepth(indxr,indx3),A.Wat.vVert(indxr,indx3));
        DABack   = VMT_LayerAveMean(A.Wat.binDepth(indxr,indx3),A.Clean.bsf(indxr,indx3));        
        %DAVeast  = nanmean(A.Wat.vEast(indxr,indx3),1)';
        %DAVnorth = nanmean(A.Wat.vNorth(indxr,indx3),1)';
        %DAVvert  = nanmean(A.Wat.vVert(indxr,indx3),1)';
        %DABack   = nanmean(A.Clean.bsf(indxr,indx3),1)';
        DAVeast  = DAVeast';
        DAVnorth = DAVnorth';
        DAVvert  = DAVvert';
        DABack   = DABack';
    elseif strcmp(ref,'hab')
        disp(['Extracting HAB Limit = ' num2str(drange) ' m'])
        i = 1;
        for j = 1:length(indx3)
            bed = nanmean(A.Nav.depth(indx3(j),:),2)';
            indxr = find(A.Wat.binDepth(:,1) >= (bed - drange(1)) & A.Wat.binDepth(:,1) <= bed);
%             DAVeast(i)  = VMT_LayerAveMean(A.Wat.binDepth(indxr,indx3(j)),A.Wat.vEast(indxr,indx3(j)));
%             DAVnorth(i) = VMT_LayerAveMean(A.Wat.binDepth(indxr,indx3(j)),A.Wat.vNorth(indxr,indx3(j)));
%             DAVvert(i)  = VMT_LayerAveMean(A.Wat.binDepth(indxr,indx3(j)),A.Wat.vVert(indxr,indx3(j)));
%             DABack(i)   = VMT_LayerAveMean(A.Wat.binDepth(indxr,indx3(j)),A.Clean.bsf(indxr,indx3(j)));
            DAVeast(i)  = nanmean(A.Wat.vEast(indxr,indx3(j)),1);
            DAVnorth(i) = nanmean(A.Wat.vNorth(indxr,indx3(j)),1);
            DAVvert(i)  = nanmean(A.Wat.vVert(indxr,indx3(j)),1);
            DABack(i)   = nanmean(A.Clean.bsf(indxr,indx3(j)),1)';
            
            i = i + 1;
        end
        
        DAVeast  = DAVeast';
        DAVnorth = DAVnorth';
        DAVvert  = DAVvert';
        DABack   = DABack';
    end
    
    % Compute the magnitude from the components
    DAVmag   = sqrt(DAVeast.^2 + DAVnorth.^2);
    
    % Compute the average direction from the velocity components
    DAVdir = 90 - (atan2(DAVnorth, DAVeast))*180/pi; %Compute the atan from the velocity componentes, convert to radians, and rotate to north axis
    qindx = find(DAVdir < 0);
    if ~isempty(qindx)
        DAVdir(qindx) = DAVdir(qindx) + 360;  %Must add 360 deg to Quadrant 4 values as they are negative angles from the +y axis
    end
        
    %Extract the Sensors
    Pitch = A.Sensor.pitch_deg(indx3);
    Roll  = A.Sensor.roll_deg(indx3);
    Heading  = A.Sensor.heading_deg(indx3);
    Temp  = A.Sensor.temp_degC(indx3);
    
    %Extract the time stamps
    MTyear      = A.Sup.year(indx3) + 2000; %works for data file after the year 2000
    MTmon       = A.Sup.month(indx3);
    MTday       = A.Sup.day(indx3);
    MThour      = A.Sup.hour(indx3);
    MTmin       = A.Sup.minute(indx3);
    MTsec       = A.Sup.second(indx3) + A.Sup.sec100(indx3)/100;
    MTdatenum   = datenum([MTyear MTmon MTday MThour MTmin MTsec]);

    %Extract Ens No
    EnsNo = A.Sup.ensNo(indx3);
    

    if comp_us %Compute normalized, bed origin profiles to prepare for log law fitting (PRJ, 8-31-12)
        d_ens   = nanmean(A.Nav.depth(indx3,:),2)';  %Average depth from the four beams for every ensemble
        z_bins  = repmat(d_ens,size(A.Wat.binDepth(:,indx3),1),1) - A.Wat.binDepth(:,indx3);  %matrix on bin depths ref to bottom
        z_norm  = z_bins./repmat(d_ens,size(A.Wat.binDepth(:,indx3),1),1);  %Matrix of normalized, bed origin bin depths     
    end 
        
        
    if 1  %Fit individual profiles to log law
        clear i j
        i = 1;
        for j = 1:length(indx3)
            dfit = nanmean(A.Nav.depth(indx3(j),:),2);
            zfit = dfit - A.Wat.binDepth(:,1);
            znorm = zfit./dfit;
            indxfr = find(znorm >= 0.2 & znorm <= 1); %takes only data above 0.2H
            ufit = A.Wat.vMag(indxfr,indx3(j))/100;
            zfit = zfit(indxfr);
            indxgd = find(~isnan(ufit));
            if ~isempty(indxgd)
                [ustar(i),zo(i),cod(i)] = fitLogLawV2(ufit(indxgd),zfit(indxgd),dfit);
                if cod(i) < 0.25 | ustar(i) < 0 | zo(i) > 1.0  %screens the results
                    ustar(i) = nan;
                    zo(i) = nan;
                end
            else
                ustar(i) = nan;
                zo(i) = nan;
                cod(i) = nan;
            end
            i = i + 1;
        end
        ustar = ustar';
        zo = zo';
        cod = cod';
    else % Fill with nans if not computing (turn off pending more testing--PRJ 6-30-11)
        ustar = nan.*ones(size(EnsNo));
        zo  = nan.*ones(size(EnsNo));
        cod = nan.*ones(size(EnsNo));
    end
    
    
    % Perform temporal averaging  (Added 3-20-12 PRJ)
    if avg_data
        disp(['Performing temporal averaging over ' num2str(tav) ' second intervals.'])
        %tav = 30; %Averaging time in seconds
        if (MTdatenum(1) + tav/(3600*24)) >= MTdatenum(end)  %uses limits of data if averaging period exceeds data limits
            tav_vec = [MTdatenum(1) MTdatenum(end)];
        else
            tav_vec = MTdatenum(1):(tav/(3600*24)):MTdatenum(end);  %Vector of serial dates representing the start and end of each averaging period
        end
        
        %Preallocate
        dnm = nan*ones(size(A.Wat.binDepth,1),length(tav_vec)-1);
        vpm = nan*ones(size(dnm));
        dUdz = nan*ones(size(A.Wat.binDepth,1)-1,length(tav_vec)-1);
        dUdz_z = nan*ones(size(dUdz));
        velprof = nan*ones(size(dnm));
        
        for i = 1:length(tav_vec)-1
            av_indx = find(MTdatenum >= tav_vec(i) & MTdatenum < tav_vec(i+1));
            numavg(i) = length(av_indx);
            EnsNo_av(i) = nanmean(ceil(EnsNo(av_indx)));
            MTdatenum_av(i) = nanmean(MTdatenum(av_indx));
            latlon_av(i,:) = nanmean(latlon(av_indx,:),1);
            Heading_av(i) = nanmean(Heading(av_indx));  %this will break down near due north
            Pitch_av(i) = nanmean(Pitch(av_indx));
            Roll_av(i) = nanmean(Roll(av_indx));
            Temp_av(i) = nanmean(Temp(av_indx));
            Depth_av(i) = nanmean(Depth(av_indx));
            BeamDepths_av(i,:) = nanmean(BeamDepths(av_indx,:),1);
            DABack_av(i) = nanmean(DABack(av_indx));
            DAVeast_av(i) = nanmean(DAVeast(av_indx));
            DAVnorth_av(i) = nanmean(DAVnorth(av_indx));
            DAVvert_av(i) = nanmean(DAVvert(av_indx));
            
            if comp_us  %Compute the shear velocity
                %Compute the mean, normalized profile (bed origin)
                %[znm,vm] = VMT_ComputeNormProf(z_norm(:,av_indx),A.Wat.vMag(:,av_indx),30);
                [binDepth,znm,Vme,Vmn,Vmv,vm,obsav,maxdepth(i)] = ComputeNormalizedProfile(nanmean(A.Nav.depth(av_indx,:),2)',A.Wat.binDepth(:,av_indx),A.Wat.vEast(:,av_indx),A.Wat.vNorth(:,av_indx),A.Wat.vVert(:,av_indx));
                znm = znm';
                binDepth = binDepth';
                vm_plot{i} = vm;
                znm_plot{i} = znm;
                binDepth_plot{i} = binDepth;
                
                %Compute the mean profile (surface origin)
                %[dnm(:,i),vpm(:,i)] = VMT_ComputeProf(A.Wat.binDepth(:,av_indx),A.Wat.vMag(:,av_indx),Depth_av(i));
                [dnm(:,i),vem,vnm,vvm,vpm(:,i),maxd(i),vdir(:,i)] = VMT_ComputeProf(A.Wat.binDepth(:,av_indx),A.Wat.vEast(:,av_indx),A.Wat.vNorth(:,av_indx),A.Wat.vVert(:,av_indx),nanmean(A.Nav.depth(av_indx,:),2)');
                
                % Determine the projected velocity (specified direction)
                psi = (vdir(:,i)-proj_dir);

                % Determine the projected velocity (U) and transverse (V)
                U(:,i) = cosd(psi).*vpm(:,i);
                V(:,i) = sind(psi).*vpm(:,i);
                                        
                %Fit the normalized profile with the log law
                gd_indx = ~isnan(vm);
                u_fit = vm(gd_indx)./100;
                z_fit = maxdepth(i) - binDepth(gd_indx); %znm(gd_indx)*maxdepth(i);               
                [ustar_av(i),zo_av(i),cod_av(i)] = fitLogLawV2(u_fit,z_fit,maxdepth(i));
            else
                ustar_av(i) = nanmean(ustar(av_indx));
                zo_av(i) = nanmean(zo(av_indx));
                cod_av(i) = nanmean(cod(av_indx));
                maxd(i) = nanmax(Depth(av_indx));
            end        
        end
        
        % Compute the magnitude and direction from the averaged
        % components

        DAVmag_av = sqrt(DAVeast_av.^2 + DAVnorth_av.^2);
        DAVdir_av = 90 - (atan2(DAVnorth_av, DAVeast_av))*180/pi; %Compute the atan from the velocity componentes, convert to radians, and rotate to north axis
        qindx = find(DAVdir_av < 0);
        if ~isempty(qindx)
            DAVdir_av(qindx) = DAVdir_av(qindx) + 360;  %Must add 360 deg to Quadrant 4 values as they are negative angles from the +y axis
        end
    else % No data averaging requested
        disp(['No spatial averaging requested. Processing all data.'])
           
        %Preallocate
        dnm = nan(size(A.Wat.binDepth,2),1);
        vpm = nan*ones(size(dnm));
        dUdz = nan*ones(size(A.Wat.binDepth,2)-1,1);
        dUdz_z = nan*ones(size(dUdz));
        velprof = nan*ones(size(dnm));
        
        for i = 1:length(dnm)
           
            if comp_us  %Compute the shear velocity
                dbstop if error
                %Compute the mean, normalized profile (bed origin)
                %[znm,vm] = VMT_ComputeNormProf(z_norm(:,av_indx),A.Wat.vMag(:,av_indx),30);
                [binDepth,znm,Vme,Vmn,Vmv,vm,obsav,maxdepth(i)] = ComputeNormalizedProfile(nanmean(A.Nav.depth(i,:),2)',A.Wat.binDepth(:,i),A.Wat.vEast(:,i),A.Wat.vNorth(:,i),A.Wat.vVert(:,i));
                znm = znm';
                binDepth = binDepth';
                vm_plot{i} = vm;
                znm_plot{i} = znm;
                binDepth_plot{i} = binDepth;
                
                %Compute the mean profile (surface origin)
                %[dnm(:,i),vpm(:,i)] = VMT_ComputeProf(A.Wat.binDepth(:,av_indx),A.Wat.vMag(:,av_indx),Depth_av(i));
                [dnm(:,i),vem,vnm,vvm,vpm(:,i),maxd(i),vdir(:,i)] = VMT_ComputeProf(A.Wat.binDepth(:,i),A.Wat.vEast(:,i),A.Wat.vNorth(:,i),A.Wat.vVert(:,i),nanmean(A.Nav.depth(i,:),2)');
                
                % Determine the projected velocity (specified direction)
                psi = (vdir(:,i)-proj_dir);

                % Determine the projected velocity (U) and transverse (V)
                U(:,i) = cosd(psi).*vpm(:,i);
                V(:,i) = sind(psi).*vpm(:,i);
                                        
                %Fit the normalized profile with the log law
                gd_indx = ~isnan(vm);
                u_fit = vm(gd_indx)./100;
                z_fit = maxdepth(i) - binDepth(gd_indx); %znm(gd_indx)*maxdepth(i);               
                [ustar_av(i),zo_av(i),cod_av(i)] = fitLogLawV2(u_fit,z_fit,maxdepth(i));
            else
                ustar_av(i) = nanmean(ustar(i));
                zo_av(i) = nanmean(zo(i));
                cod_av(i) = nanmean(cod(i));
                maxd(i) = nanmax(Depth(i));
            end        
        end
        
        % Compute the magnitude and direction from the averaged
        % components

        DAVmag_av = sqrt(DAVeast.^2 + DAVnorth.^2);
        DAVdir_av = 90 - (atan2(DAVnorth, DAVeast))*180/pi; %Compute the atan from the velocity componentes, convert to radians, and rotate to north axis
        qindx = find(DAVdir_av < 0);
        if ~isempty(qindx)
            DAVdir_av(qindx) = DAVdir_av(qindx) + 360;  %Must add 360 deg to Quadrant 4 values as they are negative angles from the +y axis
        end
        EnsNo_av = EnsNo;
    end 
    
    if plot_projected
        ylimits = [0 max(maxd)];
        xlimits = [nanmin(nanmin(U,[],1)) nanmax(nanmax(U,[],1))];
    else
        ylimits = [0 max(maxd)];
        xlimits = [0 nanmax(nanmax(vpm,[],1))];
    end
    if plot_english
        xlimits = xlimits*0.03281;  %cm/s to ft/s
        ylimits = ylimits*3.281; %m to ft
    end
    
    if avg_data
        lldist = m_lldist(latlon_av(:,2),latlon_av(:,1));
        dmg = 1000*cumsum([0; lldist]); %distance made good in meters
    else
        lldist = m_lldist(latlon(:,2),latlon(:,1));
        dmg = 1000*cumsum([0; lldist]); %distance made good in meters
    end
    
    if plot_profiles
        % Plot the profiles as panels with DMG titles
        nprofs = length(EnsNo_av);
        if nprofs > 200
            hwarn = warndlg('Averaging timstep too small.  Too many profiles to plot.  Increase averaging time or turn off profile plotting option.');
            return
        end
        if nprofs > 10 & nprofs < 20
            Nw = 10;
            Nh = ceil(nprofs/Nw);
        elseif nprofs > 20 & nprofs <= 30
            Nw = 15;
            Nh = ceil(nprofs/Nw);
        elseif nprofs > 30
            Nw = 20;
            Nh = ceil(nprofs/Nw);
        elseif nprofs < 10
            Nw = nprofs;
            Nh = 1;
        end
        h3 = figure(3); clf
        maximize(h3)
        ha = tight_subplot(Nh, Nw, [0.05 0.02], [.05 .05],[.05 .05]);
        for i = 1:nprofs
            if plot_english
                if plot_projected
                    velprof(:,i) = nanmoving_average(U(:,i)*0.03281,nma);
                    axes(ha(i)); plot(velprof(:,i), dnm(:,i)*3.281,'ko-','MarkerFaceColor','k','MarkerSize',2); hold on
                else
                    velprof(:,i) = nanmoving_average(vpm(:,i)*0.03281,nma);
                    axes(ha(i)); plot(velprof(:,i), dnm(:,i)*3.281,'ko-','MarkerFaceColor','k','MarkerSize',2); hold on
                    plot(DAVmag_av(i)*3.281, 0.0,'ro','MarkerFaceColor','r','MarkerSize',2); hold on
                end
                p1 = patch([xlimits(1) xlimits(1) xlimits(2) xlimits(2)],[ylimits(2) maxd(i)*3.281 maxd(i)*3.281 ylimits(2)], [0.7383    0.7148    0.4180],'LineStyle','-','EdgeColor','k','FaceAlpha',0.3); hold on
                %uistack2(p1,'bottom');
                plot(xlimits,[maxd(i)*3.281 maxd(i)*3.281],'k-'); hold on
                plot([0 0], [ylimits(1) maxd(i)*3.281],'k-'); hold on
                %title(['DMG = ' num2str(dmg(i)*3.281/5280,'%3.3f') ' mi'],'FontSize',8)
                title({['DMG = ' num2str(dmg(i)/1000,'%3.3f') ' km']; ['n = ' num2str(numavg(i))]},'FontSize',8)
                set(gca,'FontSize',10)
                box on
            else
                if plot_projected
                    velprof(:,i) = nanmoving_average(U(:,i),nma);
                    axes(ha(i)); plot(velprof(:,i), dnm(:,i),'ko-','MarkerFaceColor','k','MarkerSize',2); hold on
                else
                    velprof(:,i) = nanmoving_average(vpm(:,i),nma);
                    axes(ha(i)); plot(velprof(:,i), dnm(:,i),'ko-','MarkerFaceColor','k','MarkerSize',2); hold on
                    plot(DAVmag_av(i), 0.0,'ro','MarkerFaceColor','r','MarkerSize',2); hold on
                end
                %Compute dU/dz
                dUdz(:,i) = diff(velprof(:,i)/100)./diff(dnm(:,i));
                %Compute the depths for the center of the dUdz bins
                dUdz_z(:,i) = dnm(1:end-1,i) + (dnm(2:end,i) - dnm(1:end-1,i))/2;
                plot(xlimits,[maxd(i) maxd(i)],'k-'); hold on
                p1 = patch([xlimits(1) xlimits(1) xlimits(2) xlimits(2)],[ylimits(2) maxd(i) maxd(i) ylimits(2)], [0.7383    0.7148    0.4180],'LineStyle','-','EdgeColor','k','FaceAlpha',0.3); hold on
                uistack2(p1,'bottom');
                plot([0 0], [ylimits(1) maxd(i)],'k-'); hold on
                title(['DMG = ' num2str(dmg(i)/1000,'%3.3f') ' km'],'FontSize',8)
                set(gca,'FontSize',10)
                box on
            end
            set(gca,'YDir','reverse')
            xlim(xlimits)
            %ylim([0.0 ceil(maxd(i))])  %separate y scales
            ylim(ylimits)  %same y scales
            grid on
        end
        for j = i+1:length(ha)
            set(ha(j),'Visible','off')
        end
        
        %Save the dUdz profiles
        if 1
            save([outfile(1:end-4) '_dUdzProf.mat'], 'dUdz_z', 'dUdz','velprof','dnm','maxd','latlon_av');
        end
               
        %Plot the bottom-refrenced profiles with log law fits
        %Find the limits
        max_x = 0;
        for i = 1:nprofs
            max_x = max([max_x max(vm_plot{i})]);
        end

        h4 = figure(4); clf
        maximize(h4)
        haa = tight_subplot(Nh, Nw, [0.05 0.02], [.05 .05],[.05 .05]);
        for i = 1:nprofs
            axes(haa(i)); plot(vm_plot{i}, maxdepth(i) - binDepth_plot{i},'ko-','MarkerFaceColor','k','MarkerSize',2); hold on
            plot(DAVmag_av(i), 0.0,'ro','MarkerFaceColor','r','MarkerSize',2); hold on
            zeval = linspace(0,maxdepth(i));
            ueval = ustar_av(i)./0.41.*log(zeval./zo_av(i));
            plot(ueval*100,zeval,'r-'); hold on
            dum_x = get(haa(i),'XLim');
            plot([0 dum_x(2)],[maxdepth(i) maxdepth(i)],'k-'); hold on
            title({['DMG = ' num2str(dmg(i)/1000,'%3.3f') ' km']; ['u* = ' num2str(ustar_av(i),'%0.3f') ' m/s; COD = ' num2str(cod_av(i),'%0.3f')]},'FontSize',8)
            %set(gca,'YDir','reverse')
            set(gca,'FontSize',10)
            xlim([0.0 ceil(max_x + 0.05*max_x)])
            grid on
        end
        for j = i+1:length(haa)
            set(haa(j),'Visible','off')
        end

    end

    %Clear the structure
    clear A 
    
    %Save the data

    if avg_data
        outmat = [EnsNo_av' datevec(MTdatenum_av) latlon_av dmg Heading_av' Pitch_av' Roll_av' Temp_av' Depth_av' BeamDepths_av DABack_av' DAVeast_av' DAVnorth_av' DAVmag_av' DAVdir_av' DAVvert_av' ustar_av' zo_av' cod_av'];
    else
        outmat = [EnsNo MTyear MTmon MTday MThour MTmin MTsec latlon dmg Heading Pitch Roll Temp Depth BeamDepths DABack DAVeast DAVnorth DAVmag DAVdir DAVvert ustar zo cod]; 
    end
    
    % Replace nans with -9999 (ARCGIS takes nans to be zero when exporting to
    % shapefile)
    if 0  % Fill the nans  
        for i = 7:size(outmat,2)
            outmat(:,i) = inpaint_nans(outmat(:,i));
        end 
    else  %fill with -9999
        outmat(isnan(outmat)) = -9999;
    end
    
    
    
   
    if append_data & zi == 1
        %outfile = [fileName(1:end-4) '_GIS.csv'];
        firstfile = outfile;
    elseif ~append_data
        [ofile,opath] = uiputfile('*.csv','Save file name',ascii2gispath);
        outfile = [opath ofile];
        %outfile = [fileName(1:end-4) '_GIS.csv'];
    elseif append_data & zi > 1
        outfile = firstfile;
    end    
            
    
    
    if append_data & zi == 1
        ofid   = fopen(outfile, 'wt');
        outcount = fprintf(ofid,'EnsNo, Year, Month, Day, Hour, Min, Sec, Lat_WGS84, Lon_WGS84, DMG_m, Heading_deg,  Pitch_deg,  Roll_deg, Temp_C, Depth_m, B1Depth_m, B2Depth_m, B3Depth_m, B4Depth_m, BackScatter_db, DAVeast_cmps, DAVnorth_cmps, DAVmag_cmps, DAVdir_deg, DAVvert_cmps, U_star_mps, Z0_m, COD\n');
    elseif ~append_data
        ofid   = fopen(outfile, 'wt');
        outcount = fprintf(ofid,'EnsNo, Year, Month, Day, Hour, Min, Sec, Lat_WGS84, Lon_WGS84, DMG_m, Heading_deg,  Pitch_deg,  Roll_deg, Temp_C, Depth_m, B1Depth_m, B2Depth_m, B3Depth_m, B4Depth_m, BackScatter_db, DAVeast_cmps, DAVnorth_cmps, DAVmag_cmps, DAVdir_deg, DAVvert_cmps, U_star_mps, Z0_m, COD\n');
    elseif append_data & zi > 1
        ofid   = fopen(outfile, 'at');
    end
    outcount = fprintf(ofid,'%6.0f, %4.0f, %2.0f, %2.0f, %2.0f, %2.0f, %2.2f, %3.10f, %3.10f, %6.1f, %3.3f, %3.3f, %3.3f, %3.1f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.0f, %3.2f, %3.2f, %3.2f, %3.1f, %3.2f, %2.4f, %2.4f, %1.4f\n',outmat');
    fclose(ofid);
    
    if avg_data
        [Easting,Northing,utmzone] = deg2utm(latlon_av(:,1),latlon_av(:,2));
        VelOut = [VelOut; Easting Northing zeros(size(Easting)) (DAVeast_av)'./100 (DAVnorth_av)'./100];
    else
        [Easting,Northing,utmzone] = deg2utm(latlon(:,1),latlon(:,2));
        VelOut = [VelOut; Easting Northing zeros(size(Easting)) DAVeast./100 DAVnorth./100];
    end
    
%     if avg_data
%         ProfOut = [ProfOut; Easting Northing zeros(size(Easting)) (DAVeast_av)'./100 (DAVnorth_av)'./100];
%     else
%         ProfOut = [ProfOut; Easting Northing zeros(size(Easting)) DAVeast./100 DAVnorth./100];
%     end
       
    clear outmat latlon EnsNo MTyear MTmon MTday MThour MTmin MTsec latlon Heading Pitch Roll Temp Depth BeamDepths DABack DAVeast DAVnorth DAVmag DAVdir DAVvert ustar zo cod Easting Northing utmzone
    clear EnsNo_av MTdatenum_av latlon_av Heading_av Pitch_av Roll_av Temp_av Depth_av BeamDepths_av DABack_av DAVeast_av DAVnorth_av DAVmag_av DAVdir_av DAVvert_av ustar_av zo_av cod_av 
    
    waitbar(zi/z); %update the waitbar
end
delete(wbh) %remove the waitbar

msgbox('Conversion Complete','Conversion Status','help','replace');

% Save an *.anv file (for iRIC model interface)
goodrows = [];
for i = 1:size(VelOut,1)
    rowsum = sum(isnan(VelOut(i,:)));
    if rowsum == 0
        goodrows = [goodrows; i];
    end
end
%outfile = [fileName(1:end-4) '_DAV.anv'];
outfile = [outfile(1:end-4) '.anv'];
ofid    = fopen(outfile, 'wt');
outcount = fprintf(ofid,'%8.2f  %8.2f  %5.2f  %3.3f  %3.3f\n',VelOut(goodrows,:)');
fclose(ofid);

% Save an *.pro file (profile text file)
% goodrows = [];
% for i = 1:size(VelOut,1)
%     rowsum = sum(isnan(ProfOut(i,:)));
%     if rowsum == 0
%         goodrows = [goodrows; i];
%     end
% end
% %outfile = [fileName(1:end-4) '_DAV.anv'];
% outfile = [outfile(1:end-4) '.pro'];
% ofid    = fopen(outfile, 'wt');
% outcount = fprintf(ofid,'%8.2f  %8.2f  %5.2f  %3.3f  %3.3f\n',VelOut(goodrows,:)');
% fclose(ofid);


% Save the paths
% if exist('LastDir.mat') == 2
    % save('LastDir.mat','ascii2gispath','-append')
% else
    % save('LastDir.mat','ascii2gispath')
% end
else % User cancelled uigetfile dialog, return 
    VelOut = [];
    goodrows = [];
end

    

