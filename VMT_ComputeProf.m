function [dnm,vem,vnm,vvm,vmm,maxd,vdir] = VMT_ComputeProf(dn,ve,vn,vv,d)
% Computes the median profile from multiemsemble data.
%
% Inputs:
%
% dn: matrix of surface origin bin locations (depths) (#bins x #ens)
% ve:  matrix of observed east velocity (#bins x #ens)
% vn:  matrix of observed north velocity (#bins x #ens)
% vv:  matrix of observed vertical velocity (#bins x #ens)
% d:  depth of water column for each ensemble
%
% P.R. Jackson, USGS, 9-12-14

n = size(dn,1);

% Compute the max depth
maxd = nanmax(d);

% Aggregate the data for later plotting
Ve_all = [];
Vn_all = [];
Vv_all = [];
dn_all = [];
for i = 1:length(d)
    Ve_all = [Ve_all; ve(:,i)];
    Vn_all = [Vn_all; vn(:,i)];
    Vv_all = [Vv_all; vv(:,i)];
    dn_all = [dn_all; dn(:,i)];
end

% Compute median values for each bin
vem = nanmedian(ve,2);
vnm = nanmedian(vn,2);
vvm = nanmedian(vv,2);
dnm = dn(:,1);
obs = nan*ones(n,1);  %preallocate
for i = 1:n
    obs(i) = sum(~isnan(ve(i,:)));
end

% Find any cells that have < 20% of the median number of data points all
% other bins
indx = find(obs < 0.2*nanmedian(obs));

%Remove lean cells
vem(indx) = nan;
vnm(indx) = nan;
vvm(indx) = nan;

%Compute the Velocity Magnitude
vmm = sqrt(vem.^2 + vnm.^2);

%Compute the direction
vdir = 90 - (atan2(vnm, vem))*180/pi; %Compute the atan from the velocity components, convert to radians, and rotate to north axis
qindx = find(vdir < 0);
if ~isempty(qindx)
    vdir(qindx) = vdir(qindx) + 360;  %Must add 360 deg to Quadrant 4 values as they are negative angles from the +y axis
end

% Plot the data

if 0  %for debugging
    figure(1); clf
    subplot(1,4,1)
    plot(Ve_all,dn_all,'k.'); hold on
    plot(vem,dnm,'ro-')
    %ylim([0 1])
    ylabel('Depth')
    xlabel('East Velocity')
    grid on
    ylim([0 maxd])
    set(gca,'YDir','reverse')
    subplot(1,4,2)
    plot(Vn_all,dn_all,'k.'); hold on
    plot(vnm,dnm,'ro-')
    %ylim([0 1])
    xlabel('North Velocity')
    grid on
    ylim([0 maxd])
    set(gca,'YDir','reverse')
    subplot(1,4,3)
    plot(Vv_all,dn_all,'k.'); hold on
    plot(vvm,dnm,'ro-')
   % ylim([0 1])
    xlabel('Vertical Velocity')
    grid on
    ylim([0 maxd])
    set(gca,'YDir','reverse')
    subplot(1,4,4)
    plot(vmm,dnm,'ro-')
   % ylim([0 1])
    xlabel('Velocity Magnitude')
    grid on
    ylim([0 maxd])
    set(gca,'YDir','reverse')
   
end



    


