% Comparison TOPAZ - CTD within the Nordic Seas
clear all

yr  = input('Aar:');
mon = input('Mnd:');
mnd = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];
dyp = [5:5:2000];   % Interpolerer resultater til 5-m bins
dd1 = find(dyp > 6 & dyp < 101);    % 10-100 dybdeintervall
dd2 = find(dyp > 101 & dyp < 301);  % 100-300 dybdeintervall
dd3 = find(dyp > 301 & dyp < 801);  % 300-800 dybdeintervall
dd4 = find(dyp > 801 & dyp < 2001); % 800-2000 dybdeintervall
delta_s = 0.03;   % Kriterium for MLD

fpath = strcat('/home/vidarli/Prosjekt/Aktive/Copernicus/ARC_MFC/Monthly_Reports/Figures/', ...
    num2str(yr),'/M',num2str(mon,'%2.2d'),'/');

%% *** Laster inn modelldata ***
mpath = strcat('/home/vidarli/Prosjekt/Aktive/Copernicus/ARC_MFC/Data/CTD_Valid/', ...
    num2str(yr),'/M',num2str(mon),'/');
mfile = strcat(mpath,'TOPAZ_NRT_',num2str(yr),num2str(mon,'%2.2d'),'.nc');
%mfilet = strcat(mpath,'TOPAZ_NRT_T_',num2str(yr),num2str(mon,'%2.2d'),'.nc');
%mfiles = strcat(mpath,'TOPAZ_NRT_S_',num2str(yr),num2str(mon,'%2.2d'),'.nc');

mLon = ncread(mfile,'longitude',[650 200],[420 520]);
mLat = ncread(mfile,'latitude',[650 200],[420 520]);
mdep = ncread(mfile,'depth');
mtid = ncread(mfile,'time') / 86400 + datenum(1970,1,1) - .5;
mTmp = single(ncread(mfile,'thetao',[650 200 1 1],[420 520 Inf Inf]));
mSal = single(ncread(mfile,'so',[650 200 1 1],[420 520 Inf Inf]));
latM = reshape(mLat,[],1);
lonM = reshape(mLon,[],1);
lonlatM = [lonM, latM];

% Tar bort evnt. duplikater i tid
[ut,ta,tc] = unique(mtid);
mtid = mtid(ta);
mTmp = mTmp(:,:,:,ta);
mSal = mSal(:,:,:,ta);

[Xm,Ym,Dm,Tm] = size(mTmp);
if yr == 2017 && mon == 11   % Åpenbar feil i felt/dag 19 i Nov 17'
    mTmp(:,:,:,19) = mTmp(:,:,:,18);
    mSal(:,:,:,19) = mSal(:,:,:,18);
end

%% *** Laster inn CTD-data ***
opath = strcat('/home/vidarli/Prosjekt/Aktive/Copernicus/ARC_MFC/Data/CTD_Valid/', ...
    num2str(yr),'/M',num2str(mon),'/');
flist = strcat('/home/vidarli/Prosjekt/Aktive/Copernicus/ARC_MFC/Data/CTD_Valid/', ...
    num2str(yr),'/M',num2str(mon),'/filelist.txt');
fid   = fopen(flist);
F     = textscan(fid,'%s');
fclose(fid);
[f1,f2] = size(F{1});

% Initialiserer matrise
To = 0; Do = 0;
for f = 1:f1
    ofile = strcat(F{1}(f));
    tmp   = length(ncread(ofile{1},'TIME'));
    To    = To + tmp;
%    dmp   = size(ncread(ofile{1},'DEPH'),1);
    dmp   = size(ncread(ofile{1},'PRES'),1);
    if dmp > Do
        Do = dmp;
    end
end
oMat{1}(1:To)       = zeros;
oMat{2}(1:To)       = zeros;
oMat{3}(1:To)       = zeros;
oMat{4}(1:To)       = zeros;
oMat{5}(1:To)       = zeros;
oMat{6}(1:To,1:Do)  = zeros;
oMat{7}(1:To,1:Do)  = zeros;
oMat{8}(1:To,1:Do)  = zeros;
oMat{9}(1:To,1:Do)  = zeros;
oMat{10}(1:To,1:Do) = zeros;
oMat{11}(1:To,1:Do) = zeros;

teller = 1;
for f = 1:f1
    ofile = strcat(F{1}(f));
    tmp   = length(ncread(ofile{1},'TIME'));
    try
        dmp   = size(ncread(ofile{1},'DEPH'),1);
    catch exception
        dmp   = size(ncread(ofile{1},'PRES'),1);
    end
    for t = 1:tmp
        otid                   = ncread(ofile{1},'TIME',[t],[1]) + datenum(1950,1,1);
        [YY,MM,DD]             = datevec(otid);
        oMat{1}(teller)        = datenum(YY,MM,DD);
        oMat{2}(teller)        = ncread(ofile{1},'TIME_QC',[t],[1]);
        oMat{3}(teller)        = ncread(ofile{1},'LONGITUDE',[t],[1]);
        oMat{4}(teller)        = ncread(ofile{1},'LATITUDE',[t],[1]);
        oMat{5}(teller)        = ncread(ofile{1},'POSITION_QC',[t],[1]);
        % Sjekker om DEPH el. PRES:
        try
            oMat{6}(teller,1:dmp)  = ncread(ofile{1},'DEPH',[1 t],[Inf 1]);
            oMat{7}(teller,1:dmp)  = ncread(ofile{1},'DEPH_QC',[1 t],[Inf 1]);
        catch exception
            oMat{6}(teller,1:dmp)  = ncread(ofile{1},'PRES',[1 t],[Inf 1]);
            oMat{7}(teller,1:dmp)  = ncread(ofile{1},'PRES_QC',[1 t],[Inf 1]);
        end
        oMat{8}(teller,1:dmp)  = ncread(ofile{1},'TEMP',[1 t],[Inf 1]);
        oMat{9}(teller,1:dmp)  = ncread(ofile{1},'TEMP_QC',[1 t],[Inf 1]);
        % Sjekker først om CNDC er tilgj., og hvis ikke PSAL:
        try
            oMat_tmp = ncread(ofile{1},'CNDC',[1 t],[Inf 1]);
            oMat{11}(teller,1:dmp) = ncread(ofile{1},'CNDC_QC',[1 t],[Inf 1]);
            % Sjekker om S/m eller mS/cm og beregner saltholdighet fra conductivity:
            J = length(oMat_tmp(1:dmp));
            for j = 1:J
                if oMat_tmp(j) < 0
                    oMat_tmp(j) = NaN;
                elseif oMat_tmp(j) > 0 && oMat_tmp(j) < 10
                    oMat_tmp(j) = oMat_tmp(j);   % S/m
                elseif oMat_tmp(j) > 10 && oMat_tmp(j) < 50
                    oMat_tmp(j) = oMat_tmp(j) / 10;   % mS/cm -> S/m
                elseif oMat_tmp(j) > 50
                    oMat_tmp(j) = NaN;
                end
                oMat{10}(teller,j) = salinity(oMat_tmp(j),oMat{8}(teller,j),oMat{6}(teller,j));
            end
        catch exception
            oMat{10}(teller,1:dmp) = ncread(ofile{1},'PSAL',[1 t],[Inf 1]);
            oMat{11}(teller,1:dmp) = ncread(ofile{1},'PSAL_QC',[1 t],[Inf 1]);
        end
        teller = teller + 1;
    end
end

oMat{1}(oMat{1}==0)  = NaN;
oMat{2}(oMat{2}==0)  = NaN;
oMat{3}(oMat{3}==0)  = NaN;
oMat{4}(oMat{4}==0)  = NaN;
oMat{5}(oMat{5}==0)  = NaN;
oMat{6}(oMat{6}==0)  = NaN;
oMat{7}(oMat{7}==0)  = NaN;
oMat{8}(oMat{8}==0)  = NaN;
oMat{9}(oMat{9}==0)  = NaN;
oMat{10}(oMat{10}==0) = NaN;
oMat{11}(oMat{11}==0) = NaN;

% Tar bort dårlige data (QC = 4):
oMat{8}(oMat{9}==4)   = NaN;
oMat{10}(oMat{11}==4) = NaN;

% Tar manuelt bort urealistiske data som ikke er tatt bort av QC:
% Temperatur (-2 < T < 40):
err_to = find(oMat{8} > 40);
err_tu = find(oMat{8} < -2);
disp('Ant. åpenbart feile temperaturobservasjoner tatt bort:');
disp(length(err_to) + length(err_tu));
oMat{8}(oMat{8} > 40)   = NaN;
oMat{8}(oMat{8} < -2)   = NaN;

% Saltholdighet (0 < S < 50)
err_co = find(oMat{10} > 50);
err_cu = find(oMat{10} < 0);
disp('Ant. åpenbart feile saltholdighetsobservasjoner tatt bort:');
disp(length(err_co) + length(err_cu));


%% *** Sammenligner ***

Dc = length(dyp);
cTmp{1}(1:To,1:Dc)  = NaN;   % Obs temp på felles dyp
cTmp{2}(1:To,1:Dc)  = NaN;   % Mod temp på felles dyp
cSal{1}(1:To,1:Dc)  = NaN;   % Obs salt på felles dyp
cSal{2}(1:To,1:Dc)  = NaN;   % Mod salt på felles dyp
cDns{1}(1:To,1:Dc)  = NaN;   % Obs dens på felles dyp
cDns{2}(1:To,1:Dc)  = NaN;   % Mod dens på felles dyp
cmTmp{1}(1:To,1:Dm) = NaN;   % Mod temp på modelldyp
cmSal{1}(1:To,1:Dm) = NaN;   % Mod salt på modelldyp

teller = 1;
for t = 1:Tm
    tt = find(oMat{1}==mtid(t));
    if ~isempty(tt)
        for i = 1:length(tt)
            if oMat{2}(tt(i)) < 1 || oMat{2}(tt(i)) > 2 || oMat{5}(tt(i)) < 1 || oMat{5}(tt(i)) > 2 || oMat{3}(tt(i)) < -29 || oMat{3}(tt(i)) > 49 || oMat{4}(tt(i)) < 51 || oMat{4}(tt(i)) > 81
                continue
            else
                lonlatO = [oMat{3}(tt(i)), oMat{4}(tt(i))];
                pNh = dsearchn(lonlatM,lonlatO);
                [Im,Jm] = find(mLon==lonM(pNh) & mLat==latM(pNh));

                ii                  = find(~isnan(oMat{6}(tt(i),:)));
                if ~isempty(ii)
                    lon(teller)         = oMat{3}(tt(i));
                    lat(teller)         = oMat{4}(tt(i));
                    cTmp{1}(teller,:)   = interp1(oMat{6}(tt(i),ii),oMat{8}(tt(i),ii),dyp);
                    cSal{1}(teller,:)   = interp1(oMat{6}(tt(i),ii),oMat{10}(tt(i),ii),dyp);
                    cDns{1}(teller,:)   = sw_dens0(cSal{1}(teller,:),cTmp{1}(teller,:));
                    cTmp{2}(teller,:)   = interp1(mdep,squeeze(mTmp(Im,Jm,:,t)),dyp);
                    cSal{2}(teller,:)   = interp1(mdep,squeeze(mSal(Im,Jm,:,t)),dyp);
                    cDns{2}(teller,:)   = sw_dens0(cSal{2}(teller,:),cTmp{2}(teller,:));
                
                    cmTmp{1}(teller,:)  = mTmp(Im,Jm,:,t);
                    cmSal{1}(teller,:)  = mSal(Im,Jm,:,t);
                %NB! QualFlag
                
                % Beregner bias og rmse innenfor dybdeintervall
                    T_bias(1,teller) = nanmean(cTmp{2}(teller,dd1)) - nanmean(cTmp{1}(teller,dd1));
                    S_bias(1,teller) = nanmean(cSal{2}(teller,dd1)) - nanmean(cSal{1}(teller,dd1));
                    T_bias(2,teller) = nanmean(cTmp{2}(teller,dd2)) - nanmean(cTmp{1}(teller,dd2));
                    S_bias(2,teller) = nanmean(cSal{2}(teller,dd2)) - nanmean(cSal{1}(teller,dd2));
                    T_bias(3,teller) = nanmean(cTmp{2}(teller,dd3)) - nanmean(cTmp{1}(teller,dd3));
                    S_bias(3,teller) = nanmean(cSal{2}(teller,dd3)) - nanmean(cSal{1}(teller,dd3));
                    T_bias(4,teller) = nanmean(cTmp{2}(teller,dd4)) - nanmean(cTmp{1}(teller,dd4));
                    S_bias(4,teller) = nanmean(cSal{2}(teller,dd4)) - nanmean(cSal{1}(teller,dd4));
                % Bias totalt i vannsøylen 5-2000m:
                    oo = find(~isnan(cTmp{1}(teller,:))); mm = find(~isnan(cTmp{2}(teller,:)));
                    cc = intersect(oo,mm);
                    T_biast(teller)  = nanmean(cTmp{2}(teller,cc)) - nanmean(cTmp{1}(teller,cc));
                    S_biast(teller)  = nanmean(cSal{2}(teller,cc)) - nanmean(cSal{1}(teller,cc));
                
                % Beregner MLD bias
                    ds1 = find(cDns{1}(teller,:) - cDns{1}(teller,1) > delta_s);
                    ds2 = find(cDns{2}(teller,:) - cDns{2}(teller,1) > delta_s);
                    if isempty(ds1) && isempty(ds2) || isempty(cc)
                        MLD_bias(teller) = 0;
                    elseif isempty(ds1) && ~isempty(ds2)
                        MLD_bias(teller) = dyp(min(ds2)) - dyp(max(cc));
                    elseif ~isempty(ds1) && isempty(ds2)
                        MLD_bias(teller) = dyp(max(cc)) - dyp(min(ds1));
                    else
                        MLD_bias(teller) = dyp(min(ds2)) - dyp(min(ds1));
                    end
                else
                    continue
                end   % if 'isempty(ii)'
                teller = teller + 1;
            end   % if 'oMat'
        end   % for 'i'
    end   % if '~isempty(tt)'
end   % for 't'

T_rmse(1) = sqrt(nanmean(T_bias(1,:) .* T_bias(1,:)));
T_rmse(2) = sqrt(nanmean(T_bias(2,:) .* T_bias(2,:)));
T_rmse(3) = sqrt(nanmean(T_bias(3,:) .* T_bias(3,:)));
T_rmse(4) = sqrt(nanmean(T_bias(4,:) .* T_bias(4,:)));
S_rmse(1) = sqrt(nanmean(S_bias(1,:) .* S_bias(1,:)));
S_rmse(2) = sqrt(nanmean(S_bias(2,:) .* S_bias(2,:)));
S_rmse(3) = sqrt(nanmean(S_bias(3,:) .* S_bias(3,:)));
S_rmse(4) = sqrt(nanmean(S_bias(4,:) .* S_bias(4,:)));

%% *** Plotter ***
 
load /home/vidarli/Data/Nordic_4km/Nordic4km.dat;
Lonkyst = Nordic4km(:,1); Latkyst = Nordic4km(:,2);
ytick = ['  10-100 m';' 100-300 m';' 300-800 m';'800-2000 m'];
redblue_colormap2;

% Kart over bias
figure;   % Temperatur
fillseg([Lonkyst, Latkyst], [.7 .7 .7], [.7 .7 .7]); hold on;
scatter(lon,lat,[],T_biast,'filled'); colormap(rb_map); caxis([-2 2]);
cb=colorbar('peer',gca); set(get(cb,'title'),'string','\circC','fontsize',12);
title(strcat('Temperature bias',{' '},mnd(mon,:),num2str(yr)),'fontsize',12);
xlabel('Longitude [\circE]','fontsize',12);
ylabel('Latitude [\circN]','fontsize',12);
set(gca,'xlim',[-25 50],'ylim',[54 82],'fontsize',12);
set(gcf,'color','w','Renderer','painters');
figfile = strcat(fpath,'TOPAZ_CTD_Comparison_MapTemp_',num2str(yr),num2str(mon,'%2.2d'));
eval(['export_fig ' figfile, ' -eps -transparent']);

figure;   % Saltholdighet
fillseg([Lonkyst, Latkyst], [.7 .7 .7], [.7 .7 .7]); hold on;
scatter(lon,lat,[],S_biast,'filled'); colormap(rb_map); caxis([-1 1]); colorbar;
title(strcat('Salinity bias',{' '},mnd(mon,:),num2str(yr)),'fontsize',12);
xlabel('Longitude [\circE]','fontsize',12);
ylabel('Latitude [\circN]','fontsize',12);
set(gca,'xlim',[-25 50],'ylim',[54 82],'fontsize',12);
set(gcf,'color','w','Renderer','painters');
figfile = strcat(fpath,'TOPAZ_CTD_Comparison_MapSalt_',num2str(yr),num2str(mon,'%2.2d'));
eval(['export_fig ' figfile, ' -eps -transparent']);

figure;   % MLD
fillseg([Lonkyst, Latkyst], [.7 .7 .7], [.7 .7 .7]); hold on;
scatter(lon,lat,[],MLD_bias,'filled'); colormap(rb_map); caxis([-250 250]);
cb=colorbar('peer',gca); set(get(cb,'title'),'string','m','fontsize',12);
title(strcat('MLD bias',{' '},mnd(mon,:),num2str(yr)),'fontsize',12);
xlabel('Longitude [\circE]','fontsize',12);
ylabel('Latitude [\circN]','fontsize',12);
set(gca,'xlim',[-25 50],'ylim',[54 82],'fontsize',12);
set(gcf,'color','w','Renderer','painters');
figfile = strcat(fpath,'TOPAZ_CTD_Comparison_MapMLD_',num2str(yr),num2str(mon,'%2.2d'));
eval(['export_fig ' figfile, ' -eps -transparent']);


% Bias vs. Dyp:
% Temperatur:
fig = figure;
set(fig,'position',[5 5 400 600]);

% Lager legend først:
plot(nanmean(T_bias(1,:)),1,'k*','linewidth',1.5); hold on;
plot(T_rmse(1),1,'r*','linewidth',1.5);
legend('Bias','RMSE','orientation','horizontal');
% Bias:
plot(T_bias(1,:),1,'*','color',[.7 .7 .7]);
plot(nanmean(T_bias(1,:)),1,'k*','linewidth',1.5);
plot(T_bias(2,:),2,'*','color',[.7 .7 .7]);
plot(nanmean(T_bias(2,:)),2,'k*','linewidth',1.5);
plot(T_bias(3,:),3,'*','color',[.7 .7 .7]);
plot(nanmean(T_bias(3,:)),3,'k*','linewidth',1.5);
plot(T_bias(4,:),4,'*','color',[.7 .7 .7]);
plot(nanmean(T_bias(4,:)),4,'k*','linewidth',1.5);
% RMSE:
plot(T_rmse(1),1,'r*','linewidth',1.5);
plot(T_rmse(2),2,'r*','linewidth',1.5);
plot(T_rmse(3),3,'r*','linewidth',1.5);
plot(T_rmse(4),4,'r*','linewidth',1.5);
% Skriver inn verdier:
text(-1.5,1.15,num2str(nanmean(T_bias(1,:)),'%2.2f'),'color','k','fontsize',12);
text(1,1.15,num2str(T_rmse(1),'%2.2f'),'color','r','fontsize',12);
text(-1.5,2.15,num2str(nanmean(T_bias(2,:)),'%2.2f'),'color','k','fontsize',12);
text(1,2.15,num2str(T_rmse(2),'%2.2f'),'color','r','fontsize',12);
text(-1.5,3.15,num2str(nanmean(T_bias(3,:)),'%2.2f'),'color','k','fontsize',12);
text(1,3.15,num2str(T_rmse(3),'%2.2f'),'color','r','fontsize',12);
text(-1.5,4.15,num2str(nanmean(T_bias(4,:)),'%2.2f'),'color','k','fontsize',12);
text(1,4.15,num2str(T_rmse(4),'%2.2f'),'color','r','fontsize',12);

title(strcat('TOPAZ vs. CTD',{' '},mnd(mon,:),num2str(yr)),'fontsize',12);
xlabel('Temperature difference [deg C]','fontsize',12);
%ylabel('Depth range','fontsize',12);
set(gca,'xlim',[-2 2],'xtick',[-2:1:2],'ylim',[.5 4.5],'ydir','reverse', ...
    'ytick',[1 2 3 4],'yticklabel',ytick,'fontsize',12);

set(gcf,'paperpositionmode','auto');
figfile = strcat(fpath,'TOPAZ_CTD_Comparison_Temp_',num2str(yr),num2str(mon,'%2.2d'));
eval(['export_fig ' figfile, ' -eps -transparent']);

% Saltholdighet:
fig = figure;
set(fig,'position',[5 5 400 600]);

% Lager legend først:
plot(nanmean(S_bias(1,:)),1,'k*','linewidth',1.5); hold on;
plot(S_rmse(1),1,'r*','linewidth',1.5);
legend('Bias','RMSE','orientation','horizontal');
% Bias:
plot(S_bias(1,:),1,'*','color',[.7 .7 .7]);
plot(nanmean(S_bias(1,:)),1,'k*','linewidth',1.5);
plot(S_bias(2,:),2,'*','color',[.7 .7 .7]);
plot(nanmean(S_bias(2,:)),2,'k*','linewidth',1.5);
plot(S_bias(3,:),3,'*','color',[.7 .7 .7]);
plot(nanmean(S_bias(3,:)),3,'k*','linewidth',1.5);
plot(S_bias(4,:),4,'*','color',[.7 .7 .7]);
plot(nanmean(S_bias(4,:)),4,'k*','linewidth',1.5);
% RMSE:
plot(S_rmse(1),1,'r*','linewidth',1.5);
plot(S_rmse(2),2,'r*','linewidth',1.5);
plot(S_rmse(3),3,'r*','linewidth',1.5);
plot(S_rmse(4),4,'r*','linewidth',1.5);
% Skriver inn verdier:
text(-.25,1.15,num2str(nanmean(S_bias(1,:)),'%2.2f'),'color','k','fontsize',12);
text(.25,1.15,num2str(S_rmse(1),'%2.2f'),'color','r','fontsize',12);
text(-.25,2.15,num2str(nanmean(S_bias(2,:)),'%2.2f'),'color','k','fontsize',12);
text(.25,2.15,num2str(S_rmse(2),'%2.2f'),'color','r','fontsize',12);
text(-.25,3.15,num2str(nanmean(S_bias(3,:)),'%2.2f'),'color','k','fontsize',12);
text(.25,3.15,num2str(S_rmse(3),'%2.2f'),'color','r','fontsize',12);
text(-.25,4.15,num2str(nanmean(S_bias(4,:)),'%2.2f'),'color','k','fontsize',12);
text(.25,4.15,num2str(S_rmse(4),'%2.2f'),'color','r','fontsize',12);

title(strcat('TOPAZ vs. CTD',{' '},mnd(mon,:),num2str(yr)),'fontsize',12);
xlabel('Salinity difference','fontsize',12);
%ylabel('Depth range','fontsize',12);
set(gca,'xlim',[-1 1],'xtick',[-1:.5:1],'ylim',[.5 4.5],'ydir','reverse', ...
    'ytick',[1 2 3 4],'yticklabel',ytick,'fontsize',12);

set(gcf,'paperpositionmode','auto');
figfile = strcat(fpath,'TOPAZ_CTD_Comparison_Salt_',num2str(yr),num2str(mon,'%2.2d'));
eval(['export_fig ' figfile, ' -eps -transparent']);


% Prøver T-S diagram
ss=(32:.1:36);
tt=(-2:.1:15);
[TT,SS]=meshgrid(tt,ss);
DD=sw_dens0(SS,TT); DD=DD-1000;

figure;
% Lager legend først:
plot(cSal{1}(1,1),cTmp{1}(1,1),'.','color',[.7 .7 .7]); hold on;
plot(cmSal{1}(1,1),cmTmp{1}(1,1),'r.');
legend('Obs','TOPAZ','location','southwest');
plot(cSal{1},cTmp{1},'.','color',[.7 .7 .7]);
plot(cmSal{1},cmTmp{1},'r.');
[cs,ct]=contour(SS,TT,DD,[20:.5:29],'k:'); clabel(cs,ct,'fontsize',8);

title(strcat('TOPAZ vs. CTD',{' '},mnd(mon,:),num2str(yr)),'fontsize',12);
xlabel('Salinity','fontsize',12);
ylabel('Temperature [deg C]','fontsize',12);
set(gca,'xlim',[33 35.5],'ylim',[-2 12],'fontsize',12);

figfile = strcat(fpath,'TOPAZ_CTD_Comparison_TS_',num2str(yr),num2str(mon,'%2.2d'));
eval(['export_fig ' figfile, ' -eps -transparent']);



