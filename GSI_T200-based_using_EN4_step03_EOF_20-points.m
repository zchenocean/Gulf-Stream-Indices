clear all
clc
close all
%% read data ====================================================
% load data .mat : taux, lon, lat, time
file1 = './NCfile/EN4_20points_T200_time_series_1954.01-2024.09.nc';
var = 'T_po';
T9 = nc_varget(file1,var);
po_lon = nc_varget(file1,'po_lon');
po_lat = nc_varget(file1,'po_lat');
%
T9 = T9(:,1:(2024-1954)*12+9);
time = 1:1:(2024-1954)*12+9;
% get rid of nan points
% there are some times for which ALL data are NaN
nans = find(all(isnan(T9),1)==1);
T9(:,nans) = [];
time(nans) = [];
% remove positions for which ANY data are NaN
land = find(any(isnan(T9),2)==1);
ocean = find(~any(isnan(T9),2)==1);
T9(land,:) = [];
% minus the mean of the data
T9 = T9 -  repmat(nanmean(T9,2),[1 length(time)]);
%for ii=1:9
%    T9(ii,:) = detrend(T9(ii,:));
%end
%covmethod=1; %for correlation matrix
covmethod=0; %for covariance matrix
% detrend
% detrend all temperatures to focus on decadal variability not long-term trends
%Tall=detrend(([T1 T2 T3 T4 T5 T6 T7 T8 T9]));
Tall=T9(:,:)'; %(([T1 T2 T3 T4 T5 T6 T7 T8 T9]));
Teast=T9(11:20,:)';%(([T5 T6 T7 T8 T9]));
Twest=T9(1:10,:)'; %(([T1 T2 T3 T4 T5]));
if covmethod==0
   cx=cov(Tall);
else
   cx=corrcoef(Tall);%this is to pre-normalize result 
end
[v,d]=eig(cx);
var=trace(d);
eigval=diag(d);
%sort into descending variance
[eis,is]=sort(eigval);
is=flipud(is);
a=Tall*v;
gs1=a(:,is(1))-mean(a(:,is(1)));
gs2=mean(Tall') ./std(Tall');
gs3=a(:,is(2))-mean(a(:,is(2)));
%estimate amplitude (in deg.C) of fist eof at each of thenine longitudes
tamp1=sqrt((var/9)*v(:,is(1)))';
%%%%
if covmethod==0
   cxe=cov(Teast);
else
   cxe=corrcoef(Teast);%this is to pre-normalize result 
end
[ve,de]=eig(cxe);
vare=trace(de);
eigvale=diag(de);
%sort into descending variance
[eise,ise]=sort(eigvale);
ise=flipud(ise);
ae=Teast*ve;
gse=ae(:,ise(1))-mean(ae(:,ise(1)));
if covmethod==0
    cxw=cov(Twest);
else
    cxw=corrcoef(Twest);%this is to pre-normalize result 
end
[vw,dw]=eig(cxw);
varw=trace(dw);
eigvalw=diag(dw);
%sort into descending variance
[eisw,isw]=sort(eigvalw);
isw=flipud(isw);
aw=Twest*vw;
gsw=aw(:,isw(1))-mean(aw(:,isw(1)));
%y=[35.67 36.25 37.5 38.5 39 39.67 40.6 40.67 40.6];
%y=[35.67 36.25 37.5 38.5 39 39.8 40.6 40.67 40.6];%this veector corrected at 60W
%y=[36.25 37.5 38.5 39 39.8 40.4 40.67 40.6];
%y=[35.75,36.25,36.95,37.4,37.5,37.8,38.2,38.5,38.8,39.0,39.5,39.8,40.2,40.3,40.4,40.5,40.67,40.64,40.62,40.6]
y = po_lat;
%insert y values for meanwall position
%x=[-73 -70 -67 -65 -63 -60 -58 -55];
%x = [-74,-73,-72,-71,-70,-69,-68,-67,-66,-65,-64,-63,-62,-61,-60,-59,-58,-57,-56,-55];
x = po_lon;
%xe=[-63 -60 -58 -55];
%xe=[-64,-63,-62,-61,-60,-59,-58,-57,-56,-55];
xe = po_lon(11:20);
%xw=[-73 -70 -67 -65];
%xw=[-74,-73,-72,-71,-70,-69,-68,-67,-66,-65];
xw = po_lat(1:10);
figure(1)
clf
subplot(121)
plot(100*eigval(is)/var,[1:length(x)],'x');
hold on
plot(100*eigval(is)/var,[1:length(x)]);
grid
ylabel('eof #')
xlabel('percent variance')
subplot(122)
h=plot(x,v(:,is(3)),'g--',x,v(:,is(2)),'r:',x,v(:,is(1)),'k','linewidth',1.5);
set(h(3),'linewidth',3)
ylabel('eof strucure')
xlabel('longitude')
legend('eof #3','eof #2','eof #1')
grid
figure(2)
clf
subplot(121)
plot(100*eigvale(ise)/vare,[1:length(xe)],'x')
hold on
plot(100*eigvale(ise)/vare,[1:length(xe)]);
grid
ylabel('eof #, east group')
xlabel('percent variance')
subplot(122)
h=plot(xe,ve(:,ise(3)),'g--',xe,ve(:,ise(2)),'r:',xe,ve(:,ise(1)),'k','linewidth',1.5);
set(h(3),'linewidth',3)
ylabel('eof strucure')
xlabel('longitude')
legend('eof #3','eof #2','eof #1')
grid
figure(3)
clf
subplot(121)
plot(100*eigvalw(isw)/varw,[1:length(xw)],'x')
hold on
plot(100*eigvalw(isw)/varw,[1:length(xw)]);
grid
ylabel('eof #, west group')
xlabel('percent variance')
subplot(122)
h=plot(xw,vw(:,isw(3)),'g--',xw,vw(:,isw(2)),'r:',xw,vw(:,isw(1)),'k','linewidth',1.5);
set(h(3),'linewidth',3)
ylabel('eof strucure')
xlabel('longitude')
legend('eof #3','eof #2','eof #1')
grid

gs=gs1; gst=time';

%this stuff was once used to store some intermediate results
xjunk=detrend(Tall);
Tall_detrend=xjunk;
gswi=gsw/std(gsw);
gsw_detrend=detrend(gsw);
gswi_detrend=gsw_detrend/std(gsw_detrend);
gsi=gs/std(gs);
gs_detrend=detrend(gs);
gsi_detrend=gs_detrend/std(gs_detrend);
save T200_GSI_EN4_1954.01-2024.09_monthly.mat Tall Tall_detrend gsi gsi_detrend gswi gswi_detrend gst x y 
% % remove the seasonal cycle
% Da = zeros(size(Data));
% w = 2*pi/365.25;
% Ea = [ones(size(time(:))) cos(w*time(:)) sin(w*time(:))]; 
% for m=1:size(Data,1)
%   d = Data(m,:);
%   coeff = Ea\d(:);
%   dfit = Ea*coeff;
%   Da(m,:) = dfit';
% end
% Data = Data-Da;


%% EOF

%length_lon = size(Data,1);% lon length
%length_lat = size(Data,2);% lat length

