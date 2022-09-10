function [imgHISPEC]=HisSpecifBP()
% Paper: G. Jiang, S.C.F. Lin, C.Y. Wong, M.A. Rahman, T.R. Ren, 
% Ngaiming Kwok, Haiyan Shi, Ying-Hao Yu, Tonghai Wu, "Color image
% enhancement with brightness preservation using a histogram specification 
% approach," Optik 126 (2015) 5656�5664.
clc; clear all; close all; dbstop if error;
img.h=360; img.w=480; img.L=255; img.L_=linspace(0,img.L,img.L+1);
job=0; fig=[];
[fn,pn,bn]=uigetfile('Images\\*.jpg');
if bn>0,
  img.JPG=imread([pn fn]);
  [img,img.RGB]=MyImResize(img,img.JPG);
  fig(end+1)=MyFigure(); MyImshow(fig(end),img.RGB,fn,job); job=job+1;
  [img.HISPEC]=MyHISPEC(img,img.RGB);
  fig(end+1)=MyFigure(); MyImshow(fig(end),img.HISPEC,fn,job); 
  imgHISPEC=img.HISPEC;
%   mnRGB=mean(img.RGB,3); mnRGB=mean(mnRGB(:)); disp(mnRGB);
%   mnHISPEC=mean(imgHISPEC,3); mnHISPEC=mean(mnHISPEC(:)); disp(mnHISPEC);
end;
function [ADJ]=MyHISPEC(img,RGB)
HSI=rgb2hsi(RGB); GRY=HSI(:,:,3);
[hx,ix]=hist(GRY(:),img.L_/img.L);
[STR]=MyStretch(img,RGB);
HSI=rgb2hsi(STR); GRY=HSI(:,:,3); mn_gry=mean(GRY(:))*img.L;
h=ones(1,img.L+1);
if mn_gry>img.L/2,
  if mn_gry>img.L*3/4,
    h=linspace(0,1,img.L+1);
  else
    A=floor(4*mn_gry-2*img.L);
    h(1:A)=linspace(0,1,A);
  end;
elseif mn_gry<img.L/2,
  if mn_gry<img.L/4,
    h=linspace(1,0,img.L+1);
  else
    B=ceil(4*mn_gry-img.L);
    h(B+1:end)=linspace(1,0,img.L-B+1);
  end;
end; h=h/sum(h);
HEQ=histeq(GRY,h); mn_heq=mean(HEQ(:));
HSI(:,:,3)=HEQ; ADJ=hsi2rgb(HSI); ADJ=uint8(ADJ*img.L);
function [RGB]=MyStretch(img,RGB)
RGB=double(RGB)/img.L;
for k=1:size(RGB,3),
  mi=min(min(RGB(:,:,k))); RGB(:,:,k)=RGB(:,:,k)-mi;
  mx=max(max(RGB(:,:,k))); RGB(:,:,k)=RGB(:,:,k)/mx;
end;
RGB=uint8(RGB*img.L);
function [fig]=MyFigure()
mon=get(0,'MonitorPositions');
x=(rand*0.1+0.1)*mon(3);
y=(rand*0.1+0.2)*mon(4);
fig=figure('units','pixel','position',[x y 600 400]);
function []=MyImshow(fig,jmg,fn,job)
global hue sat ent grd mbr svfig DirResult;
jmg=jmg(2:end,2:end,:);% remove boundary
figure(fig); imshow(uint8(jmg));
set(gca,'position',[0 0 1 1]);
set(gcf,'name',fn); drawnow;
function [img,jmg]=MyImResize(img,jmg)
[v,u,w]=size(jmg);
if u>v,% landscape
  k=[img.h img.w];
else% portrait
  k=[img.w img.h];
end;
jmg=imresize(jmg,k,'bilinear');% resize
[img.V,img.U,img.N]=size(jmg);% image size 
img.N=img.V*img.U;
function hsi = rgb2hsi(rgb)
%RGB2HSI Converts an RGB image to HSI.
% Extract the individual component images.
rgb = im2double(rgb);
r = rgb(:, :, 1);
g = rgb(:, :, 2);
b = rgb(:, :, 3);
% Implement the conversion equations.
num = 0.5*((r - g) + (r - b));
den = sqrt((r - g).^2 + (r - b).*(g - b));
theta = acos(num./(den + eps));
H = theta;
H(b > g) = 2*pi - H(b > g);
H = H/(2*pi);
num = min(min(r, g), b);
den = r + g + b;
den(den == 0) = eps;
S = 1 - 3.* num./den;
H(S == 0) = 0;
I = (r + g + b)/3;
% Combine all three results into an hsi image.
hsi = cat(3, H, S, I);
function rgb = hsi2rgb(hsi)
%HSI2RGB Converts an HSI image to RGB.
% Extract the individual HSI component images.
H = hsi(:, :, 1) * 2 * pi;
S = hsi(:, :, 2);
I = hsi(:, :, 3);
% Implement the conversion equations.
R = zeros(size(hsi, 1), size(hsi, 2));
G = zeros(size(hsi, 1), size(hsi, 2));
B = zeros(size(hsi, 1), size(hsi, 2));
% RG sector (0 <= H < 2*pi/3).
idx = find( (0 <= H) & (H < 2*pi/3));
B(idx) = I(idx) .* (1 - S(idx));
R(idx) = I(idx) .* (1 + S(idx) .* cos(H(idx)) ./ ...
                                          cos(pi/3 - H(idx)));
G(idx) = 3*I(idx) - (R(idx) + B(idx));
% BG sector (2*pi/3 <= H < 4*pi/3).
idx = find( (2*pi/3 <= H) & (H < 4*pi/3) );
R(idx) = I(idx) .* (1 - S(idx));
G(idx) = I(idx) .* (1 + S(idx) .* cos(H(idx) - 2*pi/3) ./ ...
                    cos(pi - H(idx)));
B(idx) = 3*I(idx) - (R(idx) + G(idx));
% BR sector.
idx = find( (4*pi/3 <= H) & (H <= 2*pi));
G(idx) = I(idx) .* (1 - S(idx));
B(idx) = I(idx) .* (1 + S(idx) .* cos(H(idx) - 4*pi/3) ./ ...
                                           cos(5*pi/3 - H(idx)));
R(idx) = 3*I(idx) - (G(idx) + B(idx));
% Combine all three results into an RGB image.  Clip to [0, 1] to
% compensate for floating-point arithmetic rounding effects.
rgb = cat(3, R, G, B);
rgb = max(min(rgb, 1), 0);
