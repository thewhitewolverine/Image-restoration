function varargout = Assignment_4(varargin)
% ASSIGNMENT_4 MATLAB code for Assignment_4.fig
%      ASSIGNMENT_4, by itself, creates a new ASSIGNMENT_4 or raises the existing
%      singleton*.
%
%      H = ASSIGNMENT_4 returns the handle to a new ASSIGNMENT_4 or the handle to
%      the existing singleton*.
%
%      ASSIGNMENT_4('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ASSIGNMENT_4.M with the given input arguments.
%
%      ASSIGNMENT_4('Property','Value',...) creates a new ASSIGNMENT_4 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Assignment_4_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Assignment_4_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Assignment_4

% Last Modified by GUIDE v2.5 17-Oct-2018 23:44:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Assignment_4_OpeningFcn, ...
                   'gui_OutputFcn',  @Assignment_4_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Assignment_4 is made visible.
function Assignment_4_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Assignment_4 (see VARARGIN)

% Choose default command line output for Assignment_4
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Assignment_4 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Assignment_4_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Upload_Blur_Kernel.
function Upload_Blur_Kernel_Callback(hObject, eventdata, handles)
% hObject    handle to Upload_Blur_Kernel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im imblur acb
[path, user_cance]=imgetfile();                            
if user_cance
    magbox(sprintf('Error'),'Error','Error');
    return
end
imblur=rgb2gray(im2double(imread(path)));  
[m,n]=size(im(:,:,1));
fb=fft2(imblur,m,n);
acb=fftshift(fb);
fshb=log(1+abs(acb));
axes(handles.axes2);
imshow(imblur,[]);
axes(handles.axes4);
imshow(fshb,[]);

% --- Executes on button press in Load_Image.
function Load_Image_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im im2 iml imlen im_fr im_fb im_fg M N             
iml={};                                                     
[path, user_cance]=imgetfile();                            
if user_cance
    magbox(sprintf('Error'),'Error','Error');
    return
end
im=imread(path);                                            
im=im2double(im);
im2=im;
iml{1}=im2;
imlen=1;
axes(handles.axes1);
imshow(im);
f=fft2(rgb2gray(im));
fsh=fftshift(f);

im_r=im(:,:,1);
im_b=im(:,:,2);
im_g=im(:,:,3);
im_fr=fftshift(fft2(im_r));
im_fg=fftshift(fft2(im_g));
im_fb=fftshift(fft2(im_b));
[M,N]=size(im_r);

ac=log(1+abs(fsh));
axes(handles.axes3);
imshow(ac,[]);


% --- Executes on button press in Full_inverse_filter.
function Full_inverse_filter_Callback(hObject, eventdata, handles)
% hObject    handle to Full_inverse_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im im_fr im_fb im_fg acb M N
disp(size(im_fr));
disp(size(acb));

im_fr1 = im_fr./((abs(acb)<4)*(4) + acb);
im_fb1 = im_fb./((abs(acb)<4)*(4) + acb);
im_fg1 = im_fg./((abs(acb)<4)*(4) + acb);

imnew=im;
im_fr2=ifft2(ifftshift(im_fr1));
im_fb2=ifft2(ifftshift(im_fb1));
im_fg2=ifft2(ifftshift(im_fg1));
imnew(:,:,1)=im_fr2.*50;
imnew(:,:,2)=im_fb2.*50;
imnew(:,:,3)=im_fg2.*50;
axes(handles.axes1);
imshow(imnew);

MSE = sum(sum((imnew(:,:,1)-im(:,:,1)).^2))/(M*N);
PSNR = 10*log10(256*256/MSE);   

set(handles.PSNR,'String',PSNR);

window = fspecial('gaussian', 11, 1.5);
K = [0.05 0.05];
L = 255;
img1 = double(imnew(:,:,1));
img2 = double(im(:,:,1));

f = max(1,round(min(M,N)/256));

if(f>1)
    lpf = ones(f,f);
    lpf = lpf/sum(lpf(:));
    img1 = imfilter(img1,lpf,'symmetric','same');
    img2 = imfilter(img2,lpf,'symmetric','same');

    img1 = img1(1:f:end,1:f:end);
    img2 = img2(1:f:end,1:f:end);
end

C1 = (K(1)*L)^2;
C2 = (K(2)*L)^2;
window = window/sum(sum(window));

mu1   = filter2(window, img1, 'valid');
mu2   = filter2(window, img2, 'valid');
mu1_sq = mu1.*mu1;
mu2_sq = mu2.*mu2;
mu1_mu2 = mu1.*mu2;
sigma1_sq = filter2(window, img1.*img1, 'valid') - mu1_sq;
sigma2_sq = filter2(window, img2.*img2, 'valid') - mu2_sq;
sigma12 = filter2(window, img1.*img2, 'valid') - mu1_mu2;

if (C1 > 0 & C2 > 0)
   ssim_map = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))./((mu1_sq + mu2_sq + C1).*(sigma1_sq + sigma2_sq + C2));
else
   numerator1 = 2*mu1_mu2 + C1;
   numerator2 = 2*sigma12 + C2;
    denominator1 = mu1_sq + mu2_sq + C1;
   denominator2 = sigma1_sq + sigma2_sq + C2;
   ssim_map = ones(size(mu1));
   index = (denominator1.*denominator2 > 0);
   ssim_map(index) = (numerator1(index).*numerator2(index))./(denominator1(index).*denominator2(index));
   index = (denominator1 ~= 0) & (denominator2 == 0);
   ssim_map(index) = numerator1(index)./denominator1(index);
end

mssim = mean2(ssim_map);

set(handles.SSIM,'String',mssim);



% --- Executes on button press in Truncated_filter.
function Truncated_filter_Callback(hObject, eventdata, handles)
% hObject    handle to Truncated_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im im_fr im_fb im_fg M N acb
c= str2double(get(handles.Radius, 'String'));
per=(c/100)*(M/2);
im_fr1 = im_fr./((abs(acb)<4)*(4) + acb);
im_fb1 = im_fb./((abs(acb)<4)*(4) + acb);
im_fg1 = im_fg./((abs(acb)<4)*(4) + acb);


x = -floor(M/2):floor(M/2);
y = -floor(N/2):floor(N/2);
[xx yy] = meshgrid(x,y);
u = zeros(size(xx));
u((xx.^2+yy.^2)<per^2)=1;
u=u(1:M,1:N);
im_fr2=u.*im_fr1;
im_fb2=im_fb1.*u;
im_fg2=im_fg1.*u;
imnew=im;
axes(handles.axes3);
imshow(log(1+abs(im_fr2)),[]); 

im_fr2=ifft2(ifftshift(im_fr2));
im_fb2=ifft2(ifftshift(im_fb2));
im_fg2=ifft2(ifftshift(im_fg2));
imnew(:,:,1)=im_fr2.*50;
imnew(:,:,2)=im_fb2.*50;
imnew(:,:,3)=im_fg2.*50;
axes(handles.axes1);
imshow(imnew,[]); 

MSE = sum(sum((imnew(:,:,1)-im(:,:,1)).^2))/(M*N);
PSNR = 10*log10(256*256/MSE);   

set(handles.PSNR,'String',PSNR);

window = fspecial('gaussian', 11, 1.5);
K = [0.05 0.05];
L = 255;
img1 = double(imnew(:,:,1));
img2 = double(im(:,:,1));

f = max(1,round(min(M,N)/256));

if(f>1)
    lpf = ones(f,f);
    lpf = lpf/sum(lpf(:));
    img1 = imfilter(img1,lpf,'symmetric','same');
    img2 = imfilter(img2,lpf,'symmetric','same');

    img1 = img1(1:f:end,1:f:end);
    img2 = img2(1:f:end,1:f:end);
end

C1 = (K(1)*L)^2;
C2 = (K(2)*L)^2;
window = window/sum(sum(window));

mu1   = filter2(window, img1, 'valid');
mu2   = filter2(window, img2, 'valid');
mu1_sq = mu1.*mu1;
mu2_sq = mu2.*mu2;
mu1_mu2 = mu1.*mu2;
sigma1_sq = filter2(window, img1.*img1, 'valid') - mu1_sq;
sigma2_sq = filter2(window, img2.*img2, 'valid') - mu2_sq;
sigma12 = filter2(window, img1.*img2, 'valid') - mu1_mu2;

if (C1 > 0 & C2 > 0)
   ssim_map = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))./((mu1_sq + mu2_sq + C1).*(sigma1_sq + sigma2_sq + C2));
else
   numerator1 = 2*mu1_mu2 + C1;
   numerator2 = 2*sigma12 + C2;
    denominator1 = mu1_sq + mu2_sq + C1;
   denominator2 = sigma1_sq + sigma2_sq + C2;
   ssim_map = ones(size(mu1));
   index = (denominator1.*denominator2 > 0);
   ssim_map(index) = (numerator1(index).*numerator2(index))./(denominator1(index).*denominator2(index));
   index = (denominator1 ~= 0) & (denominator2 == 0);
   ssim_map(index) = numerator1(index)./denominator1(index);
end

mssim = mean2(ssim_map);

set(handles.SSIM,'String',mssim);



function Radius_Callback(hObject, eventdata, handles)
% hObject    handle to Radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   
% Hints: get(hObject,'String') returns contents of Radius as text
% str2double(get(hObject,'String')) returns contents of Radius as a double



% --- Executes during object creation, after setting all properties.
function Radius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Weiner.
function Weiner_Callback(hObject, eventdata, handles)
% hObject    handle to Weiner (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im im_fr im_fb im_fg acb M N
k= str2double(get(handles.k_value, 'String'));
const=(abs((abs(acb)<4)*(4) + acb).^2)/((abs((abs(acb)<4)*(4) + acb).^2)+k);

im_fr1 = im_fr./((abs(acb)<4)*(4) + acb);
im_fb1 = im_fb./((abs(acb)<4)*(4) + acb);
im_fg1 = im_fg./((abs(acb)<4)*(4) + acb);

im_fr1=(im_fr1.*const);
im_fb1=(im_fb1.*const);
im_fg1=(im_fg1.*const);

imnew=im;
im_fr2=ifft2(ifftshift(im_fr1));
im_fb2=ifft2(ifftshift(im_fb1));
im_fg2=ifft2(ifftshift(im_fg1));
imnew(:,:,1)=im_fr2.*50;
imnew(:,:,2)=im_fb2.*50;
imnew(:,:,3)=im_fg2.*50;
axes(handles.axes1);
imshow(imnew,[]); 

MSE = sum(sum((imnew(:,:,1)-im(:,:,1)).^2))/(M*N);
PSNR = 10*log10(256*256/MSE);   

set(handles.PSNR,'String',PSNR);
window = fspecial('gaussian', 11, 1.5);
K = [0.05 0.05];
L = 255;
img1 = double(imnew(:,:,1));
img2 = double(im(:,:,1));

f = max(1,round(min(M,N)/256));

if(f>1)
    lpf = ones(f,f);
    lpf = lpf/sum(lpf(:));
    img1 = imfilter(img1,lpf,'symmetric','same');
    img2 = imfilter(img2,lpf,'symmetric','same');

    img1 = img1(1:f:end,1:f:end);
    img2 = img2(1:f:end,1:f:end);
end

C1 = (K(1)*L)^2;
C2 = (K(2)*L)^2;
window = window/sum(sum(window));

mu1   = filter2(window, img1, 'valid');
mu2   = filter2(window, img2, 'valid');
mu1_sq = mu1.*mu1;
mu2_sq = mu2.*mu2;
mu1_mu2 = mu1.*mu2;
sigma1_sq = filter2(window, img1.*img1, 'valid') - mu1_sq;
sigma2_sq = filter2(window, img2.*img2, 'valid') - mu2_sq;
sigma12 = filter2(window, img1.*img2, 'valid') - mu1_mu2;

if (C1 > 0 & C2 > 0)
   ssim_map = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))./((mu1_sq + mu2_sq + C1).*(sigma1_sq + sigma2_sq + C2));
else
   numerator1 = 2*mu1_mu2 + C1;
   numerator2 = 2*sigma12 + C2;
    denominator1 = mu1_sq + mu2_sq + C1;
   denominator2 = sigma1_sq + sigma2_sq + C2;
   ssim_map = ones(size(mu1));
   index = (denominator1.*denominator2 > 0);
   ssim_map(index) = (numerator1(index).*numerator2(index))./(denominator1(index).*denominator2(index));
   index = (denominator1 ~= 0) & (denominator2 == 0);
   ssim_map(index) = numerator1(index)./denominator1(index);
end

mssim = mean2(ssim_map);

set(handles.SSIM,'String',mssim);





function k_value_Callback(hObject, eventdata, handles)
% hObject    handle to k_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k_value as text
%        str2double(get(hObject,'String')) returns contents of k_value as a double


% --- Executes during object creation, after setting all properties.
function k_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cls_filter.
function cls_filter_Callback(hObject, eventdata, handles)
% hObject    handle to cls_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im im_fr im_fb im_fg acb M N
gamma= str2double(get(handles.gamma, 'String'));
p=[0 -1 0; -1 4 -1; 0 -1 0];
Pp=fft2(p,M,N);
const=((abs((abs(acb)<4)*(4) + acb).^2)+gamma*(abs(Pp).^2));

im_fr1=(conj(((abs(acb)<4)*(4) + acb))./const).*im_fr;
im_fb1=(conj(((abs(acb)<4)*(4) + acb))./const).*im_fb;
im_fg1=(conj(((abs(acb)<4)*(4) + acb))./const).*im_fg;

imnew=im;
im_fr2=ifft2(ifftshift(im_fr1));
im_fb2=ifft2(ifftshift(im_fb1));
im_fg2=ifft2(ifftshift(im_fg1));
imnew(:,:,1)=im_fr2.*50;
imnew(:,:,2)=im_fb2.*50;
imnew(:,:,3)=im_fg2.*50;
axes(handles.axes1);
imshow(imnew,[]); 

MSE = sum(sum((imnew(:,:,1)-im(:,:,1)).^2))/(M*N);
PSNR = 10*log10(256*256/MSE);   

set(handles.PSNR,'String',PSNR);

window = fspecial('gaussian', 11, 1.5);
K = [0.05 0.05];
L = 255;
img1 = double(imnew(:,:,1));
img2 = double(im(:,:,1));

f = max(1,round(min(M,N)/256));

if(f>1)
    lpf = ones(f,f);
    lpf = lpf/sum(lpf(:));
    img1 = imfilter(img1,lpf,'symmetric','same');
    img2 = imfilter(img2,lpf,'symmetric','same');

    img1 = img1(1:f:end,1:f:end);
    img2 = img2(1:f:end,1:f:end);
end

C1 = (K(1)*L)^2;
C2 = (K(2)*L)^2;
window = window/sum(sum(window));

mu1   = filter2(window, img1, 'valid');
mu2   = filter2(window, img2, 'valid');
mu1_sq = mu1.*mu1;
mu2_sq = mu2.*mu2;
mu1_mu2 = mu1.*mu2;
sigma1_sq = filter2(window, img1.*img1, 'valid') - mu1_sq;
sigma2_sq = filter2(window, img2.*img2, 'valid') - mu2_sq;
sigma12 = filter2(window, img1.*img2, 'valid') - mu1_mu2;

if (C1 > 0 & C2 > 0)
   ssim_map = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))./((mu1_sq + mu2_sq + C1).*(sigma1_sq + sigma2_sq + C2));
else
   numerator1 = 2*mu1_mu2 + C1;
   numerator2 = 2*sigma12 + C2;
    denominator1 = mu1_sq + mu2_sq + C1;
   denominator2 = sigma1_sq + sigma2_sq + C2;
   ssim_map = ones(size(mu1));
   index = (denominator1.*denominator2 > 0);
   ssim_map(index) = (numerator1(index).*numerator2(index))./(denominator1(index).*denominator2(index));
   index = (denominator1 ~= 0) & (denominator2 == 0);
   ssim_map(index) = numerator1(index)./denominator1(index);
end

mssim = mean2(ssim_map);

set(handles.SSIM,'String',mssim);

function gamma_Callback(hObject, eventdata, handles)
% hObject    handle to gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gamma as text
%        str2double(get(hObject,'String')) returns contents of gamma as a double


% --- Executes during object creation, after setting all properties.
function gamma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PSNR_Callback(hObject, eventdata, handles)
% hObject    handle to PSNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PSNR as text
%        str2double(get(hObject,'String')) returns contents of PSNR as a double
 

% --- Executes during object creation, after setting all properties.
function PSNR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PSNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SSIM_Callback(hObject, eventdata, handles)
% hObject    handle to SSIM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SSIM as text
%        str2double(get(hObject,'String')) returns contents of SSIM as a double


% --- Executes during object creation, after setting all properties.
function SSIM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SSIM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
