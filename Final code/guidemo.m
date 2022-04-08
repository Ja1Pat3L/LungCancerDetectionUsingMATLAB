function varargout = guidemo(varargin)
% GUIDEMO MATLAB code for guidemo.fig
%      GUIDEMO, by itself, creates a new GUIDEMO or raises the existing
%      singleton*.
%
%      H = GUIDEMO returns the handle to a new GUIDEMO or the handle to
%      the existing singleton*.
%
%      GUIDEMO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIDEMO.M with the given input arguments.
%
%      GUIDEMO('Property','Value',...) creates a new GUIDEMO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guidemo_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guidemo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guidemo

% Last Modified by GUIDE v2.5 04-Jun-2018 11:52:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guidemo_OpeningFcn, ...
                   'gui_OutputFcn',  @guidemo_OutputFcn, ...
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


% --- Executes just before guidemo is made visible.
function guidemo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guidemo (see VARARGIN)

% Choose default command line output for guidemo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guidemo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = guidemo_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file path] = uigetfile('*.bmp;*.jpg','Pick an Image File');
file=strcat(path,file);
inp = imread(file);
inp = imresize(inp,[256,256]);

if size(inp,3)>1
  
  inp1 = rgb2gray(inp);
  
end
axes(handles.axes1);
imshow(inp1);
title('Input image');
handles.inp1=inp1;
handles.inp=inp;
guidata(hObject,handles);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

inp1=handles.inp1;
inp=handles.inp;
t0 = 60;Nregs = 2;
th = t0+((max(inp1(:))+min(inp1(:)))./2);

for i=1:1:size(inp1,1)
    for j=1:1:size(inp1,2)
       if inp1(i,j)>th
          sout(i,j) = 1; 
       else
          sout(i,j) = 0; 
       end    
    
    end
end

imshow(sout);
axes(handles.axes1);
title('Segmented image');


Sout = imfill(sout,'holes'); 

for i=1:1:size(inp,1)
    for j=1:1:size(inp,2)
       if (Sout(i,j)~=0)
          Lout(i,j) = inp(i,j); 
       else
          Lout(i,j) = 0; 
       end    
    
    end
end

axes(handles.axes2);
imshow(Lout,[]);
title('Backgroung removal');



%%%%%Extraction of Desired region from an input image

Iout = imfill(sout,'holes');      %%% Filling inner region
[Iout cnt] = bwlabel(Iout,8);
Iout = bwareaopen(Iout,5500);
[Iout cnt] = bwlabel(Iout,8);
prop = regionprops(Iout,'BoundingBox');
lprop = round(prop.BoundingBox);
Lrr = imsubtract(Iout,sout);
Lreg = imfill(Iout,'holes'); 
Lrr = imfill(Lrr,'holes');
axes(handles.axes3);
imshow(Lrr);
title('Regional mask');
[rr cc] = size(Lrr);
Lp = zeros(rr,cc);
for ii=1:1:size(Lrr,1)
    for jj=1:1:size(Lrr,2)
   
        if isequal(Lrr(ii,jj),1)
          Lp(ii,jj) = inp(ii,jj);
        else
          Lp(ii,jj) = 0; 
        end   
        
    end
end
Lp = uint8(Lp);
axes(handles.axes4);
imshow(Lp);
handles.Lout=Lout;
handles.Nregs=Nregs;

guidata(hObject,handles);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

inp1=handles.inp1;
cd dtwtdecomp

[Faf,Fsf] = FSfarras; % first stage filters
[af,sf] = dualfilt1;  % second stage filters
dtcomp = cplxdual2D(double(inp1),1,Faf,af);

cd .. 

tde1 = dtcomp{:,1}; tde2 = dtcomp{:,2}; 
tde11 = tde1{:,1}; tde12 = tde1{:,2};
tde21 = tde2{:,1}; tde22 = tde2{:,2};
tde211 = tde21{:,1}; tde212 = tde21{:,2}; 
tde221 = tde22{:,1}; tde222 = tde22{:,2};

tde111 = tde11{:,1}; tde112 = tde11{:,2};
tde121 = tde12{:,1}; tde122 = tde12{:,2};

tde111a = tde111{:,1}; tde111b = tde111{:,2}; tde111c = tde111{:,3};
tde112a = tde112{:,1}; tde112b = tde112{:,2}; tde112c = tde112{:,3};

tde121a = tde121{:,1}; tde121b = tde121{:,2}; tde121c = tde121{:,3};
tde122a = tde122{:,1}; tde122b = tde122{:,2}; tde122c = tde122{:,3};

%%% -----------> Feature Extraction for Training Smaple <--------------%%%

 gm1 = graycomatrix(tde111a);
gprops1 = graycoprops(gm1);
F1(1) = gprops1.Contrast; F2(1) = gprops1.Correlation;
F3(1) = gprops1.Energy; F4(1) = gprops1.Homogeneity;

gm2 = graycomatrix(tde111b);
gprops2 = graycoprops(gm2);
F1(2) = gprops2.Contrast; F2(2) = gprops2.Correlation;
F3(2) = gprops2.Energy; F4(2) = gprops2.Homogeneity;


gm3 = graycomatrix(tde111c);
gprops3 = graycoprops(gm3);
F1(3) = gprops3.Contrast; F2(3) = gprops3.Correlation;
F3(3) = gprops3.Energy; F4(3) = gprops3.Homogeneity;

gm4 = graycomatrix(tde112a);
gprops4 = graycoprops(gm4);
F1(4) = gprops4.Contrast; F2(4) = gprops4.Correlation;
F3(4) = gprops4.Energy; F4(4) = gprops4.Homogeneity;

gm5 = graycomatrix(tde112b);
gprops5 = graycoprops(gm5);
F1(5) = gprops5.Contrast;  F2(5) = gprops5.Correlation;
F3(5) = gprops5.Energy; F4(5) = gprops5.Homogeneity;

gm6 = graycomatrix(tde112c);
gprops6 = graycoprops(gm6);
F1(6) = gprops6.Contrast;  F2(6) = gprops6.Correlation;
F3(6) = gprops6.Energy; F4(6) = gprops6.Homogeneity;


gm7 = graycomatrix(tde121a);
gprops7 = graycoprops(gm7);
F1(7) = gprops7.Contrast;  F2(7) = gprops7.Correlation;
F3(7) = gprops7.Energy;  F4(7) = gprops7.Homogeneity;


gm8 = graycomatrix(tde121b);
gprops8 = graycoprops(gm8); 
F1(8) = gprops8.Contrast;  F2(8) = gprops8.Correlation;
F3(8) = gprops8.Energy; F4(8) = gprops8.Homogeneity;

pause(0.06);

Qfeat = [mean(F1) mean(F2) mean(F3) mean(F4)]';
disp(Qfeat);

helpdlg('Training is completed');

%%%%%%%%Importing the network Properties %%%%    
handles.Qfeat = Qfeat;

guidata(hObject,handles);



% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Nregs=handles.Nregs;
Lout=handles.Lout;
inp=handles.inp;
Qfeat = handles.Qfeat;
nmask = isnan(Qfeat);
Qfeat(find(nmask==1)) = 0;
netp = nnlearn;
load netp;
result = sim(netp,Qfeat); 

result= round(mean2(result));
 
tic
if result==1
        
   msgbox('Normal','Decision: ');
        
elseif result==2
        
   msgbox('Benign','Decision: ');
 
%%%%%Lloyds Clustering for Cells Segmentation
[AA1, AA2, AA3, AA4] = Lclustering(Lout);

figure('Name','Segmented Results','MenuBar','None');
subplot(2,2,1);imshow(AA1,[]);
title('Cluster: 1');
subplot(2,2,2);imshow(AA2,[]);
title('Cluster: 2');
subplot(2,2,3);imshow(AA3,[]);
title('Cluster: 3');
subplot(2,2,4);imshow(AA4,[]);
title('Cluster: 4');

%%%%%%%Postprocessing to extract the parenchymas
C1Bw = im2bw(AA3);
C2Bw  = im2bw(AA4);

Cout = bitor(C1Bw,C2Bw);
Lfout = imfill(Cout,'holes');

Lparen = bitxor(Lfout,Cout);
L1seg = Lparen;

Lparen = imfill(Lparen,'holes');
Lparen = bwareaopen(Lparen,80);
L2seg = Lparen;

figure('Name','Lung Parenchymas Mask','MenuBar','None')
imshow(Lparen,[])

%%%%%Morphological Filtering Process to locate abnormal lung region
S1e = strel('line',5,0);
S2e = strel('line',5,90);

Lparen = imdilate(Lparen,[S1e,S2e]);

[Lcout,Lcnt] = bwlabel(Lparen,8);
Lprops = regionprops(Lcout,'Area','BoundingBox');

for ii=1:1:Lcnt
    Lar(ii) = Lprops(ii).Area;
end
[Lmax,Lind]= sort(Lar,'descend');

Lcors = bwboundaries(Lcout);



for i=1:1:Nregs
    Lxycor{i} = Lcors{Lind(i)}; 
    dim = Lprops(Lind(i)).BoundingBox;  
    L1g = imcrop(L1seg,[dim(1),dim(2),dim(3),dim(4)]);
    L2g = imcrop(L2seg,[dim(1),dim(2),dim(3),dim(4)]);
    Lg = bitxor(L1g,L2g);
    Lcor = find(Lg~=0);
    Acnt(i) = length(Lcor);
    sLout{i} = imresize(Lg,[150,75]);
end   

figure('Name','Segmented Lesions from Lungs','MenuBar','None');
subplot(1,2,1); imshow(sLout{1});
subplot(1,2,2); imshow(sLout{2});

figure('Name','Tracing Lung Parenchymas','MenuBar','None')
imshow(inp);
hold on;
for i=1:1:length(Lxycor)
    bound = Lxycor{i};
    plot(bound(:,2),bound(:,1),'c','LineWidth',3);
end
hold off;
pause(0.08);

figure('Name','Tracing Abnormal Lung with Red','MenuBar','None')
imshow(inp);
hold on;
for i=1:1:length(Lxycor)
    if Lmax(i)>1920
       bound = Lxycor{i};
       plot(bound(:,2),bound(:,1),'g','LineWidth',2);
    else
       bound = Lxycor{i};
       plot(bound(:,2),bound(:,1),'r','LineWidth',2);
    end  
    
    if Acnt(i)>85
       bound = Lxycor{i};
       plot(bound(:,2),bound(:,1),'r','LineWidth',2);
    end  
end
hold off;

elseif result==3
        
   msgbox('Malignant','Decision: ');
 
%%%%%Lloyds Clustering for Cells Segmentation
[AA1, AA2, AA3, AA4] = Lclustering(Lout);

figure('Name','Segmented Results','MenuBar','None');
subplot(2,2,1);imshow(AA1,[]);
title('Cluster: 1');
subplot(2,2,2);imshow(AA2,[]);
title('Cluster: 2');
subplot(2,2,3);imshow(AA3,[]);
title('Cluster: 3');
subplot(2,2,4);imshow(AA4,[]);
title('Cluster: 4');

%%%%%%%Postprocessing to extract the parenchymas
C1Bw = im2bw(AA3);
C2Bw  = im2bw(AA4);

Cout = bitor(C1Bw,C2Bw);
Lfout = imfill(Cout,'holes');

Lparen = bitxor(Lfout,Cout);
L1seg = Lparen;

Lparen = imfill(Lparen,'holes');
Lparen = bwareaopen(Lparen,80);
L2seg = Lparen;

figure('Name','Lung Parenchymas Mask','MenuBar','None')
imshow(Lparen,[])

%%%%%Morphological Filtering Process to locate abnormal lung region
S1e = strel('line',5,0);
S2e = strel('line',5,90);

Lparen = imdilate(Lparen,[S1e,S2e]);

[Lcout,Lcnt] = bwlabel(Lparen,8);
Lprops = regionprops(Lcout,'Area','BoundingBox');

for ii=1:1:Lcnt
    Lar(ii) = Lprops(ii).Area;
end
[Lmax,Lind]= sort(Lar,'descend');

Lcors = bwboundaries(Lcout);



for i=1:1:Nregs
    Lxycor{i} = Lcors{Lind(i)}; 
    dim = Lprops(Lind(i)).BoundingBox;  
    L1g = imcrop(L1seg,[dim(1),dim(2),dim(3),dim(4)]);
    L2g = imcrop(L2seg,[dim(1),dim(2),dim(3),dim(4)]);
    Lg = bitxor(L1g,L2g);
    Lcor = find(Lg~=0);
    Acnt(i) = length(Lcor);
    sLout{i} = imresize(Lg,[150,75]);
end   

figure('Name','Segmented Lesions from Lungs','MenuBar','None');
subplot(1,2,1); imshow(sLout{1});
subplot(1,2,2); imshow(sLout{2});

figure('Name','Tracing Lung Parenchymas','MenuBar','None')
imshow(inp);
hold on;
for i=1:1:length(Lxycor)
    bound = Lxycor{i};
    plot(bound(:,2),bound(:,1),'c','LineWidth',3);
end
hold off;
pause(0.08);

figure('Name','Tracing Abnormal Lung with Red','MenuBar','None')
imshow(inp);
hold on;
for i=1:1:length(Lxycor)
    if Lmax(i)>1920
       bound = Lxycor{i};
       plot(bound(:,2),bound(:,1),'g','LineWidth',2);
    else
       bound = Lxycor{i};
       plot(bound(:,2),bound(:,1),'r','LineWidth',2);
    end  
    
    if Acnt(i)>85
       bound = Lxycor{i};
       plot(bound(:,2),bound(:,1),'r','LineWidth',2);
    end  
end
hold off;



else
   helpdlg('Database Updation Required');

end
exe_time = toc;

disp('Processing Time(in Sec) :');
disp(exe_time);



%%%%Parameters Evaluation %%%%%%total number of test samples 9
   Tp = 8; Fn = 0;  %%%%%%%after classification
   Fp = 1; Tn = 1;  %%%%%Tp --> Abnormality correctly classified as abnormal
                    %%%%%Fn --> Abnormality incorrectly classified as normal
                    %%%%%Fp --> Normal incorrectly classified as abnormal
                    %%%%%Tn --> Normal correctly classified as normal
                      
Sensitivity = (Tp./(Tp+Fn)).*100;
Specificity = (Tn./(Tn+Fp)).*100;

Accuracy = ((Tp+Tn)./(Tp+Tn+Fp+Fn)).*100;


disp('Sensitivity: '); disp(Sensitivity);
disp('Specificity: '); disp(Specificity);
disp('Accuracy:'); disp(Accuracy);



figure('Name','Performance Metrics','MenuBar','none'); 
bar3(1,Sensitivity,0.3,'m');
hold on;
bar3(2,Specificity,0.3,'r');
hold on;
bar3(3,Accuracy,0.3,'g');
hold off;

title('Performance Metrics');
xlabel('Parametrics--->');
zlabel('Value--->');
 legend('Sensitivity','Specificity','Accuracy');
