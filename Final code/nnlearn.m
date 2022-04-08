function netp = nnlearn 

Nsamples = 9;
% ldr = waitbar(0,'Please wait....');
for di = 1 : Nsamples
    
    ifile = strcat(int2str(di),'.jpg');
cd Database        
    inp = imread(ifile);
cd .. 
   inp = imresize(inp,[256,256]);

if size(inp,3)>1
  
  inp = rgb2gray(inp);
  
end


%%%%%Applying threshold to suppress background
t0 = 60;
th = t0+((max(inp(:))+min(inp(:)))./2);

for i=1:1:size(inp,1)
    for j=1:1:size(inp,2)
       if inp(i,j)>th
          sout(i,j) = 1; 
       else
          sout(i,j) = 0; 
       end    
    
    end
end

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

%%% transformation %%%

cd dtwtdecomp

        [Faf,Fsf] = FSfarras; % first stage filters

        [af,sf] = dualfilt1;  % second stage filters

        dtcomp = cplxdual2D(double(inp),1,Faf,af);

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

        Dtemp = [mean(F1) mean(F2) mean(F3) mean(F4)]';
        dfeatures(:,di) = Dtemp;  

% waitbar(di/Nsamples,ldr);

end


disp('Reference Feature Set :');
disp(dfeatures);

save dfeatures dfeatures;

nmask = isnan(dfeatures);
dfeatures(find(nmask==1)) = 0;

%%%%%Neural network creation and training 
%%%%Assigning target to each class features
Nc = 3; T=1;
save dfeatures dfeatures;
Nfeatures = size(dfeatures,1);

for dfi=1:1:Nsamples
    if Nc<1
      T = T+1;
      Nc =2;
      acti(1:Nsamples,dfi) = T;
      
    else
      acti(1:Nsamples,dfi) = T;  
      Nc = Nc-1;  
    end
end

tv = acti;
pv = dfeatures;
       
%%%%%% Initialize the backpropagation with Feed forward network
netp = newff(pv,tv);
netp = train(netp,pv,tv);

save netp netp;

helpdlg('NN Training Completed')

return;

