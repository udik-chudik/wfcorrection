
%clear all
%close all

fig = openfig('C:\Users\user0\Documents\MatlabScripts\wfcorrection\Tavrov\ga_42x42_649times_reqularization2_1000nm_RMS.fig');
%fig = openfig('C:\Users\user0\Documents\MatlabScripts\wfcorrection\Tavrov\ga_42x42_481times_reqularization2_500nm_RMS.fig');
%fig = openfig('C:\Users\user0\Documents\MatlabScripts\wfcorrection\Tavrov\ga_42x42_723times_reqularization2_100nm_RMS.fig');
%fig = openfig('C:\Users\user0\Documents\MatlabScripts\wfcorrection\Tavrov\ga_42x42_25times_reqularization2_40nm_RMS.fig');
% fig = openfig('C:\Users\user0\Documents\MatlabScripts\wfcorrection\Tavrov\ga_42x42_8times_reqularization2_10nm_RMS.fig');
h=figure(fig);
img=getimage(h);
%rectangle(h)
% rectangle(fig)
a=get(gca, 'Children')
a(1).Position;
a(2).Position;
%figure(33), imagesc(img), ...
%   R1 = rectangle('Position',[110 100 20 20]), ...   % original area to null
%   R2 = rectangle('Position',[366 100 20 20]),...   % nulled area
%    colorbar
figure(33), imagesc(img), ...
   R1 = rectangle('Position',a(2).Position), ...   % original area to null
   R2 = rectangle('Position',a(1).Position),...   % nulled area
    colorbar

% nulled region
FrameR2=img(R2.Position(2):(R2.Position(2)+R2.Position(4)),...
    R2.Position(1):(R2.Position(1)+R2.Position(4)));
figure(25), imagesc(FrameR2), colorbar
% original region
FrameR1=img(R1.Position(2):(R1.Position(2)+R1.Position(4)),...
    R1.Position(1):(R1.Position(1)+R1.Position(4)));
figure(26), imagesc(FrameR1), colorbar
figure(27), surf(FrameR2), hold on, surf(FrameR1), hold off

figure(28), hist(reshape([FrameR2;FrameR1], ...
    1, (size([FrameR2;FrameR1],1)*size([FrameR2;FrameR1],2))),21);

[counts, centers]=hist(reshape([FrameR2;FrameR1], ...
    1, (size([FrameR2;FrameR1],1)*size([FrameR2;FrameR1],2))),21);
MeanValR2=mean(mean(FrameR2)); MeanValR1=mean(mean(FrameR1)); 

StrToFig29={num2str(round(100*MeanValR2)/100), num2str(round(100*MeanValR1)/100)};
Str1ToFig29=num2str((MeanValR1-MeanValR2));
figure(29), hist(reshape([FrameR2;FrameR1], ...
    1, (size([FrameR2;FrameR1],1)*size([FrameR2;FrameR1],2))),21), ...
    hold on
    plot([MeanValR2,MeanValR2],[0,1.1*max(counts)],'-r',...
        [MeanValR1,MeanValR1],[0,1.1*max(counts)],'-r',...
        [MeanValR2,MeanValR1],[1.05*max(counts),1.05*max(counts)],'-r'), ...
    text([MeanValR2,MeanValR1],1.1.*[max(counts), max(counts)], StrToFig29),...
    text(mean([MeanValR1,MeanValR2]), max(counts), Str1ToFig29),  hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(28), surf(10.^FrameR2), hold on, surf(10.^FrameR1), hold off


% select 2x bigger region nulled
XY_FrameR2_2x = [R2.Position(1)-R2.Position(3)/2, R2.Position(1)+2*R2.Position(3), ...  % X
                R2.Position(2)-R2.Position(4)/2, R2.Position(2)+2*R2.Position(4)];      %Y
FrameR2_2x =img(XY_FrameR2_2x(3):XY_FrameR2_2x(4),XY_FrameR2_2x(1):XY_FrameR2_2x(2));
figure(35), imagesc(FrameR2_2x), colorbar
figure(36), mesh(10.^FrameR2_2x), colorbar


% select 2x bigger region original
XY_FrameR1_2x = [R1.Position(1)-R1.Position(3)/2, R1.Position(1)+2*R1.Position(3), ...  % X
                R1.Position(2)-R1.Position(4)/2, R1.Position(2)+2*R1.Position(4)];      %Y
FrameR1_2x =img(XY_FrameR1_2x(3):XY_FrameR1_2x(4),XY_FrameR1_2x(1):XY_FrameR1_2x(2));
figure(37), imagesc(FrameR1_2x), colorbar
figure(38), mesh(10.^FrameR1_2x)

% select 2x bigger region nulled 1x1 point
XY_FrameR2_2x = [R2.Position(1)-R2.Position(3)/2, R2.Position(1)+2*R2.Position(3), ...  % X
                R2.Position(2)-R2.Position(4)/2, R2.Position(2)+2*R2.Position(4)];      %Y
FrameR2_2x =img(XY_FrameR2_2x(3):XY_FrameR2_2x(4),XY_FrameR2_2x(1):XY_FrameR2_2x(2));
FrameR2_2x(1,1)=max(max(FrameR1_2x(1,1)));
figure(35), imagesc(FrameR2_2x), colorbar
figure(36), mesh(10.^FrameR2_2x), colorbar
figure(40), surf(FrameR2_2x), hold on, surf(FrameR1_2x), hold off
size(img)