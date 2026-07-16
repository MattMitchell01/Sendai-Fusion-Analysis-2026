function [CoordsOfAreaToExclude] = Choose_Area_To_Exclude(ImageWindow,ImageToShow,DefaultPathname)

%CoordsOfAreaToExclude comes out as [Xmin, Ymin, Width, Height] and 2nd
%dimension is the index of the area chosen.

[PatchFileName, PatchPathname] = uigetfile('*.tif','Select patch .tif file',...
    DefaultPathname);

ImageToShow = imread(strcat(PatchPathname,PatchFileName),1);
    MinToShow = 312;
    MaxToShow = 3000;

Choose_Another_Area = 'y';

UpdatedImage = ImageToShow;
CropWindow = figure(4);
set(0,'CurrentFigure',CropWindow);
imshow(UpdatedImage, [MinToShow, MaxToShow], 'InitialMagnification', 'fit');
drawnow;
Num_Areas = 0;

while Choose_Another_Area == 'y'
    Num_Areas = Num_Areas + 1;
    
    %Choose the region to be exlcuded.
    [~, CoordsOfAreaToExclude(1:4,Num_Areas)] = imcrop(CropWindow);
        XMin_ToRemove = round(CoordsOfAreaToExclude(1,Num_Areas));
        YMin_ToRemove = round(CoordsOfAreaToExclude(2,Num_Areas));
        Width_ToRemove = round(CoordsOfAreaToExclude(3,Num_Areas));
        Height_ToRemove = round(CoordsOfAreaToExclude(4,Num_Areas));
    
    %Black out the region just chosen
        UpdatedImage(YMin_ToRemove:YMin_ToRemove+Height_ToRemove,...
            XMin_ToRemove:XMin_ToRemove+Width_ToRemove) = 0;
        imshow(UpdatedImage, [MinToShow, MaxToShow], 'InitialMagnification', 'fit');
        drawnow;
    
    %Ask the user if they want to continue
        Prompts = {'Choose another area to exclude? (y/n)'};
        DefaultInputs = {'y',...
            };
        ChooseQuestion = inputdlg(Prompts,'Choose Another Area?', 1, DefaultInputs, 'on');

        Choose_Another_Area = char(ChooseQuestion(1,1));
    
end

end