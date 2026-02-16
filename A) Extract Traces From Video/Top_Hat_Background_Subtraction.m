function [Top_Hat_Filtered_Matrix] = Top_Hat_Background_Subtraction (Video_Matrix, Options)

Disk_Radius = Options.DiskRadius;

Disk_Element = strel('disk', Disk_Radius);   

Top_Hat_Filtered_Matrix = zeros( ...
    size(Video_Matrix,1), ...
    size(Video_Matrix,2), ...
    size(Video_Matrix,3), ...
    'uint16');

for i = 1:size(Video_Matrix, 3)
    
    Current_Frame_Matrix = Video_Matrix(:,:,i);
    
    Top_Hat_Filtered_Matrix(:,:,i) = imtophat(Current_Frame_Matrix, Disk_Element);

end