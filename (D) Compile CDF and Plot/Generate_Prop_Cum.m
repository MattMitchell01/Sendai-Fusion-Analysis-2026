function [CumX, CumY] = Generate_Prop_Cum(SortedpHtoFList)
%Changes SortedpHtoFList to proportion hemifused over time and makes cum
%dist(Not normalized)

    Y = 1:length(SortedpHtoFList);
    New_Pt = 0;
    
    for i = 1:length(Y)
       if i == length(Y)
            New_Pt = New_Pt +1;
           CumX(New_Pt) = SortedpHtoFList(i);
           CumY(New_Pt) = Y(i);
           continue
       end
       if SortedpHtoFList(i) == SortedpHtoFList(i+1)
           continue
       else
           New_Pt = New_Pt +1;
           CumX(New_Pt) = SortedpHtoFList(i);
           CumY(New_Pt) = Y(i);
       end
    end
 
end
