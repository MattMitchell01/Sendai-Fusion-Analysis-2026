function Save_Figure(FigureHandle,SavePath)
set(0,'CurrentFigure',FigureHandle)
Filename = strcat(SavePath,'CDFForPeter');
print(Filename,'-deps')
end