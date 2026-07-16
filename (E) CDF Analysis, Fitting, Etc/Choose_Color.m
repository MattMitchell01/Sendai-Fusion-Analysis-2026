function [CurrentColor] = Choose_Color(ColorFlag,InputNumber)
    if strcmp(ColorFlag,'Data') || strcmp(ColorFlag,'Single Fit')
        FileNumber = InputNumber;

            if FileNumber == 1
                CurrentColor.DataPoints = 'bo';
                CurrentColor.FitLine = 'b--';
                CurrentColor.FaceValue = [0 0.45 0.74];
            elseif FileNumber == 2
                CurrentColor.DataPoints = 'mo';
                CurrentColor.FitLine = 'm--';
                CurrentColor.FaceValue = [1 0 1];
            elseif FileNumber == 3
                CurrentColor.DataPoints = 'go';
                CurrentColor.FitLine = 'g--';
                CurrentColor.FaceValue = [0.11 0.89 0.11];
            elseif FileNumber == 4
                CurrentColor.DataPoints = 'ko';
                CurrentColor.FitLine = 'k--';
                CurrentColor.FaceValue = [0.93 0.69 0.13];
            elseif FileNumber == 5
                CurrentColor.DataPoints = 'co';
                CurrentColor.FitLine = 'c--';
                CurrentColor.FaceValue = [0 1 1];
            elseif FileNumber == 6
                CurrentColor.DataPoints = 'ro';
                CurrentColor.FitLine = 'r--';
                CurrentColor.FaceValue = [1 0 0];
            elseif FileNumber == 7
                CurrentColor.DataPoints = 'yo';
                CurrentColor.FitLine = 'y--';
                CurrentColor.FaceValue = [1 1 0];
            elseif FileNumber == 8
                CurrentColor.DataPoints = 'bx';
                CurrentColor.FitLine = 'b--';
                CurrentColor.FaceValue = [0 0 1];
            elseif FileNumber == 9
                CurrentColor.DataPoints = 'mx';
                CurrentColor.FitLine = 'm--';
                CurrentColor.FaceValue = [1 0 1];
            else
                CurrentColor.DataPoints = 'bo';
                CurrentColor.FitLine = 'b-';
                CurrentColor.FaceValue = [0 0 1];
            end
            CurrentColor.ResidualPoints = CurrentColor.DataPoints;

    elseif strcmp(ColorFlag,'MultipleFits') 

            FitNumber = InputNumber;
            
            if FitNumber == 1
        %             CurrentColor.DataPoints = 'bo';
                    CurrentColor.FitLine = 'k-';
                    CurrentColor.ResidualPoints = 'ko';
                elseif FitNumber == 2
        %             CurrentColor.DataPoints = 'mo';
                    CurrentColor.ResidualPoints = 'mo';
                    CurrentColor.FitLine = 'm-';
                elseif FitNumber == 3
        %             CurrentColor.DataPoints = 'go';
                    CurrentColor.FitLine = 'g-';
                    CurrentColor.ResidualPoints = 'go';
                elseif FitNumber == 4
        %             CurrentColor.DataPoints = 'ko';
                    CurrentColor.FitLine = 'c-';
                    CurrentColor.ResidualPoints = 'co';
            end
    else
        disp("!!!Error: Plotting colors not chosen correctly")
        StopProgramNow
    end

end