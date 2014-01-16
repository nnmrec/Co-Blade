function disp_CF = dispCorrectionFactor(zFrac, NUM_SEC, I1, I2)
% This function calculates correction factors for displacment of tapered cantilever beams.  
% Code written and validated by Matthew Trudeau, Master's Student, PSU ARL
% Trudeau, M.G., "Structural and Hydrodynamic Design Optimization Enhancements
%                 with Application to Marine Hydrokinetic Turbine Blades," Master's
%                 Thesis, The Pennsylvania State University, 2011.

Iratio1 = zeros(NUM_SEC, 1);    
Iratio2 = zeros(NUM_SEC, 1);    
for n = 1:NUM_SEC
    Iratio1(n) = I1(1) / I1(n);      
    Iratio2(n) = I2(1) / I2(n);
end

%Average the two Iratios together
Iratio = (Iratio1 + Iratio2) ./ 2;

%% This entire loop interpolates from the correction values for deflection of tapered beams.  
%  W. Young and R. Budynas, Roark's Formulas for Stress and Strain. McGraw-Hill, 7th ed., 2001.
%  this uses n = 4 in Table 8.11(d) of the above book 

interpl_Iratio1      = zeros(NUM_SEC, 1);   
interpl_Iratio2      = zeros(NUM_SEC, 1);  
interpl_y_correction = zeros(NUM_SEC, 1); 
for i = 1:NUM_SEC  
    if Iratio(i) <= 0.5  
        if zFrac(i) > 0 && zFrac(i) < 0.25
            x1  = 0;
            x2  = 0.25;
            y1  = 3.864;
            y2  = 3.532;
            y11 = 1.976;
            y22 = 1.886;
            
            x = zFrac(i);
            interpl_Iratio1(i) = y1  +   (y2-y1)/(x2-x1) * (x-x1);
            interpl_Iratio2(i) = y11 + (y22-y11)/(x2-x1) * (x-x1);
             
        elseif zFrac(i) > 0.25 && zFrac(i) < 0.5
            x1  = 0.25;
            x2  = 0.5;
            y1  = 3.532;
            y2  = 3.2;
            y11 = 1.886;
            y22 = 1.796;
            
            x = zFrac(i);
            interpl_Iratio1(i) = y1  +   (y2-y1)/(x2-x1) * (x-x1);
            interpl_Iratio2(i) = y11 + (y22-y11)/(x2-x1) * (x-x1);
        
        elseif zFrac(i) > 0.5 && zFrac(i) < 0.75
            x1  = 0.5;
            x2  = 0.75;
            y1  = 3.2;
            y2  = 2.971;
            y11 = 1.796;
            y22 = 1.727;
            
            x = zFrac(i);
            interpl_Iratio1(i) = y1  +   (y2-y1)/(x2-x1) * (x-x1);
            interpl_Iratio2(i) = y11 + (y22-y11)/(x2-x1) * (x-x1);
       
        elseif zFrac(i) > 0.75 && zFrac(i) < 1
            x1  = 0.75;
            x2  = 1;
            y1  = 2.971;
            y2  = 2.828;
            y11 = 1.727;
            y22 = 1.682;
            
            x = zFrac(i);
            interpl_Iratio1(i) = y1  +   (y2-y1)/(x2-x1) * (x-x1);
            interpl_Iratio2(i) = y11 + (y22-y11)/(x2-x1) * (x-x1);
       
        elseif zFrac(i) == 0.25
            interpl_Iratio1(i) = 3.532;
            interpl_Iratio2(i) = 1.886;
        
        elseif zFrac(i) == 0.5    
            interpl_Iratio1(i) = 3.2;
            interpl_Iratio2(i) = 1.796;
     
        elseif zFrac(i) == 0.75
            interpl_Iratio1(i) = 2.971;
            interpl_Iratio2(i) = 1.727;
           
        elseif zFrac(i) == 1    
            interpl_Iratio1(i) = 2.828;
            interpl_Iratio2(i) = 1.682;
           
        end
        
        x1 = 0;
        x2 = 0.5;
        y1 = interpl_Iratio1(i);
        y2 = interpl_Iratio2(i);
        x  = Iratio(i);
        interpl_y_correction(i) = y1 + (y2-y1)/(x2-x1) * (x-x1);
   
    elseif Iratio(i) > 0.5 && Iratio(i) <= 2
        if zFrac(i) > 0 && zFrac(i) < 0.25
            x1  = 0;
            x2  = 0.25;
            y1  = 1.976;
            y2  = 1.886;
            y11 = 0.514;
            y22 = 0.527;
            
            x = zFrac(i);
            interpl_Iratio1(i) = y1  +   (y2-y1)/(x2-x1) * (x-x1);
            interpl_Iratio2(i) = y11 + (y22-y11)/(x2-x1) * (x-x1);
             
        elseif zFrac(i) > 0.25 && zFrac(i) < 0.5
            x1  = 0.25;
            x2  = 0.5;
            y1  = 1.886;
            y2  = 1.796;
            y11 = 0.527;
            y22 = 0.553;
            
            x = zFrac(i);
            interpl_Iratio1(i) = y1  +   (y2-y1)/(x2-x1) * (x-x1);
            interpl_Iratio2(i) = y11 + (y22-y11)/(x2-x1) * (x-x1);
        
        elseif zFrac(i) > 0.5 && zFrac(i) < 0.75
            x1  = 0.5;
            x2  = 0.75;
            y1  = 1.796;
            y2  = 1.727;
            y11 = 0.553;
            y22 = 0.576;
            
            x = zFrac(i);
            interpl_Iratio1(i) = y1  +   (y2-y1)/(x2-x1) * (x-x1);
            interpl_Iratio2(i) = y11 + (y22-y11)/(x2-x1) * (x-x1);
       
        elseif zFrac(i) > 0.75 && zFrac(i) < 1
            x1  = 0.75;
            x2  = 1;
            y1  = 1.727;
            y2  = 1.682;
            y11 = 0.576;
            y22 = 0.595;
            
            x = zFrac(i);
            interpl_Iratio1(i) = y1  +   (y2-y1)/(x2-x1) * (x-x1);
            interpl_Iratio2(i) = y11 + (y22-y11)/(x2-x1) * (x-x1);    
        
        elseif zFrac(i) == 0.25
            interpl_Iratio1(i) = 1.886;
            interpl_Iratio2(i) = 0.527;
        
        elseif zFrac(i) == 0.5    
            interpl_Iratio1(i) = 1.796;
            interpl_Iratio2(i) = 0.553;
     
        elseif zFrac(i) == 0.75
            interpl_Iratio1(i) = 1.727;
            interpl_Iratio2(i) = 0.576;
           
        elseif zFrac(i) == 1    
            interpl_Iratio1(i) = 1.682;
            interpl_Iratio2(i) = 0.595;
           
        end
        
        x1 = 0.5;
        x2 = 2;
        y1 = interpl_Iratio1(i);
        y2 = interpl_Iratio2(i);
        x  = Iratio(i);
        interpl_y_correction(i) = y1 + (y2-y1)/(x2-x1) * (x-x1);
        
    elseif Iratio(i) > 2 && Iratio(i) <= 4
        if zFrac(i) > 0 && zFrac(i) < 0.25
            x1  = 0;
            x2  = 0.25;
            y1  = 0.514;
            y2  = 0.527;
            y11 = 0.249;
            y22 = 0.276;
            
            x = zFrac(i);
            interpl_Iratio1(i) = y1  +   (y2-y1)/(x2-x1) * (x-x1);
            interpl_Iratio2(i) = y11 + (y22-y11)/(x2-x1) * (x-x1);
        
         elseif zFrac(i) > 0.25 && zFrac(i) < 0.5
            x1  = 0.25;
            x2  = 0.5;
            y1  = 0.527;
            y2  = 0.553;
            y11 = 0.276;
            y22 = 0.303;
            
            x = zFrac(i);
            interpl_Iratio1(i) = y1  +   (y2-y1)/(x2-x1) * (x-x1);
            interpl_Iratio2(i) = y11 + (y22-y11)/(x2-x1) * (x-x1);
        
        elseif zFrac(i) > 0.5 && zFrac(i) < 0.75
            x1  = 0.5;
            x2  = 0.75;
            y1  = 0.553;
            y2  = 0.576;
            y11 = 0.303;
            y22 = 0.330;
            
            x = zFrac(i);
            interpl_Iratio1(i) = y1  +   (y2-y1)/(x2-x1) * (x-x1);
            interpl_Iratio2(i) = y11 + (y22-y11)/(x2-x1) * (x-x1);
       
        elseif zFrac(i) > 0.75 && zFrac(i) < 1
            x1  = 0.75;
            x2  = 1;
            y1  = 0.576;
            y2  = 0.595;
            y11 = 0.330;
            y22 = 0.354;
            
            x = zFrac(i);
            interpl_Iratio1(i) = y1  +   (y2-y1)/(x2-x1) * (x-x1);
            interpl_Iratio2(i) = y11 + (y22-y11)/(x2-x1) * (x-x1);
       
        elseif zFrac(i) == 0.25
            interpl_Iratio1(i) = 0.527;
            interpl_Iratio2(i) = 0.276;
        
        elseif zFrac(i) == 0.5    
            interpl_Iratio1(i) = 0.553;
            interpl_Iratio2(i) = 0.303;
     
        elseif zFrac(i) == 0.75
            interpl_Iratio1(i) = 0.576;
            interpl_Iratio2(i) = 0.330;
           
        elseif zFrac(i) == 1    
            interpl_Iratio1(i) = 0.595;
            interpl_Iratio2(i) = 0.354;
           
        end
        
        x1 = 2;
        x2 = 4;
        y1 = interpl_Iratio1(i);
        y2 = interpl_Iratio2(i);
        x  = Iratio(i);
        interpl_y_correction(i) = y1 + (y2-y1)/(x2-x1) * (x-x1);  

    elseif Iratio(i) > 4 && Iratio(i) <= 8
        if zFrac(i) > 0 && zFrac(i) < 0.25
            x1  = 0;
            x2  = 0.25;
            y1  = 0.249;
            y2  = 0.276;
            y11 = 0.121;
            y22 = 0.143;
            
            x = zFrac(i);
            interpl_Iratio1(i) = y1  +   (y2-y1)/(x2-x1) * (x-x1);
            interpl_Iratio2(i) = y11 + (y22-y11)/(x2-x1) * (x-x1);
             
         elseif zFrac(i) > 0.25 && zFrac(i) < 0.5
            x1  = 0.25;
            x2  = 0.5;
            y1  = 0.276;
            y2  = 0.303;
            y11 = 0.143;
            y22 = 0.165;
            
            x = zFrac(i);
            interpl_Iratio1(i) = y1  +   (y2-y1)/(x2-x1) * (x-x1);
            interpl_Iratio2(i) = y11 + (y22-y11)/(x2-x1) * (x-x1);
        
        elseif zFrac(i) > 0.5 && zFrac(i) < 0.75
            x1  = 0.5;
            x2  = 0.75;
            y1  = 0.303;
            y2  = 0.330;
            y11 = 0.165;
            y22 = 0.188;
            
            x = zFrac(i);
            interpl_Iratio1(i) = y1  +   (y2-y1)/(x2-x1) * (x-x1);
            interpl_Iratio2(i) = y11 + (y22-y11)/(x2-x1) * (x-x1);
       
        elseif zFrac(i) > 0.75 && zFrac(i) < 1
            x1  = 0.75;
            x2  = 1;
            y1  = 0.330;
            y2  = 0.354;
            y11 = 0.188;
            y22 = 0.210;
            
            x = zFrac(i);
            interpl_Iratio1(i) = y1  +   (y2-y1)/(x2-x1) * (x-x1);
            interpl_Iratio2(i) = y11 + (y22-y11)/(x2-x1) * (x-x1);
       
        elseif zFrac(i) == 0.25
            interpl_Iratio1(i) = 0.276;
            interpl_Iratio2(i) = 0.143;
        
        elseif zFrac(i) == 0.5    
            interpl_Iratio1(i) = 0.303;
            interpl_Iratio2(i) = 0.165;
     
        elseif zFrac(i) == 0.75
            interpl_Iratio1(i) = 0.330;
            interpl_Iratio2(i) = 0.188;
           
        elseif zFrac(i) == 1    
            interpl_Iratio1(i) = 0.354;
            interpl_Iratio2(i) = 0.210;
           
        end
        
        x1 = 4;
        x2 = 8;
        y1 = interpl_Iratio1(i);
        y2 = interpl_Iratio2(i);
        x  = Iratio(i);
        interpl_y_correction(i) = y1 + (y2-y1)/(x2-x1) * (x-x1); 
    
    elseif Iratio(i) > 8
        if zFrac(i) > 0 && zFrac(i) <= 0.25
            interpl_Iratio1(i) = 0.9877*(Iratio(i)^-0.997);
            interpl_Iratio2(i) = 0.9893*(Iratio(i)^-0.925);

            y1 = interpl_Iratio1(i);
            y2 = interpl_Iratio2(i);
            x1 = 0;  
            x2 = 0.25;
            x  = zFrac(i);

            interpl_y_correction(i) = y1 + (y2-y1)/(x2-x1) * (x-x1); 
          
        elseif zFrac(i) > 0.25 && zFrac(i) <= 0.5         
            interpl_Iratio1(i) = 0.9893*(Iratio(i)^-0.925);
            interpl_Iratio2(i) = 0.9878*(Iratio(i)^-0.855);

            y1 = interpl_Iratio1(i);
            y2 = interpl_Iratio2(i);
            x1 = 0.25; 
            x2 = 0.5;  
            x  = zFrac(i);

            interpl_y_correction(i) = y1 + (y2-y1)/(x2-x1) * (x-x1); 
             
        elseif zFrac(i) > 0.5 && zFrac(i) <= 0.75
            interpl_Iratio1(i) = 0.9878*(Iratio(i)^-0.855);
            interpl_Iratio2(i) = 0.9918*(Iratio(i)^-0.796);

            y1 = interpl_Iratio1(i);
            y2 = interpl_Iratio2(i);
            x1 = 0.5;   
            x2 = 0.75;
            x  = zFrac(i);
                
            interpl_y_correction(i) = y1 + (y2-y1)/(x2-x1) * (x-x1); 
       
        elseif zFrac(i) > 0.75 && zFrac(i) <= 1
            interpl_Iratio1(i) = 0.9918*(Iratio(i)^-0.796);
            interpl_Iratio2(i) = 1.0002*(Iratio(i)^-0.75);
            
            y1 = interpl_Iratio1(i);
            y2 = interpl_Iratio2(i);
            x1 = 0.75;
            x2 = 1;
            x  = zFrac(i);
            
            interpl_y_correction(i) = y1 + (y2-y1)/(x2-x1) * (x-x1); 
           
        end
    end
end

disp_CF = interpl_y_correction;

end % function dispCorrectionFactor

