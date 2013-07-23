function mi_debby_modified()
%--------------------------------XXXX--------------------------------------
%calculate area MI of bone with a composite materials model using 20 materials.
%Input variable is a text file. File should contain 4 column vectors:- 
%x,y,z coods. and HU for each pixel of CT scanned and segmented bone image. 
% No output variables. Data is saved in a text file; filename given by user
% Saved file has a header line. 
% Delimiter for header and data columns is space (" "). 
%--------------------------------XXXX--------------------------------------


% Test change
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%                              INPUT PHASE
% INPUT 1: FILENAME
% INPUT 2: CT image orientation angle compared to testing/table orientation
% INPUT 3: PIXEL SIZE
% INPUT 4: DIRECTION OF CT IMAGE ROTATION (CLOCKWISE/COUNTERCLOCKWISE)
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


filename = uigetfile('*.*','Select the file containing bone slice data');

angle = input('ENTER VALUE OF ROTATION ANGLE IN DEGREES AND PRESS THE ENTER KEY       :');

angle = pi()*angle/180;   %Convert angle from degrees to radians


A = input('PIXEL size (in mm) AND PRESS ENTER KEY      :');

A = A*A;  % area of each pixel in mm^2


disp(' ');
disp(' ');
disp(' CHECK ORIENTATION OF CT SCAN IMAGE COMPARED TO TESTING        :');
disp(' ');
disp('IF CT IMAGE IS ROTATED CLOCKWISE COMPARED TO TESTING, ENTER "1" IN NEXT STEP');
disp('IF CT IMAGE IS ROTATED COUNTERCLOCKWISE COMPARED TO TESTING, ENTER "2" IN NEXT STEP');
disp(' ');
disp(' ');
orientation = input('TYPE "1" OR "2" FOR TESTING ORIENTATION, THEN PRESS ENTER KEY       :');

disp(' ');
disp(' ');

%XXXXXXXXXXXXXXXXXXXXXXXXXX   END INPUT PHASE   XXXXXXXXXXXXXXXXXXXXXXXXXXX



B = dlmread(filename);           % data file is saved as Matrix B
SIZEB = size(B,1);           % # of rows in matrix B. 
aprho = (0.7303*B(:,4)-1.4169)/1000; %app. density calculations acc. to "UW-Mimics"
E = 2065*((aprho).^3.09);    % elastic modulus acc. to "Wirtz 2000"
[min_E] = min(E);            % identify pixel with minimum value of modulus
B(:,5) = (E)/min_E;    % new column showing relative modulus w.r.t min. E
[max_ratio] = max(B(:,5));
[min_ratio] = min(B(:,5));


%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%create basis for dividing bone into 20 materials based on relative modulus
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

ratio_range = (max_ratio-min_ratio)/20;

% divide bone into 20 materials based on relative modulus

for i = 1:SIZEB
    
   if ((B(i,5) >=1)&&(B(i,5)<(1+ratio_range)))
        
        B(i,6) = 1;
            
   elseif ((B(i,5) >=(1+ratio_range))&&(B(i,5)<(1+2*ratio_range)))
            B(i,6) = 1+ratio_range;
            
   elseif ((B(i,5) >=(1+2*ratio_range))&&(B(i,5)<(1+3*ratio_range)))
            B(i,6) = 1+2*ratio_range;
          
   elseif ((B(i,5) >=(1+3*ratio_range))&&(B(i,5)<(1+4*ratio_range)))
            B(i,6) = 1+3*ratio_range;
          
   elseif ((B(i,5) >=(1+4*ratio_range))&&(B(i,5)<(1+5*ratio_range)))
            B(i,6) = 1+4*ratio_range;
           
   elseif ((B(i,5) >=(1+5*ratio_range))&&(B(i,5)<(1+6*ratio_range)))
            B(i,6) = 1+5*ratio_range;
           
   elseif ((B(i,5) >=(1+6*ratio_range))&&(B(i,5)<(1+7*ratio_range)))
            B(i,6) = 1+6*ratio_range;
            
   elseif ((B(i,5) >=(1+7*ratio_range))&&(B(i,5)<(1+8*ratio_range)))
            B(i,6) = 1+7*ratio_range;
           
   elseif ((B(i,5) >=(1+8*ratio_range))&&(B(i,5)<(1+9*ratio_range)))
            B(i,6) = 1+8*ratio_range;
           
   elseif ((B(i,5) >=(1+9*ratio_range))&&(B(i,5)<(1+10*ratio_range)))
            B(i,6) = 1+9*ratio_range;
             
   elseif ((B(i,5) >=(1+10*ratio_range))&&(B(i,5)<(1+11*ratio_range)))
            B(i,6) = 1+10*ratio_range;
          
   elseif ((B(i,5) >=(1+11*ratio_range))&&(B(i,5)<(1+12*ratio_range)))
            B(i,6) = 1+11*ratio_range;
            
   elseif ((B(i,5) >=(1+12*ratio_range))&&(B(i,5)<(1+13*ratio_range)))
            B(i,6) = 1+12*ratio_range;
            
  elseif ((B(i,5) >=(1+13*ratio_range))&&(B(i,5)<(1+14*ratio_range)))
            B(i,6) = 1+13*ratio_range;
        
  elseif ((B(i,5) >=(1+14*ratio_range))&&(B(i,5)<(1+15*ratio_range)))
            B(i,6) = 1+14*ratio_range;
            
  elseif ((B(i,5) >=(1+15*ratio_range))&&(B(i,5)<(1+16*ratio_range)))
            B(i,6) = 1+15*ratio_range;
           
  elseif ((B(i,5) >=(1+16*ratio_range))&&(B(i,5)<(1+17*ratio_range)))
            B(i,6) = 1+16*ratio_range;
            
  elseif ((B(i,5) >=(1+17*ratio_range))&&(B(i,5)<(1+18*ratio_range)))
            B(i,6) = 1+17*ratio_range;
            
  elseif ((B(i,5) >=(1+18*ratio_range))&&(B(i,5)<(1+19*ratio_range)))
            B(i,6) = 1+18*ratio_range;
            
   else 
            B(i,6) = (1+19*ratio_range);
                  
    end
end
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%                   Calculation of MI and CS area
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

j = 0; k = 1;  % define and initialize counter variables

% j = slice counter; i,k,l,m = pixel counters


while k < SIZEB  %main loop to go thru all pixel information in all slices
    
    j = j+1;            % slice counter
    
    % following variables are for each slice
    i = k;              % define and initialise counter for centroid
                        % and area for each slice. 
    n = k;                    
    l = k;              % define counter for coordinate transformation
    z_coord = B(i,3);     % slice defined by value of 'z' coordinate
    E_sum = 0;            % initialize denominator value for centroid calc. 
    x_sum = 0;            % initialize x cood. numerator value for centroid 
    y_sum = 0;            % Initialize y-cood. numerator value for centroid
    Ix = 0;               % initialize MI value about x axis
    Iy = 0;               % initialize MI value about y axis
    Ixy_pr = 0;           % Initialize Product of Inertia 
    ymaxplus = 0;         % intialize  (+) 'C' value along y direction
    ymaxminus = 0;        % initialize (-) 'C' value along y direction
    ymax_prime_plus = 0;  % initialize (+) 'C' value along yprime direction
    x_prime_plus = 0;     % initialize corresponding xprime coordinate
    ymax_prime_minus = 0; % initialize (-) 'C' value along yprime direction
    x_prime_minus = 0;    % initialize corresponding xprime coordinate
    csarea = 0;           % initialize Cross Sectional area 
    
        
    % loop calculates centroid coordinates and C.S. area for each slice
    
    while (B(i,3) == z_coord)  
        
        E_sum = E_sum+B(i,6);
        x_sum = x_sum+(B(i,1)*B(i,6));
        y_sum = y_sum+(B(i,2)*B(i,6));
        csarea = csarea+1;
        
        if i == SIZEB
            break;
        end
        i = i+1;    
    end                    % end of loop calculating C.S. area and centroid

    cxx = x_sum/E_sum;     % x coordinate of centroid  
    cyy = y_sum/E_sum;     % y coordinate of centroid
    
    % loop calculates MI about x and y axes, Product of Inertia (xy)
    % and 'C' values for each slice
    
    while (B(k,3) == z_coord)
        
              
        Ix = Ix + B(k,6)*(B(k,2)-cyy)^2;
        Iy = Iy + B(k,6)*(B(k,1)-cxx)^2;
        Ixy_pr = Ixy_pr+B(k,6)*((B(k,1)-cxx)*(B(k,2)-cyy));                     
        
        if k == SIZEB
            break;
        end
        
        k = k+1; 
        
    end             % end of loop calculating MI
    
    
 % Do not need this commented portion below. Jump to line 249   
 %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
     % loop calculates +/- C values and corresponding n values along y dirn.
   %{ 
    while (B(n,3)==z_coord)
        
        if (B(n,2)-cyy) > 0 % calculate(+)'C'value along y direction
            
            tmp_ymaxplus = B(n,2)-cyy; 
                      
            if tmp_ymaxplus > ymaxplus
                
               ymaxplus = tmp_ymaxplus;
               n_plus = 0;              %initialise modulus ratio 'n'
               %ctr = n;                 %new counter variable for next loop
               
               while ((B(n,2)-cyy) == ymaxplus)
                   
                   if B(n,6) > n_plus
                       
                       n_plus = B(n,6);
                                              
                   end
                   
                   if n == SIZEB
                       break;
                   end
                   
                   n= n + 1;
                   
               end
               
               n = n-1;   %re-adjust counter'n' after last 'ctr'increment
            end
            
        else   % calculate (-) 'C' value along y direction
            
            tmp_ymaxminus = B(n,2)-cyy;
                        
            if tmp_ymaxminus < ymaxminus
                
                ymaxminus = tmp_ymaxminus;
                n_minus = 0;
                %ctr = n;
                
                while ((B(n,2)-cyy) == ymaxminus)
                   
                   if B(n,6) > n_minus
                       
                       n_minus = B(n,6);
                                              
                   end
                   
                   if n == SIZEB
                       break;
                   end
                   
                   n = n + 1;
                   
                end
                n = n - 1;
            end
            
        end
        
        if n == SIZEB
            break;
        end
       
        n = n+1; 
        
    end
         
    %A = 0.124; % area of each pixel in mm^2
    %}
% Do not need above commented portion     
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx%
    
    
% XXXXXXX    Calculating various geometry values for each slice    XXXXXXX
    
    
    c_xx(j,1) = cxx;                     % x-cood of centroid for slice 'j'
    c_yy(j,1) = cyy;                     % y-cood of centroid for slice 'j'
    cs_area(j,1) = A*csarea;             % Cross Sectional area in mm^2
    %Ixx(j,1) = A*Ix;                    % MI about x axis in mm^4
    %Iyy(j,1) = A*Iy;                    % MI about y axis in mm^4
    
    % Ixx = MI about x axis in mm^4. Iyy = MI about y axis in mm^4
    
    Ixx(j,1)= A*((Ix+Iy)/2 + (Ix-Iy)*cos(2*angle)/2 - Ixy_pr*sin(2*angle));
    Iyy(j,1)= A*((Ix-Iy)/2 - (Ix-Iy)*cos(2*angle)/2 + Ixy_pr*sin(2*angle)); 
    
    Ixy(j,1) = A*Ixy_pr;          % Product of Inertia (xy) in mm^4
    ymax(j,1) = ymaxplus;         % (+) 'C' value in y dir. for "jth" slice
    ymin(j,1) = ymaxminus;        % (-) 'C' value in y dir. for "jth" slice
    nplus(j,1) = n_plus;
    nminus(j,1) = n_minus;
    
   % theta(j,1) = atan(Ixy(j,1)/Ixx(j,1)); % Orientation of N.Axis w.r.t 
                                          % Bending moment axis: (x-axis)
                                          
   %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX%
   
  %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx   
  %x  Transforming coordinates according to orientation for Testing       x     
  %x  X axis should be nominally medial-lateral. Y axis is nominally      x
  %x  Antero-posterior. Rotation angle "angle" is based on angle recorded x
  %x  during segmentation in Mimics                                       x
  %x  Counter-clockwise rotation of angle "Angle" is required to orient    x
  %x  bone from scanned position to testing position 
  %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx   
        
        
    while (B(l,3) == z_coord)  
        
        %x_prime  and y_prime coordinates after shifting origin to centroid
        
        x_prime = (B(l,1)-cxx);
        y_prime = (B(l,2)-cyy);
        
       % B(l,7) = x_prime*cos(theta(j,1)) - y_prime*sin(theta(j,1));  % x_prime coordinate
       % B(l,8) = x_prime*sin(theta(j,1)) + y_prime*cos(theta(j,1));  % y_prime coordinate
       
       if orientation == 1
           % CT image rotated clockwise compared to testing/table
           % orientation
           B(l,7) = x_prime*cos(angle) - y_prime*sin(angle);  % x coordinate oriented for "angle"
           B(l,8) = x_prime*sin(angle) + y_prime*cos(angle);  % y coordinate oriented for "angle"
           
       elseif orientation == 2
           % CT image rotated counterclockwise compared to testing/table
           % orientation
           B(l,7) = x_prime*cos(angle) + y_prime*sin(angle);  % x coordinate oriented for "angle"
           B(l,8) = -x_prime*sin(angle) + y_prime*cos(angle);  % y coordinate oriented for "angle"
       end
       
        if B(l,8) > 0          % calculate(+)'C'value along new Y direction
            
            tmp_ymax_prime_plus = B(l,8); 
                      
            if tmp_ymax_prime_plus > ymax_prime_plus
                
               ymax_prime_plus = tmp_ymax_prime_plus;
               x_prime_plus = B(l,7);
               
            end
            
        else                   % calculate(+)'C'value along new Y direction

            tmp_ymax_prime_minus = B(l,8);
                        
            if tmp_ymax_prime_minus < ymax_prime_minus
                
                ymax_prime_minus = tmp_ymax_prime_minus;
                x_prime_minus = B(l,7);
            end
            
        end
               
        if l == SIZEB
            break;
        end
        
        l = l+1; 
        
    end            % end of loop calculating transformed coordinates                                            
      
    ymax_prime(j,1) = ymax_prime_plus;  % (+) 'C' value in y_prime dir. 
                                         % for "jth" slice
    
    ymin_prime(j,1) = ymax_prime_minus;  % (-) 'C' value in y_prime dir.
                                         % for "jth" slice                                  
                 
                  
    
    xmax_prime(j,1) = x_prime_plus;      % Corresponding xprime coordinate
    xmin_prime(j,1) = x_prime_minus;     % Corresponding xprime coordinate
    
    %nmax_prime(j,1) = n_prime_plus;
    %nmin_prime(j,1) = n_prime_minus;
               
end  % end of outer loop. all calculations completed for all slices 

%XXXXXXXXXXXXXXXX End of Calculation of MI and CS area  XXXXXXXXXXXXXXXXXXX


% Old matrix is commented out. New matrix is below
%      1    2   3   4   5    6   7       8          9          10        11
%                                                                                   12    13     14       15
%M = [c_xx c_yy Ixx Iyy Ixy ymax ymin xmax_prime ymax_prime xmin_prime ymin_prime cs_area nplus nminus  theta]; 

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%      Create matrix with all geometry variables as column vectors        x
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

M = [c_xx c_yy Ixx Iyy ymax_prime ymin_prime cs_area];

% save matrix M as a text file. File name input by user without extension

filename = input('TYPE FILE NAME (WITHOUT EXTENSION) TO SAVE FILE TO (e.g. file1)       :  ','s');
ext = '.txt';
filename = strcat(filename,ext); % creates filename with ".txt" extension

%dlmwrite(filename,M); % data saved as a text file.

% Data saved as a text file. Includes "header line" describing each column.
% When opening in excel, use space (" ") as the delimiter

fid = fopen(filename, 'w');
fprintf(fid, 'Centroid_x_cood Centroid_y_cood MI_x_axis MI_y_axis Cmax Cmin CS_Area\n\n');
fprintf(fid, '%f %f %f %f %f %f %f\n', M);
fclose(fid);

clear all

%xxxxxxxxxxxxxxxxxxxxx       end of function       xxxxxxxxxxxxxxxxxxxxxxxx
