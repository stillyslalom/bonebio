classdef Brep2D
    % Brep Version 2013.01
    % Version 2013.02
    % - Added code to compute area
    % - Added myArea private variable 
    % Version 2013.03
    % - Added distOfPointsToBrep
    properties(GetAccess = 'public', SetAccess = 'private')
        myBrep; % See brepformat.txt
        myNumDivisionsPerArcSegment;
        myMinArcSegments;
        myPDEGeom;
        myDebug;
        mySegmentMapping;
        myArea;
    end
    methods(Static)
        function brep = readBrep2D(fileName)
            fid = fopen(fileName,'r');
            nPoints = fscanf(fid,' %d ',1);
            brep.vertices = fscanf(fid,' %f ',[2, nPoints]);
            nSegments = fscanf(fid,' %d ',1);
            brep.segments = fscanf(fid,' %f ',[4, nSegments]);
            fclose(fid);
        end
    end
    methods
        function obj = Brep2D(brep)
            % brep can be a brep structure or a file name
            if (ischar(brep))
                obj.myBrep =  obj.readBrep2D(brep);
            else
                obj.myBrep = brep;
            end
            obj.myArea = obj.brepArea();
            obj.myNumDivisionsPerArcSegment = 50; % affects meshing on arcs, so create a fine division
            obj = obj.convertBreptoPdeGeom(obj.myNumDivisionsPerArcSegment);
        end
        function obj = resetBrep(obj,brep)
            obj.myBrep = brep; 
            obj.myArea = obj.brepArea();
            obj = obj.convertBreptoPdeGeom(obj.myNumDivisionsPerArcSegment);
        end
        function obj =  convertBreptoPdeGeom(obj,minArcSegments)
            % convert from brep format to Matlab's pdegeom format
            if (nargin == 1)
                minArcSegments = 2;
            end
            brep = obj.myBrep;
            v = brep.vertices;
            nSegments = size(brep.segments,2);
            pdeGeom = zeros(11,1);
            segmentMapping = [];
            npdeGeomSegments = 0;
            dl = 1; dr = 0;
            for i = 1:nSegments
                breptype = brep.segments(1,i);
                vs = brep.segments(2,i);
                ve = brep.segments(3,i);
                
                xs = v(1,vs);xe = v(1,ve);
                ys = v(2,vs);ye = v(2,ve);
                if (breptype == 1) % line
                    pdetype = 2;
                    npdeGeomSegments = npdeGeomSegments+1;
                    pdeGeom(:,npdeGeomSegments) = [pdetype xs xe ys ye dl dr 0 0 0 0];
                    segmentMapping(npdeGeomSegments) = i;
                elseif (breptype == 2) % arc 
                    pdetype = 1;
                    vc = abs(brep.segments(4,i));
                    xc = v(1,vc); yc  = v(2,vc);
                    R = sqrt((xc-xs)^2+(yc-ys)^2);
                    % # of segments generated will depend on hmax
                    d = sqrt((xe-xs)^2+(ye-ys)^2);
                    theta = real(2*asin(d/2/R)); % due to floating point errors, sometimes d > 2*R resulitng in complex values
                    arcLength = R*theta;
                    if (nargin == 2)
                        numArcSegments = minArcSegments; % use prescribed arc length
                    else
                        hmax = obj.myElemSize;
                        numArcSegments = max(ceil(0.5*arcLength/hmax),minArcSegments);
                    end                   
                    vecStart = [xs-xc ys-yc];
                    vecStart = vecStart/norm(vecStart);
                    vecEnd = [xe-xc ye-yc];
                    vecEnd = vecEnd/norm(vecEnd);
                    thetaStart = atan2(vecStart(2),vecStart(1));
                    thetaEnd = atan2(vecEnd(2),vecEnd(1));
                    if (brep.segments(4,i) > 0) % cw
                        if (thetaEnd > thetaStart) % not allowed
                            thetaEnd = thetaEnd - 2*pi;
                        end
                    else % ccw
                        if (thetaEnd < thetaStart) % not allowed
                            thetaEnd = thetaEnd + 2*pi;
                        end
                    end
                    for j = 1:numArcSegments
                        npdeGeomSegments = npdeGeomSegments+1;
                        segmentMapping(npdeGeomSegments) = i; %#ok<*AGROW>
                        theta0 = thetaStart+(j-1)*(thetaEnd-thetaStart)/numArcSegments;
                        theta1= thetaStart+j*(thetaEnd-thetaStart)/numArcSegments;
                        xs = xc + R*cos(theta0); xe = xc + R*cos(theta1);
                        ys = yc + R*sin(theta0); ye = yc + R*sin(theta1);
                        if (brep.segments(4,i) > 0) % clockwise
                            pdeGeom(:,npdeGeomSegments) = [pdetype xe xs ye ys dr dl xc yc R R];
                        else % anti-clockwise
                            pdeGeom(:,npdeGeomSegments) = [pdetype xs xe ys ye dl dr xc yc R R];
                        end
                    end
                elseif (breptype == 3) % circle ... start and end vertices are same
                    pdetype = 1; % arc
                    vc = abs(brep.segments(4,i));
                    xc = v(1,vc); yc  = v(2,vc);
                    R = sqrt((xc-xs)^2+(yc-ys)^2);
                    
                    % # of segments generated will depend on hmax
                    arcLength = R*2*pi;
                    if (nargin == 2) % use prescribed minarcsegments
                        numArcSegments = 4*minArcSegments;
                    else
                        hmax = obj.myElemSize;
                        numArcSegments = max(ceil(arcLength/hmax),4*minArcSegments);
                    end
                    for j = 1:numArcSegments
                        npdeGeomSegments = npdeGeomSegments+1;
                        segmentMapping(npdeGeomSegments) = i;
                        thetaStart = (j-1)*2*pi/numArcSegments;
                        thetaEnd = j*2*pi/numArcSegments;
                        
                        xs = xc + R*cos(thetaStart); xe = xc + R*cos(thetaEnd);
                        ys = yc + R*sin(thetaStart); ye = yc + R*sin(thetaEnd);
                        pdeGeom(:,npdeGeomSegments) = [pdetype xs xe ys ye dl dr xc yc R R];
                    end
                end
            end
            obj.myPDEGeom = pdeGeom;
            obj.mySegmentMapping = segmentMapping;    
        end
        function [dMin,closestPts] = distOfPointsToBrep(obj,pts)
            brep = obj.myBrep;
            v = brep.vertices;
            nSegments = size(brep.segments,2);
            nPts = size(pts,2);
            dMin = 1e12*ones(1,nPts);
            closestPts = pts; % default
            for i = 1:nSegments
                breptype = brep.segments(1,i);
                vs = brep.segments(2,i);
                ve = brep.segments(3,i);
                if (breptype == 1) % line
                    [dSeg,closestPtsSeg] = obj.distOfPointsToLineSegment(pts,v(:,vs),v(:,ve));
                elseif (breptype == 2) % arc 
                    vc = abs(brep.segments(4,i));
                    [dSeg,closestPtsSeg] = obj.distOfPointsToArcSegment(pts,v(:,vs),v(:,ve),v(:,vc));
                    closestPtsSeg
                end
                for n = 1:nPts
                    if (dSeg(n) < dMin(n))
                        dMin(n) = dSeg(n);
                        closestPts(:,n) = closestPtsSeg(:,n);       
                    end
                end
            end
        end
        function [d,closestPts] = distOfPointsToLineSegment(~,pts,lineStart,lineEnd)
            % vectorized operation
            ABS_TOL = 1e-10;
            L = norm(lineStart-lineEnd);
            closestPts = pts;
            if (L < ABS_TOL) 
                d = norm(pts-lineStart);
                closestPts(1,:) = lineStart(1);
                closestPts(2,:) = lineStart(2);
                return;
            end
            lineTangent = (lineEnd-lineStart)/L;
            v1(1,:) = pts(1,:) - lineStart(1);
            v1(2,:) = pts(2,:) - lineStart(2);      
            v2(1,:) = pts(1,:) - lineEnd(1);
            v2(2,:) = pts(2,:) - lineEnd(2);
            numer = (lineEnd-lineStart)'*v1;   
            u = numer/L^2;
            distToStart = norm(v1);
            distToEnd= norm(v2);
            distToLineSeg = abs((-v1(1,:)*lineTangent(2) + v1(2,:)*lineTangent(1)));
            d = (u < 0).*distToStart + (u >1).*distToEnd + (u>=0).*(u<= 1).*distToLineSeg;
            closestPts(1,:) = (u < 0).*lineStart(1) + (u >1).*lineEnd(1) + ...
                            (u>=0).*(u<= 1).*((1-u).*lineStart(1) + u.*lineEnd(1));  
            closestPts(2,:) = (u < 0).*lineStart(2) + (u >1).*lineEnd(2) + ...
                            (u>=0).*(u<= 1).*((1-u).*lineStart(2) + u.*lineEnd(2));   
        end
        function [d,closestPts] = distOfPointsToArcSegment(~,pts,arcStart,arcEnd,arcCenter)
            % vectorized operation
            ABS_TOL = 1e-10;
            L = norm(arcStart-arcEnd);
            if (L < ABS_TOL)
                d = norm(pts-arcStart);
                closestPts(1,:) = arcStart(1);
                closestPts(2,:) = arcStart(2);
                return;
            end
            startAngle = atan2(arcStart(2)-arcCenter(2),arcStart(1)-arcCenter(1));
            endAngle = atan2(arcEnd(2)-arcCenter(2),arcEnd(1)-arcCenter(1));
            angle = atan2(pts(2,:)-arcCenter(2),pts(1,:)-arcCenter(1));
            r = norm(arcCenter-arcStart);
            closestPts = pts;
            if ((startAngle >= 0) && (endAngle < 0))
                endAngle = endAngle + 2*pi;
                angle = (angle+2*pi).*(angle < 0) + angle.*(angle >= 0); 
            end
            % need to  loop over all pts
            for n = 1:size(pts,2);
                diffAngle = (angle(n)-startAngle)*(angle(n)-endAngle);
                if (diffAngle < 0) % pt lies inbetween
                    distToCenter = norm(arcCenter-pts(:,n));
                    if (distToCenter >= r)
                        distToArc = distToCenter - r;
                    else
                        distToArc = r-distToCenter;
                    end
                    unitVec = (pts(:,n)-arcCenter)/distToCenter;
                    closestPts(:,n) = arcCenter + unitVec*r;
                else %// find the shorter of two distances
                    distToStart = norm(arcStart-pts(:,n));
                    distToEnd = norm(arcEnd-pts(:,n));
                    if (distToStart <= distToEnd)
                        distToArc = distToStart;
                        closestPts(:,n) = arcStart;
                    else
                        distToArc = distToEnd;
                        closestPts(:,n) = arcEnd;
                    end
                end
                d(n) = distToArc;
            end
        end
        function adjustFigScale(~)
            axis auto;axis equal;axis tight;
            v = axis;
            xRange = v(2)-v(1);
            yRange = v(4)-v(3);
            deltaX = 0.1*xRange;
            deltaY = 0.1*yRange;
            v = [v(1)-deltaX v(2)+deltaX  v(3)-deltaY v(4)+deltaY];
            axis(v);
        end
        function plotGeometry(obj,noNumbering)
            if (nargin == 1)
                noNumbering = 0;
            end
            pdeGeom = obj.myPDEGeom;
            segmentMapping = obj.mySegmentMapping;
            pdegplot(pdeGeom);
            obj.adjustFigScale();
            if (noNumbering),return;end
            
            segments = 1:size(pdeGeom,2);
            for seg = segments
                brepSegment = segmentMapping(seg);
                xmid = mean(pdeGeom(2:3,seg));
                ymid = mean(pdeGeom(4:5,seg));
                text(xmid,ymid,num2str(brepSegment),'fontsize',10);
            end
            hold on;
            p = obj.myBrep.vertices;
            for i = 1:size(p,2)
                plot(p(1,i),p(2,i),'b*');
                text(p(1,i),p(2,i),num2str(i),'fontsize',12);
            end
            obj.adjustFigScale();
        end
        function saveFigureToWordDoc(~,fileName)
            [fpath,fname,fext] = fileparts(fileName);
            if isempty(fpath); fpath = pwd; end
            if isempty(fext); fext = '.doc'; end
            fileName = fullfile(fpath,[fname,fext]);
            print -dbitmap;
            word = actxserver('Word.Application');
            if ~exist(fileName,'file');
                op = invoke(word.Documents,'Add');
            else
                op = invoke(word.Documents,'Open',fileName);
            end
            end_of_doc = get(word.activedocument.content,'end');
            set(word.application.selection,'Start',end_of_doc);
            set(word.application.selection,'End',end_of_doc);
            invoke(word.Selection,'Paste');
            if ~exist(fileName,'file')
                invoke(op,'SaveAs',fileName,1);
            else
                invoke(op,'Save');
            end
            invoke(op,'Close');
            invoke(word,'Quit');
            delete(word);
        end
        function area = areaUnderLineSegment(~,startPt,endPt)
            base = (startPt(1)-endPt(1));
            a1 = base*startPt(2);
            a2 = 0.5*(base)*(endPt(2)-startPt(2));
            area = a1 + a2;
        end
        
        function area = areaUnderArcSegment(obj,startPt,endPt,centerPt)
            a1 = obj.areaUnderLineSegment(centerPt,endPt);
            a2 = obj.areaUnderLineSegment(centerPt,startPt);
            R = norm(endPt-centerPt);
            d = norm(endPt-startPt);
            theta = real(2*asin(d/2/R));
            a3 = theta*R^2/2;
            area = a1 - a2 - a3;
        end
        
        function area = brepArea(obj)
            brep = obj.myBrep;
            v = brep.vertices;
            nSegments = size(brep.segments,2);
            area = 0;
            for i = 1:nSegments
                breptype = brep.segments(1,i);
                vs = brep.segments(2,i);
                ve = brep.segments(3,i);
                if (breptype == 1) % line
                    area = area + obj.areaUnderLineSegment(v(:,vs),v(:,ve));
                elseif (breptype == 2) % arc 
                    vc = abs(brep.segments(4,i));
                    area = area + obj.areaUnderArcSegment(v(:,vs),v(:,ve),v(:,vc));
                end
            end
        end
        function Area = trimeshArea(~,p,t)
            nTriangles = size(t,2);
            Area = 0;
            for elem = 1:nTriangles
                nodes = t(1:3,elem)';
                xNodes = p(1,nodes);
                yNodes = p(2,nodes);
                invJ = [(-yNodes(1)+yNodes(3)) (-yNodes(2)+yNodes(1)); ...
                (-xNodes(3)+xNodes(1)) (-xNodes(1)+xNodes(2))];
                dJ = invJ(1,1)*invJ(2,2)-invJ(1,2)*invJ(2,1);
                Area = Area + dJ/2;
            end
        end
        function obj = setDebugOn(obj)
            obj.myDebug = 1;
        end
        function obj = setDebugOff(obj)
            obj.myDebug = 0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function drawArrow(~,start,stop,clr,thickness)
            if (nargin == 3)
                clr = 'r';
                thickness = 1.0;
            end
            for i = 1:size(start,1)
                x0 = start(i,1);
                y0 =start(i,2);
                x1 = stop(i,1);
                y1 = stop(i,2);
                
                h = plot([x0 x1],[y0, y1],clr);
                set(h,'linewidth',thickness);
                hold on;
                p = stop(i,:)-start(i,:);
                alpha = 0.1;  % Size of arrow head relative to the length of the vector
                beta = 0.35;  % Width of the base of the arrow head relative to the length
                hu = [x1-alpha*(p(1)+beta*(p(2)+eps)); x1; x1-alpha*(p(1)-beta*(p(2)+eps))];
                hv = [y1-alpha*(p(2)-beta*(p(1)+eps)); y1; y1-alpha*(p(2)+beta*(p(1)+eps))];
                h = plot(hu(:),hv(:),clr);  % Plot arrow head
                set(h,'linewidth',thickness);
            end
        end
        function [u, niter, flag,residualNorm] = CG(~,A,f,tolRes,maxiter,u0)
            %  Conjugate Gradient method.
            % Solve Au = f, with initial guess of u0
            % Terminate if any of the following conditions are satisfied
            %(1) if relative residual norm is less than tolRes
            %(2) if iter exceeds maxiter
            if (nargin == 3)
                tolRes = 1e-13;
                maxiter = 2000;
                u0 = 0*f;
            end
            if (nargin == 4), u0 = 0*f; end
            u = u0;         % Set u_0 to the start vector s
            r = f - A*u;   % Compute first residuum
            D = diag(A);
            a = r./D;
            p = a;
            rho = a'*r;
            niter = 1;     % Init counter for number of iterations
            normf = norm(f);
            residualNorm = zeros(1,maxiter);
            while (1)   % Test break condition
                residualNorm(niter+1) = sqrt(r'*r);
                a = A*p;
                alpha = rho/(a'*p);
                u = u + alpha*p;
                r = r - alpha*a;
                a = r./D;
                rhonew = a'*r;
                beta = rhonew/rho;
                p = a + beta * p;
                rho = rhonew;
                niter = niter + 1;
                normr = norm(r);
                
                if (normr/normf < tolRes)
                    flag = 1;
                    break
                end
                if (niter >= 2) && isnan(normr)
                    flag = -1;
                    return
                end
                if (niter == maxiter) || (normr/normf > 1e5)
                    flag = -1;                   % is reached, break.
                    break;
                end
            end
            residualNorm = residualNorm(1:niter);
        end
    end
end