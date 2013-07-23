classdef QuadElasticity < QuadMesher
   % QuadElasticity 
   % Version 2013.01
   % Version 2013.02
   % - resetBrepAndMesh 'bug'
   % Version 2013.03
   % - Line 442: von Mises stress bug
   % - Added setForceVector
    properties(GetAccess = 'public', SetAccess = 'private')
        myTitle;
        myBCtype;
        myBCvalue;
        myDOFPerNode;
        myNumDOF;
        myDOFPerElem;
        myClass;
        myE;
        myNu;
        myD;
        myXi;
        myEta;
        myWt;
        myGradN;
        myK;
        myF;
        myC;
        myfC;
        mySol;
        myU;
        myV;
        myStrainElems;
        myStressElems;
        myVonMisesElems;
        myVonMisesNodes;
        myCompliance;
        myElementArea;
        myFixedDOF; % fixed degrees of freedom
        myFreeDOF; % all free dof
        myForcedNodes; % with non-zero forces
        mySolverMethod;  % solver methods
        myMaxDelta; % maximum displacement
        myMaxStress; % max von Mises stress in the domain at any instance
    end
    methods
        % default: use the parent constructor
        function obj = QuadElasticity(brepFileName,nElements,shape,class)
            obj = obj@QuadMesher(brepFileName,nElements,shape); % call superclass 
            % set the default boundary Conditions
            obj.myBCtype = zeros(obj.myNumBoundarySegments,2);% u v
            obj.myBCvalue = zeros(obj.myNumBoundarySegments,2);% u v
            obj.myDOFPerNode = 2; % (u, v)
            obj.myDOFPerElem = obj.myDOFPerNode*obj.myNodesPerElement;
            obj.myNumDOF = obj.myDOFPerNode*obj.myNumNodes;   
            obj.myE = 2e11;
            obj.myNu = 0.33;  
            E = obj.myE;
            nu = obj.myNu;
            obj.myClass = class;
            obj.myTitle = 'Fem2D';
            obj.mySolverMethod = 2; 
            if (strcmp(obj.myClass,'PlaneStrain'))
                obj.myD = E/((1+nu)*(1-2*nu))*[1-nu nu 0; nu 1-nu 0;0 0 (1-2*nu)/2];
            elseif (strcmp(obj.myClass,'PlaneStress'))
                obj.myD = E/(1-nu^2)*[1 nu 0; nu 1 0;0 0 (1-nu)/2];
            elseif (strcmp(obj.myClass,'AxiSymmetric'))
                obj.myD = E/((1+nu)*(1-2*nu))*[1-nu nu nu 0; nu 1-nu nu 0; nu nu 1-nu 0;0 0 0 (1-2*nu)/2];
            else
                disp('Class undefined ');
                return;
            end
        end
                
        function obj = resetBrepAndSolve(obj,brep)
            if (obj.myDebug), disp('QuadElasticity: resetBrepAndSolve '); end
            obj = obj.resetBrepAndMesh(brep);
            obj.myNumDOF = obj.myDOFPerNode*obj.myNumNodes;
            obj = obj.solveFEProblem(); % Solve Primary FEA problem
        end
        
        function obj = setTitle(obj,title)
            obj.myTitle = title;
        end
        function obj = setYoungsModulus(obj,E)
            obj.myE = E;
        end
        function obj = setPoissonsRatio(obj,nu)
            obj.myNu = nu;
        end
        function obj = setSolverMethod(obj,method)
            obj.mySolverMethod = method;
        end 
        function obj = applyXForceOnEdge(obj,boundaryEdges,force)
            obj.myBCvalue(boundaryEdges,1) = force;
        end
        function obj = applyYForceOnEdge(obj,boundaryEdges,force)
            obj.myBCvalue(boundaryEdges,2) = force;
        end
        function obj = fixXOfEdge(obj,boundaryEdges)
            obj.myBCtype(boundaryEdges,1) = 1;
        end
        function obj = fixYOfEdge(obj,boundaryEdges)
            obj.myBCtype(boundaryEdges,2) = 1;
        end
        function obj = fixEdge(obj,boundaryEdges)
            obj.myBCtype(boundaryEdges,1:2) = 1;
        end
        function obj = setForceVector(obj,f)
            obj.myF = f;
        end
        function obj  = assembleK(obj)
            % Assemble K and f (without boundary conditions) for 2D Linear Elasticity
            % mapping is needed when a subset of elements need to be assembled.
            nDOF = obj.myNumDOF;
            if (strcmp(obj.myShape,'Linear'))
                [xi_GQ,eta_GQ,wt_GQ]= obj.GaussQuad();
            elseif (strcmp(obj.myShape,'Quadratic'))
                disp('Not implemented');
            end
            obj.myXi = xi_GQ;
            obj.myEta = eta_GQ;
            obj.myWt = wt_GQ;
           
            obj.myDOFPerElem = obj.myDOFPerNode*obj.myNodesPerElement;
            nElements = obj.myNumElems;
            nzmax = obj.myDOFPerElem^2*nElements;
            RowTriplets = zeros(nzmax,1);
            ColTriplets = zeros(nzmax,1);
            EntryTriplets = zeros(nzmax,1);
            f = zeros(nDOF,1);
            NCell = cell(1,length(xi_GQ));
            gradNCell = cell(1,length(xi_GQ));
            for i = 1:length(xi_GQ)
                [NCell{i},gradNCell{i}] = obj.QuadShapeFunction(xi_GQ(i),eta_GQ(i));
            end
            obj.myGradN = gradNCell;
            index = 1;
            for elem = 1:nElements
                nodes = obj.myMesh.q(1:obj.myNodesPerElement,elem)';
                [KElem,fElem,AElem] = obj.integrateKOverElem(NCell,gradNCell,elem);     
                if (obj.myPseudoDensity(elem) == 0) % during topology optimization
                    KElem = eye(obj.myDOFPerElem); %avoids singularity in matrix
                    fElem = zeros(obj.myDOFPerElem,1);
                else
                    KElem = obj.myPseudoDensity(elem)^3*KElem;
                    fElem = obj.myPseudoDensity(elem)*fElem;
                end
                obj.myElementArea(elem)= AElem;
                dof = [2*nodes-1; 2*nodes];
                dof = reshape(dof,1,obj.myDOFPerElem);
                temp = dof(ones(1,obj.myDOFPerElem),:);
                colIndex = reshape(temp',1,obj.myDOFPerElem^2);
                rowIndex = reshape(temp,1,obj.myDOFPerElem^2);
                entries = reshape(KElem',1,obj.myDOFPerElem^2);
                RowTriplets(index:index+obj.myDOFPerElem^2-1,1) = rowIndex';
                ColTriplets(index:index+obj.myDOFPerElem^2-1,1) = colIndex';
                EntryTriplets(index:index+obj.myDOFPerElem^2-1,1) = entries';
                index = index+obj.myDOFPerElem^2;
                f(dof) = f(dof) + fElem;
            end
            obj.myK = sparse(RowTriplets,ColTriplets,EntryTriplets,nDOF,nDOF);
            obj.myF = f;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [KElem,fElem,AElem] = integrateKOverElem(obj,NCell,gradNCell,elem)
            nodes = obj.myMesh.q(1:obj.myNodesPerElement,elem)';
            xNodes = obj.myMesh.p(1,nodes);
            KElem = zeros(obj.myDOFPerElem,obj.myDOFPerElem);
            fElem = zeros(obj.myDOFPerElem,1);
            Z = zeros(1,obj.myDOFPerElem/2);
            xi_GQ = obj.myXi;
            eta_GQ = obj.myEta; 
            wt_GQ = obj.myWt;
            D = obj.myD;
            AElem = 0;
            for g = 1:length(wt_GQ)
                N = NCell{g}';
                gradN = gradNCell{g};
                x = xNodes*N;% note: for axisymmetric, x is the radius
                bx = 0;
                by = 0;
                J = obj.Jacobian(elem,xi_GQ(g),eta_GQ(g));
                dJ = (det(J));
                AElem = AElem + dJ;
                T = J'\gradN;
                
                if (strcmp(obj.myClass,'PlaneStrain'))
                    B = [T(1,:) Z; Z T(2,:); T(2,:) T(1,:)];
                    KElem = KElem + wt_GQ(g)*dJ*B'*D*B;
                    fElem = fElem + wt_GQ(g)*dJ*[N*bx;N*by];
                elseif (strcmp(obj.myClass,'PlaneStress'))
                    B = [T(1,:) Z; Z T(2,:); T(2,:) T(1,:)];
                    KElem = KElem + wt_GQ(g)*dJ*B'*D*B;
                    fElem = fElem + wt_GQ(g)*dJ*[N*bx;N*by];
                elseif (strcmp(obj.myClass,'AxiSymmetric'))
                    B = [T(1,:) Z; Z T(2,:);N'/x Z; T(2,:) T(1,:)]; % note: for axisymmetric, x is the radius
                    KElem = KElem + x*wt_GQ(g)*dJ*B'*D*B;
                    fElem = fElem + x*wt_GQ(g)*dJ*[N*bx;N*by];
                end
            end
            if (obj.myDOFPerElem == 8)
                order = reshape([1:4;5:8],1,8); % alternate u and v
            end
            KElem = KElem(order,order);
            fElem = fElem(order);
        end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = assembleBC(obj)
            if (strcmp(obj.myShape,'Linear'))
                [xi_GQ, wt_GQ] = obj.GaussQLine(1);
                N1D = cell(1,length(xi_GQ));
                for i = 1:length(xi_GQ)
                    N1D{i} = obj.edgeShapeFunction(xi_GQ(i));
                end
            elseif (strcmp(obj.myShape,'Quadratic'))
                disp('Not implemented');pause
            end 
            nDOF = obj.myNumDOF;
            % Assemble surface force (Neumann data)
            fBoundary = zeros(nDOF,1);
            obj.myForcedNodes = [];
            for geomEdge = 1:size(obj.myBrep.segments,2)
                typeu = obj.myBCtype(geomEdge,1);
                typev = obj.myBCtype(geomEdge,2);
                valueu = obj.myBCvalue(geomEdge,1);
                valuev = obj.myBCvalue(geomEdge,2);
                if (typeu == 0 ) && (abs(valueu) > 0)
                    boundarySegments = find(obj.myMesh.e(5,:) == geomEdge);
                    for seg = boundarySegments
                        nodes = obj.myMesh.e(1:obj.myNodesPerEdge,seg);
                        udof = 2*nodes-1;
                        fBoundaryElem = obj.integrateOverBoundary(geomEdge,seg,wt_GQ,N1D,1);
                        fBoundary(udof) = fBoundary(udof) +fBoundaryElem;
                        obj.myForcedNodes = unique([obj.myForcedNodes ;nodes]);
                    end
                end
                if (typev == 0 ) && (abs(valuev) > 0)
                    boundarySegments = find(obj.myMesh.e(5,:) == geomEdge);
                    for seg = boundarySegments
                        nodes = obj.myMesh.e(1:obj.myNodesPerEdge,seg);
                        vdof = 2*nodes;
                        fBoundaryElem = obj.integrateOverBoundary(geomEdge,seg,wt_GQ,N1D,2);
                        fBoundary(vdof) = fBoundary(vdof) + fBoundaryElem;
                        obj.myForcedNodes = unique([obj.myForcedNodes ;nodes]);
                    end
                end
            end
            % Gather Dirichlet boundary conditions
            isDirichlet = zeros(nDOF,1);
            dirValue = zeros(nDOF,1);
            for geomEdge = 1:size(obj.myBrep.segments,2)
                typeu = obj.myBCtype(geomEdge,1);
                typev = obj.myBCtype(geomEdge,2);
                valueu = obj.myBCvalue(geomEdge,1);
                valuev = obj.myBCvalue(geomEdge,2);
                if (typeu == 1)
                    boundarySegments = find(obj.myMesh.e(5,:) == geomEdge);
                    for seg = boundarySegments
                        nodes = obj.myMesh.e(1:obj.myNodesPerEdge,seg);
                        nodes = unique(nodes(:));
                        udof = 2*nodes-1;
                        isDirichlet(udof) = 1;
                        dirValue(udof) =  valueu;
                    end
                end
                if (typev == 1)
                    boundarySegments = find(obj.myMesh.e(5,:) == geomEdge);
                    for seg = boundarySegments
                        nodes = obj.myMesh.e(1:obj.myNodesPerEdge,seg);
                        nodes = unique(nodes(:));
                        vdof = 2*nodes;
                        isDirichlet(vdof) = 1;
                        dirValue(vdof) =  valuev;
                    end
                end
            end
            dirichletDOF = find(isDirichlet(:) == 1);
            obj.myFixedDOF = dirichletDOF;
            nDirichletDOF = length(dirichletDOF);
            C = zeros(nDirichletDOF,nDOF);
            for i = 1:nDirichletDOF
                C(i,dirichletDOF(i)) = 1;
            end
            C = sparse(C);
            fDirichlet = dirValue(dirichletDOF);
            obj.myF = obj.myF + fBoundary;
            obj.myC = C;
            obj.myfC = fDirichlet;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function fBoundaryElem = integrateOverBoundary(obj,geomEdge,seg,wt_GQ,N1D,dof)
            nodes = obj.myMesh.e(1:obj.myNodesPerEdge,seg);
            xNodes = obj.myMesh.p(1,nodes);
            yNodes = obj.myMesh.p(2,nodes);
            dx = xNodes(2)-xNodes(1);
            dy = yNodes(2)-yNodes(1);
            L = sqrt(dx^2 + dy^2);
            vec = [dx dy 0]/L;
            zVec = [0 0 1];
            normal = cross(vec,zVec);
            nx = normal(1); %#ok<NASGU>
            ny = normal(2); %#ok<NASGU>
            fBoundaryElem = zeros(numel(N1D{1}),1);
            for g = 1:length(wt_GQ)
                N = N1D{g};
                x = xNodes*N; % radius for axisymmetric problems
                f =  obj.myBCvalue(geomEdge,dof); % dof is either 1 (u) or 2 (v)
                if (strcmp(obj.myClass,'PlaneStrain')) || (strcmp(obj.myClass,'PlaneStress'))
                    fBoundaryElem = fBoundaryElem + wt_GQ(g)*(L/2)*N*f;
                elseif (strcmp(obj.myClass,'AxiSymmetric'))
                    fBoundaryElem = fBoundaryElem + x*wt_GQ(g)*(L/2)*N*f;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj,success] = solveFEProblem(obj)
            if (obj.myDebug), disp('QuadElasticity: solveFEProblem '); end
            obj = obj.assembleK();
            obj = obj.assembleBC();
            [obj,success] = obj.solveLinearSystem();
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj,success] = solveLinearSystem(obj)
            allDOF = 1:obj.myNumDOF;
            method = obj.mySolverMethod;
            success = 2;
            if (method == 1) % direct solve with Lagrange multipliers
                scale = max(max(obj.myK));
                obj.myC = scale*obj.myC;
                obj.myfC = scale*obj.myfC;
                [nDirichletDOF] = size(obj.myC,1);
                % Solution of Algebraic Problem
                Z = spalloc(nDirichletDOF,nDirichletDOF,1);
                KBar = sparse([obj.myK obj.myC'; obj.myC Z]);
                fBar = sparse([obj.myF;  obj.myfC]);
                soln = KBar \ fBar;
            elseif (method == 2)% direct solve with elimination  of fixed dof
                obj.myFreeDOF  = setdiff(allDOF,obj.myFixedDOF);
                % the useful part of the KMatrix
                KTilde = obj.myK(obj.myFreeDOF,obj.myFreeDOF);
                fTilde = obj.myF(obj.myFreeDOF);
                % now subtract all the dirichlet values from rhs
                for i = 1:numel(obj.myFixedDOF)
                    dof = obj.myFixedDOF(i);
                    if (abs(obj.myfC(i)) > 0)
                        fTilde = fTilde - obj.myK(obj.myFreeDOF,dof)*obj.myfC(i);
                    end
                end
                soln = zeros(obj.myNumDOF,1);
                soln(obj.myFixedDOF) = obj.myfC; % fixed values
                soln(obj.myFreeDOF) =  KTilde\fTilde;
            elseif (method == 3)% iterative solve with elimination of fixed dof
                obj.myFreeDOF  = setdiff(allDOF,obj.myFixedDOF);
                % the useful part of the KMatrix
                KTilde = obj.myK(obj.myFreeDOF,obj.myFreeDOF);
                fTilde = obj.myF(obj.myFreeDOF);
                % now subtract all the dirichlet values from rhs
                for i = 1:numel(obj.myFixedDOF)
                    dof = obj.myFixedDOF(i);
                    if (abs(obj.myfC(i)) > 0)
                        fTilde = fTilde - obj.myK(obj.myFreeDOF,dof)*obj.myfC(i);
                    end
                end
                soln = zeros(obj.myNumDOF,1);
                soln(obj.myFixedDOF) = obj.myfC; % fixed values
                [soln(obj.myFreeDOF),~, flag] =  obj.CG(KTilde,fTilde);
                if (flag == -1)
                    success = -1;
                    return;
                end
            end
            nDOF = size(obj.myK,1);
            obj.mySol = full(soln(1:nDOF));
            obj.myU = full(soln(1:2:nDOF));
            obj.myV = full(soln(2:2:nDOF));
            obj.myMaxDelta = max(sqrt(obj.myU.^2 + obj.myV.^2));
            obj.myCompliance = (obj.mySol)'*obj.myF;
            obj = obj.computeStresses();
            obj.myMaxStress = max(obj.myVonMisesElems);
        end
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plotDeformation(obj)
            scale = 0.1*obj.myModelScale/obj.myMaxDelta;
            X = zeros(obj.myNumElems,5);
            Y = zeros(obj.myNumElems,5);
            for count = 1:obj.myNumElems
                if (obj.myPseudoDensity(count) == 0), continue;end
                nodes = obj.myMesh.q(:,count);
                X(count,:) = obj.myMesh.p(1,[nodes' nodes(1)]) + scale*(obj.myU([nodes' nodes(1)]))';
                Y(count,:) = obj.myMesh.p(2,[nodes' nodes(1)])+ scale*(obj.myV([nodes' nodes(1)]))';
            end
            plot( X',Y','g','LineWidth',0.2); hold on;% quad
            axis equal; axis on;view(2);
            obj.adjustFigScale();
            hold on;
            title(['Max Deformation = ' num2str(obj.myMaxDelta)]);
        end     
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plotStress(obj)
            X = zeros(obj.myNumElems,5);
            Y = zeros(obj.myNumElems,5);
            Z = zeros(obj.myNumElems,5);
            nodalField = obj.myVonMisesNodes;
            nodalField(isnan(nodalField)) = 0;   
            for count = 1:obj.myNumElems
                if (obj.myPseudoDensity(count) == 0), continue;end
                nodes = obj.myMesh.q(:,count);
                X(count,:) = obj.myMesh.p(1,[nodes' nodes(1)]);
                Y(count,:) = obj.myMesh.p(2,[nodes' nodes(1)]);
                Z(count,:) = nodalField([nodes' nodes(1)]);
                %Z(count,:) = obj.myVonMisesElems(count);
            end
            fill( X', Y',Z','EdgeColor','none') % surface
            axis equal; axis on;view(2);
            obj.adjustFigScale();
            hold on;
            title(['Max Stress = ' num2str(max(nodalField))]);
        end    
         function obj = computeStresses(obj)
            % Compute stresses at the center of element
            xi = 0;
            eta = 0;
            nElements = obj.myNumElems;
            nNodes = obj.myNumNodes;
            [~,gradN] = obj.QuadShapeFunction(xi,eta); % gradient at center
            obj.myStrainElems = zeros(nElements,2,2);
            obj.myStressElems = zeros(nElements,2,2);
            obj.myVonMisesElems = zeros(1,nElements);
            obj.myVonMisesNodes = zeros(1,nNodes);
            nElemsConnectedToNode  = zeros(1,nNodes);
            D = obj.myD;
            for elem = 1:nElements
                if (obj.myPseudoDensity(elem) == 0)
                    continue;
                end
                nodes = obj.myMesh.q(1:obj.myNodesPerElement,elem)';
                nElemsConnectedToNode(nodes) = nElemsConnectedToNode(nodes) + 1;
                uvalue = obj.myU(nodes);
                vvalue = obj.myV(nodes);
                J = obj.Jacobian(elem,xi,eta);
                B = J'\gradN;
                gradu = B*uvalue;
                gradv = B*vvalue;
                ux = gradu(1);
                uy = gradu(2);
                vx = gradv(1);
                vy = gradv(2);       
                obj.myStrainElems(elem,:,:) = [ux (uy+vx)/2; (uy+vx)/2 vy];
                sxx = D(1,1)*obj.myStrainElems(elem,1,1) + D(1,2)*obj.myStrainElems(elem,2,2);
                syy = D(2,1)*obj.myStrainElems(elem,1,1) + D(2,2)*obj.myStrainElems(elem,2,2);
                sxy = 2*D(3,3)*obj.myStrainElems(elem,1,2);
                obj.myStressElems(elem,:,:) = [sxx sxy; sxy syy];      
                obj.myVonMisesElems(elem) =  sqrt(sxx*sxx + syy*syy - sxx*syy + 3*sxy*sxy);
                obj.myVonMisesNodes(nodes) = obj.myVonMisesNodes(nodes) + obj.myVonMisesElems(elem);
            end
            obj.myVonMisesNodes =  obj.myVonMisesNodes./nElemsConnectedToNode;    
         end
        function plotMesh(obj)
            obj.plotMesh@QuadMesher();
            hold on;
            % mark all fixed nodes
            index = (rem(obj.myFixedDOF,2)==1);
            fixedXNodes = (obj.myFixedDOF(index)-1)/2+1;        
            plot(obj.myMesh.p(1,fixedXNodes),obj.myMesh.p(2,fixedXNodes),'xk');
            index = (rem(obj.myFixedDOF,2)==0);
            fixedYNodes = (obj.myFixedDOF(index))/2;
            plot(obj.myMesh.p(1,fixedYNodes),obj.myMesh.p(2,fixedYNodes),'ok');  
            
            % mark all forced nodes
            nodes = obj.myForcedNodes;
            Fx = obj.myF(2*nodes-1);
            Fy = obj.myF(2*nodes);
            xRange = max(obj.myMesh.p(1,:))-min(obj.myMesh.p(1,:));
            yRange = max(obj.myMesh.p(2,:))-min(obj.myMesh.p(2,:));
            scale = 0.1*max(xRange,yRange);
            normF = sqrt(Fx.^2 + Fy.^2);
            if (normF > 0)
                Fx = Fx./normF;
                Fy = Fy./normF;
                Fx = Fx(:)';
                Fy = Fy(:)';
                plot(obj.myMesh.p(1,nodes),obj.myMesh.p(2,nodes),'dr');
                start = [obj.myMesh.p(1,nodes); obj.myMesh.p(2,nodes)];
                stop = start + scale*[Fx; Fy];
                obj.drawArrow(start',stop');
            end
            axis auto;axis equal;axis tight;
            obj.adjustFigScale();
            str = [obj.myTitle '; #Elems = ' num2str(obj.myNumElems)];
            title(str); 
        end
    end
end