classdef Bonebio5_7_13 < QuadElasticity
    properties(GetAccess = 'public', SetAccess = 'private')
        myAge % model age
        mya_P_OBp % "preosteoblastic proliferation fraction"
        myD_Pivonka_OBu % maximum differentiation rate of osteoblast progenitor cells by Pivonka et al.
        myD_OBu % differentiation rate of osteoblast precursor cells
        myC_OBu % concentration of uncommitted osteoblasts
        myC_OBp % concentration of osteoblast precursors
        myC_OBa % concentration of active osteoblasts
        myC_OBa_Nodes
        myC_OCp % Concentation of osteoclast precursors
        myC_OCa % Conentration of active osteoclasts
        myA_OBa % rate of active osteoblast apoptosis
        myD_OBp % max differentiation rate of osteoblast precursors
%         myP_PTH_d %PTH dosage term
%         myP_PTH_d_Nodes %PTH dosage at node
        
        myP_TGF_d_Nodes %PTH dosage at node        
        myPhi_Nodes       
        mySED
        mymaxSED
        mySED_Nodes 
        myK_res%resorption rate
        myK_form%formation rate
%         myRho% inital density
        mylambda% anabolic strength parameter
        myPhi% current porosity
%         myPi_PTH_ACT_OB% activatior function of PTH
%         myPi_PTH_REP_OB% repressor equlibrium constant
%         myB_PTH % intrinsic PTH production rate
%         myD_PTH % constant degradation rate of PTH
%         myK_PTH_ACT_OB % RANKL production-relevant equilibrium dissociation 
%         %constant related to binding of PTH to its receptors expressed on osteoblasts
%         myK_PTH_REP_OB % OPG production-relevant equilibrium dissociation constant 
%         %related to binding of PTH to its receptors expressed on osteoblasts
%         
        myPi_TGF_ACT_OBu% activatior function of PTH
        myPi_TGF_REP_OBp% repressor equlibrium constant
        myD_TGF_beta % constant degradation rate of PTH
        myK_TGF_ACT_OBu % equilibrium dissociation constant related to binding of TGF to its receptors expressed on uncommitted osteoblasts
        myK_TGF_REP_OBp % equilibrium dissociation constant related to binding of TGF to its receptors expressed on osteoblasts precursors
        myMin_Pi_MECH_act
        myP_TGF_d

    end
    methods
        function obj = Bonebio5_7_13(brepFileName,nElements,shape,class) % <insert pithy comment here>
             obj = obj@QuadElasticity(brepFileName,nElements,shape,class);
        end
        function obj=setCellConc(obj,C_OBuo,C_OBpo,C_OBao,C_OCpo,C_OCao) %set initial bone cell concentrations
            obj.myC_OBu(1:obj.myNumElems) = C_OBuo; % initial molar concentration of uncommited osteoblast cells
            obj.myC_OBp(1:obj.myNumElems) = C_OBpo; % initial molar concentration of osteoblast precursor cells
            obj.myC_OBa(1:obj.myNumElems) = C_OBao; % initial molar concentration of active osteoblast cells
            obj.myC_OCp(1:obj.myNumElems) = C_OCpo; % initial molar concentration of osteoclast precursor cells
            obj.myC_OCa(1:obj.myNumElems) = C_OCao; % initial molar concentration of active osteoclast cells
         end
        function obj=setRates(obj,D_OBp,A_OBa,K_res,Min_Pi_MECH_act,lambda) % set differentiation and reaction rates
            obj.myD_OBp(1:obj.myNumElems) = D_OBp; % max differentiation rate of osteoblast precursor cells
            obj.myA_OBa(1:obj.myNumElems) = A_OBa; % apoptosis rate of active osteoblasts 
            obj.myK_res(1:obj.myNumElems) = K_res; % bone resorption rate 
            obj.myMin_Pi_MECH_act(1:obj.myNumElems) = Min_Pi_MECH_act;
            obj.mylambda(1:obj.myNumElems) = lambda;
        end
 
%         function obj=setPTH(obj,P_PTH_d,B_PTH,D_PTH,K_PTH_ACT_OB,K_PTH_REP_OB) %grand master PTH concentration control
% %             if (nargin == 6)
% %                 members = 1:obj.myNumElems;
% %             else
% %                 assert(max(members) <= obj.myNumElems);
% %                 assert(min(members) >=  1);
% %             end
%             obj.myB_PTH(1:obj.myNumElems)=B_PTH; %intrinsc PTH production rate
%             obj.myD_PTH(1:obj.myNumElems)=D_PTH; %constant degradation rate of PTH
%             obj.myK_PTH_ACT_OB(1:obj.myNumElems)=K_PTH_ACT_OB; %RANKL production-relevant 
%             %equilibrium dissociation constant related to binding of PTH to its receptors expressed on osteoblasts
%             obj.myK_PTH_REP_OB(1:obj.myNumElems)=K_PTH_REP_OB; %OPG production-relevant 
%             %equilibrium dissociation constant related to binding of PTH to its receptors expressed on osteoblasts
%             obj.myP_PTH_d(1:obj.myNumElems) =P_PTH_d.*rand(1,obj.myNumElems); %PTH dosage
%             nElements=obj.myNumElems;
%             for elem=1:nElements %assigning PTH concentrations to each element
%             C_PTHElem=(obj.myB_PTH(elem)+obj.myP_PTH_d(elem))/obj.myD_PTH(elem);% molar concentration of PTH
%             obj.myPi_PTH_ACT_OB(elem)=C_PTHElem/(obj.myK_PTH_ACT_OB(elem)+C_PTHElem); %activator function of carrying capacity of RANKL due to PTH
%             obj.myPi_PTH_REP_OB(elem)=C_PTHElem/(obj.myK_PTH_REP_OB(elem)+C_PTHElem); %repressor function of OPG production following PTH
%             end
%         end
        
        function obj=setTGFact(obj,a_P_OBp,D_Pivonka_OBu,D_TGF_beta,K_TGF_ACT_OBu,K_TGF_REP_OBp,P_TGF_d) %grand master TGF concentration control
            obj.mya_P_OBp(1:obj.myNumElems) = a_P_OBp; % "preosteoblastic proliferation fraction"
            obj.myD_Pivonka_OBu(1:obj.myNumElems) = D_Pivonka_OBu; % maximum differentiation rate of osteoblast progenitor cells by Pivonka et al.
            obj.myD_OBu(1:obj.myNumElems) = (1-obj.mya_P_OBp).*obj.myD_Pivonka_OBu; % maximum differentiation rate of uncommited osteoblasts
            obj.myD_TGF_beta(1:obj.myNumElems) = D_TGF_beta; %constant degradation rate of PTH
            obj.myK_TGF_ACT_OBu(1:obj.myNumElems) = K_TGF_ACT_OBu; % equilibrium dissociation constant related to binding of TGF-beta to its receptors expressed on uncommitted osteoblasts
            obj.myK_TGF_REP_OBp(1:obj.myNumElems) = K_TGF_REP_OBp; % equilibrium dissociation constant related to binding of TGF-beta to its receptors expressed on osteoblasts precursors
            obj.myP_TGF_d(1:obj.myNumElems) = P_TGF_d.*rand(1,obj.myNumElems); %TGF dosage
            nElements = obj.myNumElems;
            for elem=1:nElements %assigning TGF concentrations to each element
            C_TGF_betaElem = obj.myP_TGF_d(elem); % molar concentration of PTH **SIMPLIFICATION TO REMOVE C_OCa DEPENANCE**
            obj.myPi_TGF_ACT_OBu(elem) = C_TGF_betaElem/(obj.myK_TGF_ACT_OBu(elem)+C_TGF_betaElem); %activator function of carrying capacity of RANKL due to PTH
            obj.myPi_TGF_REP_OBp(elem) = obj.myK_TGF_REP_OBp(elem)/(obj.myK_TGF_REP_OBp(elem)+C_TGF_betaElem); %repressor function of OPG production following PTH
            end
        end       
        function obj=setPorosity(obj,Phi) %change porosity of bone
            obj.myPhi(1:obj.myNumElems)=Phi; %bone porosity
        end
        function obj = simulateBone(obj)
            i=0;
            while i<=20
                obj=obj.solveFEProblem();
                obj=obj.solveSEDProblem();
                figure(1)
                obj.plotPhi(); 
                pause(10)
                close(1)
                figure(2)
                obj.plotStress();
                pause(10)
                close(2)
                figure(3)
                obj.plotSED();
                pause(10)
                close(3)
                figure(4)
                obj.plotC_OBa();
                pause(10)
                close(4)
                obj=obj.assemblePhi;                         
                i=i+1;
            end
        end
        function obj = solveSEDProblem(obj)
            StressX = obj.myStressElems(:,1,1);
            StressY = obj.myStressElems(:,2,2);
            StressXY = obj.myStressElems(:,1,2);
            strainX = obj.myStrainElems(:,1,1);
            strainY = obj.myStrainElems(:,2,2);
            strainXY = obj.myStrainElems(:,1,2);
            obj.mySED = ((1/2)*((StressX.*strainX)+(StressY.*strainY)...
                    +(StressXY.*strainXY)))'./obj.myElementArea;
            obj.mymaxSED = max(obj.mySED);
%             devStrainX = (2/3)*strainX-(1/3)*strainY;
%             devStrainY =(-1/3)*strainX+(2/3)*strainX;
%             gamXY = 2*(strainXY);           
%             strain = 2/3*sqrt(1.5*(devStrainX.^2+devStrainY.^2)+.75*gamXY.^2);
        end
        function obj = assemblePhi(obj)
            nElements = obj.myNumElems;
            for elem = 1:nElements
                C_OBuElem = obj.myC_OBu(elem); %concentration of uncommitted osteoblast cells
                P_OBpElem = (obj.myD_Pivonka_OBu(elem)*C_OBuElem*obj.myPi_TGF_ACT_OBu(elem))/(obj.myC_OBp(elem)*obj.myMin_Pi_MECH_act(elem));
                Pi_MECH_ACT = obj.myMin_Pi_MECH_act(elem)*(1+(obj.mylambda(elem)*((obj.mySED(elem)/obj.mymaxSED)-1)));
                deltaC_OBpElem = (obj.myD_OBu(elem)*C_OBuElem*obj.myPi_TGF_ACT_OBu(elem)) + (P_OBpElem*obj.myC_OBp(elem)*Pi_MECH_ACT) - (obj.myD_OBp(elem)*obj.myC_OBp(elem)*obj.myPi_TGF_REP_OBp(elem));
                C_OBpElem = obj.myC_OBp(elem) + deltaC_OBpElem; %decay of molar concentration of osteoblast precursor cells in each element           
                deltaC_OBaElem = (obj.myD_OBp(elem).*C_OBpElem.*obj.myPi_TGF_REP_OBp(elem)) - obj.myA_OBa(elem).*obj.myC_OBa(elem); %change in active osteoblasts
                obj.myC_OBa(elem) = obj.myC_OBa(elem) + deltaC_OBaElem; %new value of active osteoblast concentration
                if obj.mySED(elem) <= .2*obj.mymaxSED
                    obj.myC_OCa(elem) = 1.05*obj.myC_OCa(elem); %concentration of active osteoclasts in each element
                else
                    obj.myC_OCa(elem) = obj.myC_OCa(elem);
                end
                obj.myK_form(elem) = (obj.myK_res(elem)*obj.myC_OCa(elem))/obj.myC_OBa(elem); %bone formation rate
                deltaPhiElem = obj.myK_res(elem).*obj.myC_OCa(elem) - 6*obj.myK_form(elem).*obj.myC_OBa(elem); %change in porosity
                obj.myPhi(elem) = obj.myPhi(elem)+deltaPhiElem; %new porosity
                if obj.myPhi(elem)<0
                    obj.myPhi(elem)=0;
                elseif obj.myPhi(elem)>1
                    obj.myPhi(elem)=1;
                else 
                    obj.myPhi(elem)=obj.myPhi(elem);
                end
                obj.myPseudoDensity(elem) = 1-obj.myPhi(elem); %bone density as a result of porosity                    
            end            
        end
%         function plotPTH(obj) %plot PTH 
%             obj.myP_PTH_d_Nodes = zeros(1,obj.myNumNodes); %PTH at note
%             for elem = 1:obj.myNumElems
%                if (obj.myPseudoDensity(elem) == 0)
%                     continue;
%                end
%                nodes = obj.myMesh.q(1:obj.myNodesPerElement,elem)'; %mesh
%                obj.myP_PTH_d_Nodes(nodes) = obj.myP_PTH_d_Nodes(nodes) + obj.myP_PTH_d(elem); 
%             end
%          end
            
          function plotTGF(obj) %plot TGF 
            obj.myP_TGF_d_Nodes = zeros(1,obj.myNumNodes); %TGF at nodes
            nElemsConnectedToNode  = zeros(1,obj.myNumNodes);
            for elem = 1:obj.myNumElems
               if (obj.myPseudoDensity(elem) == 0)
                    continue;
               end
               nodes = obj.myMesh.q(1:obj.myNodesPerElement,elem)'; %mesh
               nElemsConnectedToNode(nodes) = nElemsConnectedToNode(nodes) + 1;
               obj.myP_TGF_d_Nodes(nodes) = (obj.myP_TGF_d_Nodes(nodes) + obj.myP_TGF_d(elem)); 
            end
            obj.myP_TGF_d_Nodes=obj.myP_TGF_d_Nodes./nElemsConnectedToNode;
            X = zeros(obj.myNumElems,5);
            Y = zeros(obj.myNumElems,5);
            Z = zeros(obj.myNumElems,5);
%             nodalField = obj.myP_PTH_d_Nodes;
%             nodalField(isnan(nodalField)) = 0;   
            
            nodalField = obj.myP_TGF_d_Nodes;
            nodalField(isnan(nodalField)) = 0;
            
            for count = 1:obj.myNumElems
                if (obj.myPseudoDensity(count) == 0), continue;end
                nodes = obj.myMesh.q(:,count);
                X(count,:) = obj.myMesh.p(1,[nodes' nodes(1)]);
                Y(count,:) = obj.myMesh.p(2,[nodes' nodes(1)]);
                Z(count,:) = nodalField([nodes' nodes(1)]);

            end
            fill( X', Y',Z','EdgeColor','none') % surface
            axis equal; axis on;view(2);
            obj.adjustFigScale();
            hold on;
            title(['TGF Conc = ' num2str(max(nodalField))]);
        end    
        function plotSED(obj) %plot SED 
            obj.mySED_Nodes = zeros(1,obj.myNumNodes); %TGF at nodes
            nElemsConnectedToNode  = zeros(1,obj.myNumNodes);
            for elem = 1:obj.myNumElems
               if (obj.myPseudoDensity(elem) == 0)
                    continue;
               end
               nodes = obj.myMesh.q(1:obj.myNodesPerElement,elem)'; %mesh
               nElemsConnectedToNode(nodes) = nElemsConnectedToNode(nodes) + 1;
               obj.mySED_Nodes(nodes) = (obj.mySED_Nodes(nodes) + obj.mySED(elem)); 
            end
            obj.mySED_Nodes=obj.mySED_Nodes./nElemsConnectedToNode;
            X = zeros(obj.myNumElems,5);
            Y = zeros(obj.myNumElems,5);
            Z = zeros(obj.myNumElems,5);
%             nodalField = obj.myP_PTH_d_Nodes;
%             nodalField(isnan(nodalField)) = 0;   
            
            nodalField = obj.mySED_Nodes;
            nodalField(isnan(nodalField)) = 0;
            
            for count = 1:obj.myNumElems
                if (obj.myPseudoDensity(count) == 0), continue;end
                nodes = obj.myMesh.q(:,count);
                X(count,:) = obj.myMesh.p(1,[nodes' nodes(1)]);
                Y(count,:) = obj.myMesh.p(2,[nodes' nodes(1)]);
                Z(count,:) = nodalField([nodes' nodes(1)]);

            end
            fill( X', Y',Z','EdgeColor','none') % surface
            axis equal; axis on;view(2);
            obj.adjustFigScale();
            hold on;
            title(['Max Strain Energy Density = ' num2str(max(nodalField))]);
        end 
        function plotPhi(obj) %plot Phi 
            obj.myPhi_Nodes = zeros(1,obj.myNumNodes); %Phi at nodes
            nElemsConnectedToNode  = zeros(1,obj.myNumNodes);
            for elem = 1:obj.myNumElems
               if (obj.myPseudoDensity(elem) == 0)
                    continue;
               end
               nodes = obj.myMesh.q(1:obj.myNodesPerElement,elem)'; %mesh
               nElemsConnectedToNode(nodes) = nElemsConnectedToNode(nodes) + 1;
               obj.myPhi_Nodes(nodes) = (obj.myPhi_Nodes(nodes) + obj.myPhi(elem)); 
            end
            obj.myPhi_Nodes=obj.myPhi_Nodes./nElemsConnectedToNode;
            X = zeros(obj.myNumElems,5);
            Y = zeros(obj.myNumElems,5);
            Z = zeros(obj.myNumElems,5);
  
            nodalField = obj.myPhi_Nodes;
            nodalField(isnan(nodalField)) = 0;
            
            for count = 1:obj.myNumElems
                if (obj.myPseudoDensity(count) == 0), continue;end
                nodes = obj.myMesh.q(:,count);
                X(count,:) = obj.myMesh.p(1,[nodes' nodes(1)]);
                Y(count,:) = obj.myMesh.p(2,[nodes' nodes(1)]);
                Z(count,:) = nodalField([nodes' nodes(1)]);
            end
            fill( X', Y',Z','EdgeColor','none') % surface
            axis equal; axis on;view(2);
            obj.adjustFigScale();
            hold on;
            title(['Phi Max = ' ,num2str(max(nodalField)),'Phi Min = ', num2str(min(nodalField))]);
        end
        function plotC_OBa(obj) %plot concentration of active osteoblasts 
            obj.myC_OBa_Nodes = zeros(1,obj.myNumNodes); %C_OBa at nodes
            nElemsConnectedToNode  = zeros(1,obj.myNumNodes);
            for elem = 1:obj.myNumElems
               if (obj.myPseudoDensity(elem) == 0)
                    continue;
               end
               nodes = obj.myMesh.q(1:obj.myNodesPerElement,elem)'; %mesh
               nElemsConnectedToNode(nodes) = nElemsConnectedToNode(nodes) + 1;
               obj.myC_OBa_Nodes(nodes) = (obj.myC_OBa_Nodes(nodes) + obj.myC_OBa(elem)); 
            end
            obj.myC_OBa_Nodes=obj.myC_OBa_Nodes./nElemsConnectedToNode;
            X = zeros(obj.myNumElems,5);
            Y = zeros(obj.myNumElems,5);
            Z = zeros(obj.myNumElems,5);
%             nodalField = obj.myP_PTH_d_Nodes;
%             nodalField(isnan(nodalField)) = 0;   
            
            nodalField = obj.myC_OBa_Nodes;
            nodalField(isnan(nodalField)) = 0;
            
            for count = 1:obj.myNumElems
                if (obj.myPseudoDensity(count) == 0), continue;end
                nodes = obj.myMesh.q(:,count);
                X(count,:) = obj.myMesh.p(1,[nodes' nodes(1)]);
                Y(count,:) = obj.myMesh.p(2,[nodes' nodes(1)]);
                Z(count,:) = nodalField([nodes' nodes(1)]);

            end
            fill( X', Y',Z','EdgeColor','none') % surface
            axis equal; axis on;view(2);
            obj.adjustFigScale();
            hold on;
            title(['Max Concentration of Active Osteoblasts = ' num2str(max(nodalField))]);
        end 
    end
end

