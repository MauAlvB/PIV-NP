!************************************************************************************************************************************
!************************************************************************************************************************************
!                                    PARTICLE IMAGE VELOCIMETRY - NUMERICAL PARTICLE - SATURATION (PIV-NP-Sr) (v.2024.04.17)
!                                           
! Copyright (c) 2017 Núria M. Pinyol and Mauricio Alvarado. Centre de Metodes Numerics en Enginyeria (CIMNE)
!
! All rights reserved
!
! LICENSE BSD4
! 
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met: 
!
! 1. Redistributions of source code must retain the above copyright notice, this list of 
! conditions and the following disclaimer. 
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of 
! conditions and the following disclaimer in the documentation and/or other materials 
! provided with the distribution. 
!
! 3. All advertising materials mentioning features or use of this software must display the 
! following acknowledgement: 
!
! This product includes PIV-NP published in Pinyol, N.M. & Alvarado, M. (2017) Novel 
! PIV-based analysis for large displacement. Canadian Geotechnical Journal 54(7): 933-944.
! 
! 4. Neither the name of the copyright holder nor the names of its contributors may be used to 
! endorse or promote products derived from this software without specific prior written 
! permission. 
! 
! THIS SOFTWARE IS PROVIDED BY COPYRIGHT HOLDER "AS IS" AND ANY 
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL COPYRIGHT 
! HOLDER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
! USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
! AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN 
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
! POSSIBILITY OF SUCH DAMAGE.

!************************************************************************************************************************************
!************************************************************************************************************************************
!                                    PARTICLE IMAGE VELOCIMETRY - NUMERICAL PARTICLE (PIV-NP)
!      
! This code provide the accumulated displacements, strains and other variables associated to Numerical particles from the data provided 
! by PIV-LAB. PIV-LAB data consists in the velocity measured in cornes of a structured rectangular cells comparing subsequent images elaspsed 
! an increment of time. 
	
	INCLUDE 'common_PIV-NP.for'
	CHARACTER*80 ARCH_PAR,ARCH_REC,ARCH
	INTEGER HORA(8),CARGA
	CHARACTER (LEN=12) RCK(3)
	
	ARCH_REC=TRIM(ARCH)//'.REC'
	
	CALL PIVLAB_DATA                              !Subroutine to define the case (PREPROCESS)
	
	IF (IREC.EQ.1) THEN
	    IP=0
          TIEMPO=0.0D0
          CALL IMPRES_GiD
	ENDIF
	
	DO IETT=1,TOTAL_STEPS                         !Loop for each PIVLAB file
	  IP=IP+1                                     !Time step number
	  TIEMPO=TIEMPO+DT
	  
	  CALL VELOCIDADES                            !Subroutine to load PIVLAB data
	  CALL CONTOUR                                !Subroutine to correct contour velocity
	  CALL SOLMOV                                 !Subroutine to update position, velocity, deformation, etc
	  
	  MODULO=MOD(IP,IMPPAS)                       !Subroutine to print results in GID (POSTPROCESS)
	  IF (IP.EQ.1.OR.MODULO.EQ.0) THEN
	  	CALL IMPRES_GiD
 	  END IF
      ENDDO
	
      CALL RECOM                                    !Subroutine to save data for a new calculation
	
      CALL DATE_AND_TIME(RCK(1),RCK(2),RCK(3),HORA)
	WRITE (*,2000) HORA(5),':',HORA(6),':',HORA(7),', ANALYSIS FINISHED'
      
2000	FORMAT(I2,A,I2,A,I2,A)
3000	FORMAT(I2,A,I2,A,I2,A,A)
	CLOSE (1)
	STOP                                                             
	END
      
!************************************************************************************************************************************
!                                  Subroutine to read input data from external file
!************************************************************************************************************************************      
	
	SUBROUTINE PIVLAB_DATA

	!Data read: 
      !Block 2  NC:Number of PIV-Lab cells; NN: Number of PIV-Lab nodes; NPC: Number of particle per cell; NFIL: Number of rows;
      !         AXC: Width (x-direction) of PIV-Lab cell (m); AYC: Height (y-direction) of PIV-Lab cell (m)
      !Block 3  DT: Time elapsed betweeen images; TOTAL_STEPS: Total number of PIV-Lab files; IMPPAS: Steps between printed results;
      !         MOISTER: 1 = Read moister files, 0 = no info about moister
      !         V.PIVNP: 1 = Mesh PIV-NP equivalent to PIV-Lab mesh (always calculated), 2 = Nodes of PIV-Lab cells are the central points of the elements of PIV-NP mesh 
      !         V.PIVLAB: 1 = 4 Columns files (OLD PIVLAB VERSION), 5 Columns files (NEW PIVLAB VERSION)
      !         IREC: 1 = load previous data and continue analysis, 0 = new analysis
      !         PTV: 1 = Generate additional numerical particle for tracking purposes (to compare wiht tracking particles in laboratori tests), 0 = No additional particles are generated
      !Block 4  S_DENSITY: Soils density in kg/m3; POROSITY: Soils initial porosity
      !         
	
	
	INCLUDE 'common_PIV-NP.for'
	REAL*8 GAUSS(10)
	CHARACTER*80 TEXTO_USUARIO
	CHARACTER*80 ARCH_PAR,ARCH,ARCH_REC,ARCH_POST_MSH,ARCH_POST_RES
	INTEGER HORA(8),NCH,NCH2,IFI
	CHARACTER (LEN=12) RCK(3)
	
	OPEN (2,FILE='PIV-NP.TXT', FORM='FORMATTED')    !Read in PIV-NP.txt the name of the file of data
	
	READ (2,*) ARCH
	
	ARCH_PAR=TRIM(ARCH)//'.PAR'                     !Trim the file name	
	ARCH_REC=TRIM(ARCH)//'.REC'                     !Trim the file name
	
	CLOSE (2)
	
	OPEN (5,FILE=ARCH_PAR, FORM='FORMATTED')        !Read file of data
	CALL DATE_AND_TIME(RCK(1),RCK(2),RCK(3),HORA)
	WRITE (*,3000) HORA(5),':',HORA(6),':',HORA(7),', LEYENDO DATOS...'
	
	READ (5,2000) TIT                                                      !Block 1
	READ (5,2000) TEXTO_USUARIO                                            !Block 2 (Details of mesh)
	READ(5,*) NC,NN,NPC,NFIL,AXC,AYC
	NCH=NC/NFIL
	
	READ (5,2000) TEXTO_USUARIO                                            !Block 3 (Analysis data)
	READ(5,*) DT,TOTAL_STEPS,IMPPAS,MOISTER,IVERSION,IPIVLAB,ICONTOUR,
     * IREC,ITR
	READ (5,2000) TEXTO_USUARIO                                            !Block 4 (Soil parameters)
	READ(5,*) S_DENSITY,POROSITY
	
	IF (IVERSION.EQ.1) THEN  
	  NP = NC*NPC*NPC
	ELSE !(IVERSION.EQ.2)  
	  NFIL2 = NFIL + 1
	  NCH2 = NCH+1
	  NC2 = NCH2*NFIL2
	  NN2 = (NCH2+1)*(NFIL2+1)
	  NP = NC2*NPC*NPC
	ENDIF
	
	NP0=NP
	
	IF (ITR.NE.0) THEN                              !Additional numerical particles are generated
	  READ (5,2000) TEXTO_USUARIO                                          !Block 5 (Tracking PTV)
	  READ(5,*) DXT,DYT,PTVX1,PTVY1,PTVX2,PTVY2,PTVX3,PTVY3
	  
	  XP(NP+1,1)=DXT*PTVX1
	  XP(NP+1,2)=DYT*PTVY1
	  XP(NP+2,1)=DXT*PTVX2
	  XP(NP+2,2)=DYT*PTVY2
	  XP(NP+3,1)=DXT*PTVX3
	  XP(NP+3,2)=DYT*PTVY3
	  NP = NP + 3
      ENDIF
	  
!	PIV-NP mesh definition from PIV-Lab mesh
	
	NCF(1)=1                        !Id number of the first cell from the bottom left
	NNF(1)=1                        !Id number of the first node from the bottom left
	NNFA(1)=NNF(1)+NCH+1            !Id of the fist node of the second row
	XF(1)=0.d0                      !Coordinate "X" of the first row of nodes (left) 
	YF1 = 0.d0                      !Coordinate "Y" of the first row of nodes (bottom)
	
	DO I=2,NFIL                     
	  NCF(I)=NCF(I-1)+NCH           !Id number of the firts cell of each row
	  NNF(I)=NNF(I-1)+NCH+1         !Id number of the firts node (botom left) of each row
	  NNFA(I)=NNFA(I-1)+NCH+1       !Id number of the firts node (top left) of each row
	  XF(I)=XF(1)                   !Coordinate "X" of each row of nodes
      ENDDO
	
	
	NCF(NFIL+1)=NC+1
	DO I=1,NFIL+1
	  YF(I)=YF1+(I-1)*AYC             !Coordenada "Y" de cada fila de nodos
	ENDDO
	AXP=AXC/NPC
	AYP=AYC/NPC
	
	IF (IVERSION.EQ.2) THEN          !Additional mesh generated
	  NCF2(1)=1
	  NNF2(1)=1
	  NNFA2(1)=NNF2(1)+NCH2+1
	  XF2(1)= - AXC/2                  !Coordenada "X" de la primera fila de nodos (Izquierda)
	  YFB2 = - AYC/2                   !Coordenada "Y" de la primera fila de nodos (Abajo)
	
	  DO I=2,NFIL2
	      NCF2(I)=NCF2(I-1)+NCH2
	      NNF2(I)=NNF2(I-1)+NCH2+1
	      NNFA2(I)=NNFA2(I-1)+NCH2+1
	      XF2(I)=XF2(1)                     !Coordenada "X" de cada fila de nodos
	  ENDDO
	
	  NCF2(NFIL2+1)=NC2+1
	  DO I=1,NFIL2+1
	      YF2(I)=YFB2+(I-1)*AYC             !Coordenada "Y" de cada fila de nodos
	  ENDDO
	  NB1(1) = 1
	  NB2(1) = NC2
	ENDIF
	
!	Assign material type to PIV-NP mesh elements

	N1(1) = 1                         !Celda donde empieza a asignar material
	N2(1) = NC                        !Celda donde termina de asignar material	
	
!	Calculation of a vector that relates id of PIV-Lab nodes and PIV-NP nodes (PIVLAB --> PIV-NP)
	DO I=1,NCH+1
	  DO J=1,NFIL+1
	      ICONECTIVIDAD(I*(NFIL+2-J)+(J-1)*(I-1))=I+(J-1)*(NCH+1)
	  ENDDO
	ENDDO	

!	Generation and localization of numerical particles depending of the number of particles per element (NPC)

	SELECT CASE(NPC)
	  CASE(1)
	      GAUSS(1)=0.
	  CASE(2)
	      GAUSS(1)=-.5 !-.577350269189626
	      GAUSS(2)=.5 !.577350269189626
	  CASE(3)
	      GAUSS(1)=-.66666666666667 !-.774596669241483
	      GAUSS(2)=0.
	      GAUSS(3)=0.66666666666667 !.774596669241483
	  CASE(4)
	      GAUSS(1)=-.861136311594053
	      GAUSS(2)=-.339981043584856
	      GAUSS(3)=.339981043584856
	      GAUSS(4)=.861136311594053
	  CASE(5)
	      GAUSS(1)=-.906179845988664
	      GAUSS(2)=-.538469310105683
	      GAUSS(3)=0.
	      GAUSS(4)=.538469310105683
	      GAUSS(5)=.906179845988664
	  CASE(6)
	      GAUSS(1)=-.932469514203152
	      GAUSS(2)=-.661209386466265
	      GAUSS(3)=-.238619186083197
	      GAUSS(4)=.238619186083197
	      GAUSS(5)=.661209386466265
	      GAUSS(6)=.932469514203152
	END SELECT
	
!	GENERACION DE PARTÍCULAS	 
	IF (IVERSION.EQ.1.AND.IREC.EQ.0) THEN
	  K=1
	  KP=0
	  DO I=1,NFIL                         !Recorre verticalmente las celdas
		    DO J=NCF(I),NCF(I+1)-1        !Recorre horizontalmente las celdas
	          DO II=1,NFIL
	              DO JJ=N1(II),N2(II)
	                  IF( JJ.EQ.J) THEN
	                      GOTO 40
	                  ENDIF
	              ENDDO
	          ENDDO
	          GOTO 50

40	          DO II=1,NPC
	              DO JJ=1,NPC
	                  KP=KP+1
	                  IF (NPC.LE.10) THEN
	                  XP(KP,1)=XF(I)+(J-NCF(I))*AXC+AXC/2.+GAUSS(JJ)*AXC/2.
	                      XP(KP,2)=YF(I)+AYC/2.+GAUSS(II)*AYC/2.
	                  ELSE
	                      XP(KP,1)=XF(I)+(J-NCF(I))*AXC+(JJ-1)*AXP+ AXP/2.
	                      XP(KP,2)=YF(I)+(II-1)*AYP+ AYP/2.
	                  ENDIF
	              ENDDO
	          ENDDO
50	      ENDDO
	  ENDDO
	ELSEIF (IVERSION.EQ.2.AND.IREC.EQ.0) THEN
	  K=1
	  KP=0
	  KP2=0
	  
	  DO I=1,NFIL2                       !Recorre verticalmente las celdas
		    DO J=NCF2(I),NCF2(I+1)-1        !Recorre horizontalmente las celdas
	          DO II=1,NFIL2
	              DO JJ=NB1(II),NB2(II)
	                  IF( JJ.EQ.J) THEN
	                      GOTO 41
	                  ENDIF
	              ENDDO
	          ENDDO
	          GOTO 51

41	          KP2=KP2+1
	          XP2(KP2,1)=XF2(I)+(J-NCF2(I))*AXC+AXC/2.
	          XP2(KP2,2)=YF2(I)+AYC/2.
	          
	          DO II=1,NPC
	              DO JJ=1,NPC
	                  KP=KP+1
	                  IF (NPC.LE.10) THEN
	                XP(KP,1)=XF2(I)+(J-NCF2(I))*AXC+AXC/2.+GAUSS(JJ)*AXC/2.
	                XP(KP,2)=YF2(I)+AYC/2.+GAUSS(II)*AYC/2.
	                  ELSE
	                XP(KP,1)=XF2(I)+(J-NCF2(I))*AXC+(JJ-1)*AXP+ AXP/2.
	                XP(KP,2)=YF2(I)+(II-1)*AYP+ AYP/2.
	                  ENDIF
	              ENDDO
	          ENDDO
51	      ENDDO
	  ENDDO
	ELSE
!	  RETURN
	ENDIF
	
!	Mass particle and initial energy (CURRENTLY NOT USED. MASS OF EACH PARTICLE IMPOSED TO BE 1)
	DO I=1,NP
	  VVP(I)=AXC*AYC/(NPC*NPC)
	  AMP(I)= 1.0d0   !S_DENSITY*(1-POROSITY)*VVP(I)                 
	  E_Potential(I)=AMP(I)*9.81d0*XP(I,2)
	  E_Cinetic_x(I)=0.d0
	  E_Cinetic_y(I)=0.d0
	  E_Total(I)=E_Potential(I)+E_Cinetic_x(I)+E_Cinetic_y(I)
	ENDDO

!	Initialization
	DO I=1,NP
	  DO J=1,4
	      SIG(I,J)=0.
	      EPS(I,J)=0.
	      DEPS(I,J)=0.
	  ENDDO
	  VP(I,1)=0.d0
	  VP(I,2)=0.d0
	  VEL_X_NODO_V2(I)=0.d0
	  VEL_Y_NODO_V2(I)=0.d0
	  VEL_X_NODO_V(I)=0.d0
	  VEL_Y_NODO_V(I)=0.d0
	  VEL_X_NODO(I)=0.d0
	  VEL_Y_NODO(I)=0.d0
	  UP(I,1)=0.d0
	  UP(I,2)=0.d0
	ENDDO
	CLOSE (5)
!	DO I=1,NN
!	    Node_NaN(I)=1
!    ENDDO
	
!     Loading data if there is a previous case (IREC=1)
	
	OPEN (3,FILE=ARCH_REC, FORM='UNFORMATTED')
	IF (IREC.EQ.1) THEN
	
	    REWIND 3
	    READ(3) NP1
	    READ(3) IVERSION
	    READ(3) ((XP(I,J),J=1,2),I=1,NP1)
	    READ(3) ((UP(I,J),J=1,2),I=1,NP1)
	    READ(3) ((EPS(I,J),J=1,4),I=1,NP1)
	    READ(3) (EPSEQ(I),I=1,NP1)
	    READ(3) (NaN_P(I),I=1,NP1)
	    DO I = 1,NP1
	        UPO(I,1)=UP(I,1)
	        UPO(I,2)=UP(I,2)
	    ENDDO
!	    CLOSE (3)
	ELSE
	    DO I=1,NP1
	        DO J=1,4
	            EPS(I,J)=0.0d0
	            DEPS(I,J)=0.0d0
	        ENDDO
	        ACEP(I,1)=0.0d0
	        ACEP(I,2)=0.0d0
	        EPSEQ(I)=0.0d0
	        NaN_P(I)=0
	    ENDDO
	ENDIF
	
2000  FORMAT(A80)
3000	FORMAT(I2,A,I2,A,I2,A)
	END
	
!************************************************************************************************************************************
!                                              Subrutine to read PIVLAB data
!************************************************************************************************************************************ 
! Read position of the points (PIV-Lab nodes) and velocity vector
      
	SUBROUTINE VELOCIDADES
	USE, INTRINSIC :: IEEE_ARITHMETIC          !Fichero que idenficica "NaN"
	INCLUDE 'common_PIV-NP.for'

	REAL*8 EPSIL,VEL_X(200000),VEL_Y(200000)
	REAL*8 POSCX(200000), POSCY(200000) !Read position of the points (PIV-Lab nodes). Not used in the code
	CHARACTER*80 TEXTO_USUARIO,c_aux,c_aux3
	INTEGER c_aux2,SUM_NAN,IVECTOR(200000)
	
	DIMENSION XN(4),YN(4)
	
	DATA XN /-1.,1.,-1.,1./
	DATA YN /-1.,-1.,1.,1./

	EPSIL=EPSILON(EPSIL)
	
	SUM_NAN = 0
	IF (IVERSION.EQ.1) THEN	
	  DO I=1,NN
	      AM(I)=0

	  ENDDO
	ELSE !(IVERSION.EQ.2)
	  DO I=1,NN2
	      AM(I)=0
	  ENDDO
	ENDIF

!     LAZO SOBRE LAS PARTICULAS
!	CALCULA MASA DEL NODO
	
	DO I=1,NP
	  IVERSIONCASO = 1      !Mesh generated incase of Iversion=1 is selected in subroutine UCELDA
	  CALL UCELDA(I)
	  IF (IDONDE(I).NE.-1) THEN
!     CALCULA LAS FUNCIONES DE FORMA PARA LA PARTÍCULA
	      XCM=XC+AXC/2.d0
	      YCM=YC+AYC/2.d0

	      XX=2.d0*(XP(I,1)-XCM)/AXC
	      YY=2.d0*(XP(I,2)-YCM)/AYC
	      
	      DO J=1,4
	          FN=(1.d0+XX*XN(J))*(1.d0+YY*YN(J))/4.d0      !Función de forma que depende de la posición de la partícula respecto al nodo
	          JJ=IN(J)                                     !Indica los nodos del elemento al que pertenece la partícula
	          AM(JJ)=AM(JJ)+AMP(I)*FN                      !Masa del Nodo JJ	      
	      ENDDO
	  ENDIF
	ENDDO
	
!	LECTURA DE ARCHIVOS TXT CON INFORMACION DE PIV_LAB
	WRITE (c_aux,*)IETT
	c_aux=TRIM(c_aux)
	c_aux=ADJUSTL(c_aux)
	c_aux2=LEN(ADJUSTL(TRIM(c_aux)))-1
	c_aux3(:7)='datos ('
	c_aux3(8:(8+c_aux2))=c_aux
	c_aux3((8+c_aux2+1):)=').TXT'
	
	OPEN (7,FILE=c_aux3, FORM='FORMATTED')
	READ (7,2000) TIT
	READ (7,2000) TEXTO_USUARIO
	READ (7,2000) TEXTO_USUARIO
	IF (IPIVLAB.EQ.1) THEN
		READ (7,*)(POSCX(I),POSCY(I),VEL_X(I),VEL_Y(I),I=1,NN)
	ELSE
		DO I=1,NN
	    	READ (7,*) POSCX(I),POSCY(I),VEL_X(I),VEL_Y(I),IVECTOR(I)
	    ENDDO
	ENDIF
	CLOSE (7)
	
	IF (MOISTER.EQ.1) THEN
		WRITE (c_aux,*)IETT

		c_aux=TRIM(c_aux)
		c_aux=ADJUSTL(c_aux)
		c_aux2=LEN(ADJUSTL(TRIM(c_aux)))-1
		c_aux3(:6)='Moist_'
		c_aux3(7:(7+c_aux2))=c_aux
		c_aux3((7+c_aux2+1):)='.TXT'
	
		OPEN (8,FILE=c_aux3, FORM='FORMATTED')
		READ (8,2000) TIT
		READ (8,*)(POSCX(I),POSCY(I),SMOISTURE(I),SATURATION(I),I=1,NN)
		CLOSE (8)
	ENDIF
	
!	CONVIERTE EL NUMERO DE LOS NODOS DE PIV_LAB A PIV-NP	
	DO I=1,NN
		VEL_X_NODO_V(I)=VEL_X_NODO(I)     !Almacena Velocidad X del paso anterior.
		VEL_Y_NODO_V(I)=VEL_Y_NODO(I)     !Almacena Velocidad Y del paso anterior.
	  
		IF (ieee_is_nan(VEL_X(I)).OR.ieee_is_nan(VEL_Y(I))) THEN      !In case of reading velocity equal to NaN in PIV-Lab, the velocity is equal to zero
	    	VEL_X_NODO(ICONECTIVIDAD(I))=0.d0
	        VEL_Y_NODO(ICONECTIVIDAD(I))=0.d0
	        Node_NaN(ICONECTIVIDAD(I))=1                                    !Indica si el nodo contiene NaN (Si = 1, No = 0)
	    ELSE
	    	VEL_X_NODO(ICONECTIVIDAD(I))=VEL_X(I)
	        VEL_Y_NODO(ICONECTIVIDAD(I))=-VEL_Y(I)
	        Node_NaN(ICONECTIVIDAD(I))=0
		ENDIF
	
		IF (ieee_is_nan(SMOISTURE(I)).OR.ieee_is_nan(SATURATION(I))) THEN      !In case of reading velocity equal to NaN in PIV-Lab, the velocity is equal to zero
	    	SMOISTURE_N2(ICONECTIVIDAD(I))=0.d0
	        SATURATION_N2(ICONECTIVIDAD(I))=0.d0
	        Node_NaN_M(ICONECTIVIDAD(I))=1                                    !Indica si el nodo contiene NaN (Si = 1, No = 0)
		ELSE
	    	IF (SMOISTURE(I).LT.0.0) THEN
	        	SMOISTURE(I)=0.D0
	        ENDIF
			SMOISTURE_N2(ICONECTIVIDAD(I))=SMOISTURE(I)
	        SATURATION_N2(ICONECTIVIDAD(I))=SATURATION(I)
	        Node_NaN_M(ICONECTIVIDAD(I))=0
		ENDIF
	ENDDO
	
	IF (IVERSION.EQ.1) THEN
		DO I=1,NN                                             !Over number of nodes mesh1
!			ICOUNT_NO_NAN(I) = 0
	        ICOUNT_NAN(I) = 0
	        DO J=1,8
	            IF (Node_NaN(I).EQ.0.AND.NODE_CONECT(I,J).NE.0) THEN
	                IF (Node_NaN(NODE_CONECT(I,J)).EQ.1) THEN
	                    ICOUNT_NAN(I) = ICOUNT_NAN(I) + 1       !COUNT NUMBER OF NaN AROUND THE NODE I
	                    Node_NaN2(I) = 2
	                    Node_NaN2(NODE_CONECT(I,J)) = 3
	                ENDIF
	            ENDIF
	        ENDDO
	        
	        AM(I) = 1.0D0
	        IF (Node_NaN(I).EQ.0) THEN
	        	PV(I,1)=VEL_X_NODO(I)*AM(I)                     !Calcula cantidad de movimiento (m*v)
	            PV(I,2)=VEL_Y_NODO(I)*AM(I)
	            APV(I,1)=VEL_X_NODO(I)-VEL_X_NODO_V(I)          
	            APV(I,2)=VEL_Y_NODO(I)-VEL_Y_NODO_V(I)
	            APV(I,1)=APV(I,1)*AM(I)                         !Calcula incremento de cantidad de movimiento [m*(v2-v1)]
	            APV(I,2)=APV(I,2)*AM(I)
	            SMOISTURE_N(I)=SMOISTURE_N2(I)*AM(I)
	            SATURATION_N(I)=SATURATION_N2(I)*AM(I)
	        ENDIF
	    ENDDO
	    
	ELSE                                                  !If IVERSION is not equal to 1
		DO I=1,NN2                                            !Over number of nodes mesh2
	    	AMASSNODE(I)=0.0D0
	        PV(I,1)=0.0D0
	        PV(I,2)=0.0D0
	        APV(I,1)=0.0D0
	        APV(I,2)=0.0D0
	        ICOUNT_NO_NAN(I)=0
	        ICOUNT_NO_NAN_OLD(I)=0
	        SMOISTURE_N(I)=0.0D0
	        SATURATION_N(I)=0.0D0
		ENDDO
	    DO I=1,NN !NC2 = NN                                   !Over number of elements mesh2
			IVERSIONCASO = 2 !Mesh generated in case of Iversion=2 is selected in subroutine UCELDA
	        CALL UCELDA(I)

	
	        IF(IDONDE(I).NE.-1.AND.Node_NaN(I).EQ.0) THEN
	        	DO J=1,4
	            	JJ=IN(J)                                  ! JJ nodes of mesh2
	                F1=0.25d0*1.0D0 !AMASSINI(I) !AXC*AYC*S_DENSITY*(1-POROSITY) !AM(JJ)                      !Cuidado con la masa
	                AMASSNODE(JJ)=1.0D0 !AMASSNODE(JJ)+0.25d0*AMASSINI(I)
	                PV(JJ,1)= PV(JJ,1)+VEL_X_NODO(I)*F1
	                PV(JJ,2)= PV(JJ,2)+VEL_Y_NODO(I)*F1
	                APV(JJ,1)=APV(JJ,1)+(VEL_X_NODO(I)-VEL_X_NODO_V(I))*F1
	                APV(JJ,2)=APV(JJ,2)+(VEL_Y_NODO(I)-VEL_Y_NODO_V(I))*F1
	                ICOUNT_NO_NAN(JJ)=ICOUNT_NO_NAN(JJ)+1                                              !Count active elements per node
	                SMOISTURE_N(JJ)=SMOISTURE_N(JJ)+SMOISTURE_N2(I)*F1
	                SATURATION_N(JJ)=SATURATION_N(JJ)+SATURATION_N2(I)*F1
	            ENDDO
	        ENDIF
		ENDDO
	    
	ENDIF
	
	MODULO2=MOD(IP,10)
	IF (IP.EQ.1.OR.MODULO2.EQ.0) THEN
	    WRITE (*,*) 'ANALYZING ', c_aux3
	ENDIF
	
2000  FORMAT(A80)
	END
	
!************************************************************************************************************************************
!                                              Subroutine to calculate displacements and deformation
!************************************************************************************************************************************ 
	
	SUBROUTINE SOLMOV
	USE, INTRINSIC :: IEEE_ARITHMETIC          !Fichero que idenficica "NaN"
	INCLUDE 'common_PIV-NP.for'

	REAL*8 EPSIL,FN,E_TOTAL_FULL,MULTIPLICADOR
	INTEGER NODO_NAN,count_active,CONT,SUM_NAN
	DIMENSION XN(4),YN(4)
	DATA XN /-1.,1.,-1.,1./
	DATA YN /-1.,-1.,1.,1./

	EPSIL=EPSILON(EPSIL)
	E_TOTAL_FULL = 0.0D0
	count_active = 0
	
!	Update of displacement and velocity of numerical particles
	
	DO I=1,NP
		IVERSIONCASO = 1  !Mesh generated in case of Iversion=1 is selected in subroutine UCELDA
		CALL UCELDA(I)
		IF(IDONDE(I).NE.-1) THEN
	        
	        XCM=XC+AXC/2.d0
	    	YCM=YC+AYC/2.d0

	    	XX=2*(XP(I,1)-XCM)/AXC    !Posición x local de la partícula respecto al centro geométrico
	    	YY=2*(XP(I,2)-YCM)/AYC    !Posición y local de la partícula respecto al centro geométrico
	
	    	!Inizialization
	    	XPP(I,1)=0.d0
	    	XPP(I,2)=0.d0
	    	VP(I,1)=0.d0
	    	VP(I,2)=0.d0
	    	ACEP(I,1)=0.d0
	    	ACEP(I,2)=0.d0
	    	UPO(I,1)=0.d0
	    	UPO(I,2)=0.d0
	        NaN_P2(I)=0
	    	IF (IP.EQ.1.AND.IREC.EQ.0) THEN       
	    		NaN_P(I)=0        
	      	ENDIF
	      	IF (IP.EQ.1.AND.IREC.EQ.1) THEN
	        	UPO(I,1)=UP(I,1)
	        	UPO(I,2)=UP(I,2)
	      	ENDIF

	      	SUM_NAN = 0
	      	SMOIST_NP(I)=0.d0
	      	SATURA_NP(I)=0.d0
	      
	      ! Compute velocity of numerical particles from node velocity using shape functions
			DO J=1,4 !Local nodes 
	        	FN=(1.d0+XX*XN(J))*(1.d0+YY*YN(J))/4.d0      !Shape functions
	          	JJ=IN(J)
	          
	          	IF(AM(JJ).GE.EPSIL) THEN
	            	F1=FN!/AM(JJ)
	              	F2=FN*DT!/AM(JJ)
	          	ELSE
	              	F1=0.0d0
	              	F2=0.0d0
	          	ENDIF
	            
! ACTUALIZA VELOCIDAD 	
	          	DO K=1,2
	              VP(I,K)=VP(I,K)+PV(JJ,K)*F1                 !Velocidad de la particula
	              ACEP(I,K)=ACEP(I,K)+APV(JJ,K)*F1/DT         !Aceleración de la particula
	              XPP(I,K)=XPP(I,K)+PV(JJ,K)*F2               !Posición local de la particula
	              UPO(I,K)=UPO(I,K)+PV(JJ,K)*F2               !Desplazamiento instantaneo de la particula
	              UP(I,K)=UP(I,K)+PV(JJ,K)*F2                 !Desplazamiento acumulado de la particula
	          	ENDDO
	         
	          	IF (MOISTER.EQ.1) THEN
	              SMOIST_NP(I)=SMOIST_NP(I)+SMOISTURE_N(JJ)*F1
	              SATURA_NP(I)=SATURA_NP(I)+SATURATION_N(JJ)*F1
	          	ENDIF
	            
	            IF (IVERSION.EQ.1) THEN  !Evalúa cuántos nodos tienen asignado NaN 
	                SUM_NAN = SUM_NAN + Node_NaN(JJ)
	            ELSE !(IVERSION.EQ.2)     !Asigna NaN a las partículas ubicadas en elementos que su centro es un punto PIV-Lab igual a NaN 
	                IF (Node_NaN(INDC).EQ.1) THEN       !
	                    IF (IP.EQ.1.AND.IREC.EQ.0) THEN
	                        NaN_P(I)=1
	                    ENDIF
	                    NaN_P2(I)=1
	                ENDIF
	            ENDIF               
	        ENDDO
	    
	        IF (IVERSION.EQ.1.AND.SUM_NAN.EQ.4) THEN     !Asigna NaN a las partículas ubicadas en elementos en el que los cuatro nodos tienen asignado NaN
	            IELEMENT_ACTIVE(INDC) = 0
	            IF (IP.EQ.1.AND.IREC.EQ.0) THEN
	                NaN_P(I)=1
	            ELSE
	                NaN_P2(I)=1
	            ENDIF      
	        ENDIF      
	          
	        IF (NaN_P2(I).EQ.0) THEN
	          count_active = count_active + 1
              ENDIF
	    ENDIF
	    IF (NaN_P2(I).EQ.1) THEN
	        IELEMENT_ACTIVE(INDC) = 0
	    ELSE
	        IELEMENT_ACTIVE(INDC) = 1
	    ENDIF
	    IF (NaN_P(I).EQ.1) THEN       !Posiciona las partículas asignadas como NaN en el Step 1 a las coordenadas (-1,-1) 
!	      XP(I,1)=-0.01
!	      XP(I,2)=-0.01
	      XPP(I,1)=0.0
	      UP(I,1)=0.0
	      XPP(I,2)=0.0
	      UP(I,2)=0.0
	      IDONDE(I)=-1
	    ENDIF
	ENDDO
	
	IF (ITR.NE.0) THEN
	  NaN_P(NP)=2
	  NaN_P(NP-1)=2
	  NaN_P(NP-2)=2
	ENDIF   

!	CALCULA DEFORMACIONES
	CONT=1
	DO I=1,NP
	  IVERSIONCASO = 1
	  CALL UCELDA(I)
	  IF (IDONDE(I).NE.-1) THEN

	      XCM=XC+AXC/2.d0
	      YCM=YC+AYC/2.d0

	      XX=2.d0*(XP(I,1)-XCM)/AXC
	      YY=2.d0*(XP(I,2)-YCM)/AYC

	      DEPS(I,1)=0.d0
	      DEPS(I,2)=0.d0
	      DEPS(I,3)=0.d0

	      DO J=1,4
	          JJ=IN(J)
	          DX=XN(J)*0.5d0*(1+YY*YN(J))/AXC !XN(J)*(1+YY*YN(J))/4.d0*2.d0/AXC
	          DY=YN(J)*0.5d0*(1+XX*XN(J))/AYC !YN(J)*(1+XX*XN(J))/4.d0*2.d0/AYC
	          DX0=XN(J)*0.5d0/AXC !XN(J)*(1./4.)*2./AXC
	          DY0=YN(J)*0.5d0/AYC !YN(J)*(1./4.)*2./AYC
	          IF(AM(JJ).GE.EPSIL) THEN
	              F1=DX*DT/AM(JJ)
	              F2=DY*DT/AM(JJ)
	              F10=DX0*DT/AM(JJ)
	              F20=DY0*DT/AM(JJ)
	          ELSE
	              F1=0.0
	              F2=0.0
	              F10=0.0
	              F20=0.0
	          ENDIF
	          DEPS(I,1)=DEPS(I,1)+PV(JJ,1)*F10                        !Incremento de deformacion x  !Anterior version: DEPS(I,1)=DEPS(I,1)+PV(JJ,1)*F10  
	          DEPS(I,2)=DEPS(I,2)+PV(JJ,2)*F20                        !Incremento de deformacion y  !Anterior version: DEPS(I,2)=DEPS(I,2)+PV(JJ,2)*F2
	          DEPS(I,3)=DEPS(I,3)+PV(JJ,1)*F20 +PV(JJ,2)*F10           !Incremento de deformacion xy !Anterior version: DEPS(I,3)=DEPS(I,3)+PV(JJ,1)*F2 +PV(JJ,2)*F10
	      ENDDO
	
	      DO J=1,3
	          IF(ABS(DEPS(I,J)).LT.EPSIL.AND.DEPS(I,J).NE.0.) THEN
	              DEPS(I,J)=0.d0
	          ENDIF
	      ENDDO
!              CONT=CONT+1
!            ENDIF
	      
	      EPS(I,1)=EPS(I,1)+DEPS(I,1)                                 !Deformacion x total
	      EPS(I,2)=EPS(I,2)+DEPS(I,2)                                 !Deformacion y total
	      EPS(I,3)=EPS(I,3)+DEPS(I,3)                                 !Deformacion xy total
	      
	      EVOLUMETRIC(I)=EPS(I,1)+EPS(I,2)                            !Deformacion volumetrica total
	      EVOL_IN(I)=DEPS(I,1)+DEPS(I,2)                              !Deformacion volumetrica instantanea
	      E_Potential(I)=AMP(I)*9.81*XP(I,2)
	      E_Cinetic_x(I)=0.5*AMP(I)*VP(I,1)*VP(I,1)
	      E_Cinetic_y(I)=0.5*AMP(I)*VP(I,2)*VP(I,2)
	      E_Total(I)=E_Potential(I)+E_Cinetic_x(I)+E_Cinetic_y(I)
	      
	      IF (NaN_P2(I).EQ.0) THEN
	          E_TOTAL_FULL = E_TOTAL_FULL+E_Total(I)
	      ENDIF

	      IF (SMOIST_NP(I).LT.0.0) THEN
	        SMOIST_NP(I)=0.0D0
	      ENDIF
	      
!	ACTUALIZA POSICIÓN DE LA PARTÍCULA

	      DO K=1,2
	          XP(I,K)=XP(I,K)+XPP(I,K)
	      ENDDO
	      
	      CALL INVAR2(EPS(I,1),EPS(I,2),EPS(I,4),EPS(I,3)/2.)
	      EPSEQ(I)=2.*Q/3.                                    !Calcula la deformacion total de corte equivalente
	      
	      CALL INVAR2(DEPS(I,1),DEPS(I,2),DEPS(I,4),DEPS(I,3)/2.)
	      EPSEQ2(I)=2.d0*Q/3.d0                                    !Calcula el incremento de deformacion de corte equivalente
	  ENDIF
	ENDDO
	E_TOTAL_FULL = E_TOTAL_FULL/count_active
	
!	WRITE (14,*) TIEMPO,E_TOTAL_FULL

	END
	
!************************************************************************************************************************************
!                                              Subroutine for GID post-process
!************************************************************************************************************************************ 
	
	SUBROUTINE IMPRES_GiD

	INCLUDE 'common_PIV-NP.for'
	CHARACTER*80 ARCH,ARCH_POST_MSH,ARCH_POST_RES
	
	OPEN (2,FILE='PIV-NP.TXT', FORM='FORMATTED')
	READ (2,*) ARCH
	CLOSE(2)
	
	ARCH_POST_MSH=TRIM(ARCH)//'.POST.MSH'
	ARCH_POST_RES=TRIM(ARCH)//'.POST.RES'
	
	IF (IP.EQ.1.AND.IREC.EQ.0) THEN
	  OPEN (11,FILE=ARCH_POST_MSH, STATUS ='UNKNOWN')
	  WRITE (11,100)
	  WRITE (11,200) '"Moved mesh"', 2, 'Point',1
	  WRITE (11,300)
	  WRITE (11,400) 
	  DO I=1,NP
	      WRITE (11,500) I, XP(I,1), XP(I,2)
	  ENDDO
	  WRITE (11,600)  
	  WRITE (11,700)
	  WRITE (11,800)
	  DO I=1,NP
	      WRITE (11,900) I, I, (NaN_P(I)+1)     !Asigna un id de material diferente para partículas NaN en el step 1
        ENDDO
	  WRITE(11,1000)
	  CLOSE(11)
      ELSEIF (IP.EQ.0.AND.IREC.EQ.1) THEN
        OPEN (11,FILE=ARCH_POST_MSH, STATUS ='UNKNOWN')
	  WRITE (11,100)
	  WRITE (11,200) '"Moved mesh"', 2, 'Point',1
	  WRITE (11,300)
	  WRITE (11,400) 
	  DO I=1,NP
	      WRITE (11,500) I, XP(I,1), XP(I,2)
	  ENDDO
	  WRITE (11,600)  
	  WRITE (11,700)
	  WRITE (11,800)
	  DO I=1,NP
	      WRITE (11,900) I, I, (NaN_P(I)+1)     !Asigna un id de material diferente para partículas NaN en el step 1
        ENDDO
	  WRITE(11,1000)
	  CLOSE(11)  
	ENDIF 
   
100   FORMAT('# MESH OF POINTS')  
200   FORMAT('MESH',1x,a15,1x,'dimension',1x,i5,1x,
     .        'ElemType',1x,a15,1x,'Nnode',1x,i5)
300   FORMAT('Coordinates')
400   FORMAT('# PARTICLE_NUMBER',3x,'COORDINATE_X',1x,'COORDINATE_Y') 
500   FORMAT(i14,2e14.6)
600   FORMAT('End Coordinates')
700   FORMAT('Elements')
800   FORMAT('# ELEMENTS',1x, 'MATERIAL') 
900   FORMAT(3i9)        
1000  FORMAT('End Elements')            
	
	
	OPEN (22,FILE=ARCH_POST_RES, STATUS ='UNKNOWN')
	
	IF (IP.EQ.1.AND.IREC.EQ.0) THEN
	    WRITE(22,10)
	ELSEIF (IP.EQ.0.AND.IREC.EQ.1) THEN
	    WRITE(22,10)
	ENDIF
	
	WRITE(22,11) 'Displacement',TIEMPO,'OnNodes'          !Desplazamiento total acumulado
	WRITE(22,12)
	DO I=1,NP
	  IVERSIONCASO = 1
	  CALL UCELDA(I)
	  IF(IDONDE(I).NE.-1) THEN
	  IF (NaN_P(I).EQ.0) THEN
	      WRITE(22,13) I,UP(I,1),UP(I,2)
	  ENDIF
	  ENDIF
	ENDDO
	WRITE(22,14)
	
	WRITE(22,11) 'Inst_displacement',TIEMPO,'OnNodes'     !Desplazamiento instantaneo
	WRITE(22,12)
	DO I=1,NP
	  IVERSIONCASO = 1
	  CALL UCELDA(I)
	  IF(IDONDE(I).NE.-1) THEN
	  IF (NaN_P(I).EQ.0) THEN
	      WRITE(22,13) I,UPO(I,1),UPO(I,2)
	  ENDIF
	  ENDIF
	ENDDO
	WRITE(22,14) 
	
	WRITE(22,18) 'NaNs',TIEMPO,'OnNodes'           !NaNS
	WRITE(22,12)
	DO I=1,NP
	  IVERSIONCASO = 1
	  CALL UCELDA(I)
	  IF(IDONDE(I).NE.-1) THEN
!	  IF (NaN_P(I).EQ.0) THEN
	      WRITE(22,*) I,ICOUNT_NO_NAN(I)
!	  ENDIF
	  ENDIF
	ENDDO
	WRITE(22,14)
	
	WRITE(22,11) 'Velocity',TIEMPO,'OnNodes'              !Velocidad de particulas
	WRITE(22,12)
	DO I=1,NP
	  IVERSIONCASO = 1
	  CALL UCELDA(I)
	  IF(IDONDE(I).NE.-1) THEN
	  IF (NaN_P(I).EQ.0) THEN
	      WRITE(22,13) I,VP(I,1),VP(I,2)
	  ENDIF
	  ENDIF
	ENDDO
	WRITE(22,14)
	
	WRITE(22,11) 'Acceleration',TIEMPO,'OnNodes'          !Aceleración de particulas
	WRITE(22,12)
	DO I=1,NP
	  IVERSIONCASO = 1
	  CALL UCELDA(I)
	  IF(IDONDE(I).NE.-1) THEN
	  IF (NaN_P(I).EQ.0) THEN
	      WRITE(22,13) I,ACEP(I,1),ACEP(I,2)
	  ENDIF
	  ENDIF
	ENDDO
	WRITE(22,14)  
	
	WRITE(22,11) 'Total_strain',TIEMPO,'OnNodes'          !Deformaciones totales acumuladas
	WRITE(22,12)
	DO I=1,NP
	  IVERSIONCASO = 1
	  CALL UCELDA(I)
	  IF(IDONDE(I).NE.-1) THEN
	  IF (NaN_P(I).EQ.0) THEN
	      WRITE(22,15) I,EPS(I,1),EPS(I,2),EPS(I,3)
	  ENDIF
	  ENDIF
	ENDDO
	WRITE(22,14)
	
	WRITE(22,11) 'Inc_strain',TIEMPO,'OnNodes'            !Incremento de deformación total
	WRITE(22,12)
	DO I=1,NP
	  IVERSIONCASO = 1
	  CALL UCELDA(I)
	  IF(IDONDE(I).NE.-1) THEN
	  IF (NaN_P(I).EQ.0) THEN
	      WRITE(22,15) I,DEPS(I,1),DEPS(I,2),DEPS(I,3)
	  ENDIF
	  ENDIF
	ENDDO
	WRITE(22,14)
	
	WRITE(22,18) 'Equi_strain',TIEMPO,'OnNodes'           !Deformacion total de corte equivalente
      WRITE(22,12)
      DO I=1,NP
	  IVERSIONCASO = 1
	  CALL UCELDA(I)
	  IF(IDONDE(I).NE.-1) THEN
	  IF (NaN_P(I).EQ.0) THEN
            WRITE(22,17) I,EPSEQ(I)
	  ENDIF
	  ENDIF
      ENDDO
      WRITE(22,14)
      
      WRITE(22,18) 'Vol_strain',TIEMPO,'OnNodes'           !Deformacion volumetrica total
      WRITE(22,12)
      DO I=1,NP
	  IVERSIONCASO = 1
	  CALL UCELDA(I)
	  IF(IDONDE(I).NE.-1) THEN
	  IF (NaN_P(I).EQ.0) THEN
            WRITE(22,17) I,EVOLUMETRIC(I)
	  ENDIF
	  ENDIF
      ENDDO
      WRITE(22,14)
      
      WRITE(22,18) 'Ins_vol_strain',TIEMPO,'OnNodes'           !Deformacion volumetrica instantanea
      WRITE(22,12)
      DO I=1,NP
	  IVERSIONCASO = 1
	  CALL UCELDA(I)
	  IF(IDONDE(I).NE.-1) THEN
	  IF (NaN_P(I).EQ.0) THEN
            WRITE(22,17) I,EVOL_IN(I)
	  ENDIF
	  ENDIF
      ENDDO
      WRITE(22,14)
      
      WRITE(22,18) 'In_E_strain',TIEMPO,'OnNodes'           !Incremento de Deformacion de corte equivalente
      WRITE(22,12)
      DO I=1,NP
	  IVERSIONCASO = 1
	  CALL UCELDA(I)
	  IF(IDONDE(I).NE.-1) THEN
	  IF (NaN_P(I).EQ.0) THEN
            WRITE(22,17) I,EPSEQ2(I)
	  ENDIF
	  ENDIF
      ENDDO
      WRITE(22,14)
      
      WRITE(22,18) 'E_potential',TIEMPO,'OnNodes'           !Energia potencial
      WRITE(22,12)
      DO I=1,NP
	  IVERSIONCASO = 1
	  CALL UCELDA(I)
	  IF(IDONDE(I).NE.-1) THEN
	  IF (NaN_P(I).EQ.0) THEN
            WRITE(22,17) I,E_Potential(I)
	  ENDIF
	  ENDIF
      ENDDO
      WRITE(22,14)
      
      WRITE(22,18) 'E_kinetic',TIEMPO,'OnNodes'           !Energia cinetica
      WRITE(22,12)
      DO I=1,NP
	  IVERSIONCASO = 1
	  CALL UCELDA(I)
	  IF(IDONDE(I).NE.-1) THEN
	  IF (NaN_P(I).EQ.0) THEN
            WRITE(22,13) I,E_Cinetic_x(I),E_Cinetic_y(I)
	  ENDIF
	  ENDIF
      ENDDO
      WRITE(22,14)
      
      WRITE(22,18) 'E_total',TIEMPO,'OnNodes'           !Energia total
      WRITE(22,12)
      DO I=1,NP
	  IVERSIONCASO = 1
	  CALL UCELDA(I)
	  IF(IDONDE(I).NE.-1) THEN
	  IF (NaN_P(I).EQ.0) THEN
            WRITE(22,17) I,E_Total(I)
	  ENDIF
	  ENDIF
      ENDDO
      WRITE(22,14)
      
      IF (MOISTER.EQ.1) THEN
        WRITE(22,18) 'Moisture',TIEMPO,'OnNodes'           !Humedad
        WRITE(22,12)
        DO I=1,NP
	      IVERSIONCASO = 1
	      CALL UCELDA(I)
	      IF(IDONDE(I).NE.-1) THEN
	      IF (NaN_P(I).EQ.0) THEN
                WRITE(22,17) I,SMOIST_NP(I)
	      ENDIF
	      ENDIF
        ENDDO
        WRITE(22,14)
      
        WRITE(22,18) 'Saturation',TIEMPO,'OnNodes'           !Saturacion
        WRITE(22,12)
        DO I=1,NP
	      IVERSIONCASO = 1
	      CALL UCELDA(I)
	      IF(IDONDE(I).NE.-1) THEN
	      IF (NaN_P(I).EQ.0) THEN
	          WRITE(22,17) I,SATURA_NP(I)
	      ENDIF
	      ENDIF
	  ENDDO
	  WRITE(22,14)
	ENDIF
	            
10    FORMAT('GiD Post Results File 1.0')
11    FORMAT('Result',1x,a15,1x,'Isochrones',1x,e14.6,1x,
     .       'Vector',1x,a15)
12    FORMAT('Values')
13    FORMAT(i14,2e14.6)
14    FORMAT ('End values')
15    FORMAT(i14,3e14.6)
16    FORMAT(i14,4e14.6)
17    FORMAT(i14,e14.6)
18    FORMAT('Result',1x,a15,1x,'Isochrones',1x,e14.6,1x,
     .       'Scalar',1x,a15)
	END
	  
!************************************************************************************************************************************
!                                              Subroutine for localization of numerical particles
!************************************************************************************************************************************ 
	
	SUBROUTINE UCELDA(IPP)
	INCLUDE 'common_PIV-NP.for'
	
	IFI=0
!   BUSQUEDA DE LA CELDA
!	PRIMERO POR FILA

	IF (IVERSION.EQ.1) THEN
	  DO J=1,NFIL
	      IF (XP(IPP,2).GE.YF(J).AND. XP(IPP,2).LT.YF(J+1)) THEN
	          IFI=J
	          GOTO 100
	      ENDIF
	  ENDDO
	ELSE
	  IF (IVERSIONCASO.EQ.1) THEN
	      DO J=1,NFIL2
	          IF (XP(IPP,2).GE.YF2(J).AND. XP(IPP,2).LT.YF2(J+1)) THEN
	              IFI=J
	              GOTO 100
	          ENDIF
	      ENDDO
	  ELSE
	      DO J=1,NFIL2
	          IF (XP2(IPP,2).GE.YF2(J).AND. XP2(IPP,2).LT.YF2(J+1)) THEN
	              IFI=J
	              GOTO 100
	          ENDIF
	      ENDDO
	  ENDIF
	ENDIF

100	IF(IFI.EQ.0) THEN
!	WRITE(15,*)' NO SE ENCUENTRA LA FILA DONDE ESTÁ LA PARTÍCULA ',IPP
!	WRITE(*,*)' NO SE ENCUENTRA LA FILA DONDE ESTÁ LA PARTÍCULA ',IPP
!     * ,' XP= ',XP(IPP,1),XP(IPP,2)
        IDONDE(IPP)=-1
        RETURN

	ENDIF

!	LUEGO LA CELDA DE LA FILA
	INDC=0
	
	IF (IVERSION.EQ.1) THEN
	  YC=YF(IFI)
	  DO IC=NCF(IFI),NCF(IFI+1)-1
	      XC=XF(IFI)+AXC*(IC-NCF(IFI))
	      XC1=XC+AXC
	
! CALCULA LOS NUMEROS DE LOS NODOS INVOLUCRADOS	
	
	      IF (XP(IPP,1).GE.XC.AND.XP(IPP,1).LT.XC1)  THEN
	          INDC=IC
	          IN(1)=NNF(IFI)+IC-NCF(IFI)
	          IN(2)=IN(1)+1
	          IN(3)=NNFA(IFI)+IC-NCF(IFI)
	          IN(4)=IN(3)+1
	          IF(IC.NE.IDONDE(IPP).AND.IP.NE.1) THEN
!	WRITE(15,*) " CAMBIO DE CELDA PARTICULA ",IPP,IC,IDONDE(IPP)
	              IDONDE(IPP)=IC
	          ENDIF
	          RETURN
	      ENDIF
	  ENDDO
	ELSE
	  YC=YF2(IFI)
	  DO IC=NCF2(IFI),NCF2(IFI+1)-1
	      XC=XF2(IFI)+AXC*(IC-NCF2(IFI))
	      XC1=XC+AXC
	
! CALCULA LOS NUMEROS DE LOS NODOS INVOLUCRADOS	
	      IF (IVERSIONCASO.EQ.1) THEN
	          IF (XP(IPP,1).GE.XC.AND.XP(IPP,1).LT.XC1)  THEN
	              INDC=IC
	              IN(1)=NNF2(IFI)+IC-NCF2(IFI)
	              IN(2)=IN(1)+1
	              IN(3)=NNFA2(IFI)+IC-NCF2(IFI)
	              IN(4)=IN(3)+1
	              IF(IC.NE.IDONDE(IPP).AND.IP.NE.1) THEN
!	WRITE(15,*) " CAMBIO DE CELDA PARTICULA ",IPP,IC,IDONDE(IPP)
	                  IDONDE(IPP)=IC
	              ENDIF
	              RETURN
	          ENDIF
	      ELSE
	          IF (XP2(IPP,1).GE.XC.AND.XP2(IPP,1).LT.XC1)  THEN
	              INDC=IC
	              IN(1)=NNF2(IFI)+IC-NCF2(IFI)
	              IN(2)=IN(1)+1
	              IN(3)=NNFA2(IFI)+IC-NCF2(IFI)
	              IN(4)=IN(3)+1
	              IF(IC.NE.IDONDE(IPP).AND.IP.NE.1) THEN
!	WRITE(15,*) " CAMBIO DE CELDA PARTICULA ",IPP,IC,IDONDE(IPP)
	                  IDONDE(IPP)=IC
	              ENDIF
	              RETURN
	          ENDIF
	      ENDIF
	  ENDDO
	ENDIF

	IF (INDC.EQ.0) THEN
!	WRITE(15,*) ' NO SE ENCUENTRA LA CELDA DONDE ESTÁ LA PARTÍCULA',IPP
!	WRITE(*,*) ' NO SE ENCUENTRA LA CELDA DONDE ESTÁ LA PARTÍCULA',IPP
!     * ,' XP= ',XP(IPP,1),XP(IPP,2) 
	  IDONDE(IPP)=-1
	ENDIF
	END
	
!************************************************************************************************************************************
!                                              Subroutine to calculate invariants
!************************************************************************************************************************************
	
	SUBROUTINE INVAR2(SIGX,SIGY,SIGZ,SIGXY)
	INCLUDE 'common_PIV-NP.for'
	REAL*8 SIGX,SIGY,SIGZ,SIGXY,UNO
            
	UNO=1.
	SIGM=(SIGX+SIGY+SIGZ)/3.
	SX=SIGX-SIGM
	SY=SIGY-SIGM
	SZ=SIGZ-SIGM
	SXY=SIGXY

!    	2º INVARIANTE 'RJ2'
	RJ2=(SX*SX+SY*SY+SZ*SZ)/2.+SXY*SXY
	
	IF (RJ2.GT.1.E-10) THEN
	  Q=SQRT(3.*RJ2)                      !DESVIADOR 'Q'
	ELSE                                  !CASO ESPECIAL J2=0
	  Q=0.
	ENDIF
	
      END
	
!************************************************************************************************************************************
!                                              Subroutine to save data for new calculation
!************************************************************************************************************************************
	SUBROUTINE RECOM
	INCLUDE 'common_PIV-NP.for'
	CHARACTER*80 ARCH,ARCH_REC
	
	REWIND 3
	WRITE(3) NP0
	WRITE(3) IVERSION
	WRITE(3) ((XP(I,J),J=1,2),I=1,NP0)
	WRITE(3) ((UP(I,J),J=1,2),I=1,NP0)
	WRITE(3) ((EPS(I,J),J=1,4),I=1,NP0)
	WRITE(3) (EPSEQ(I),I=1,NP0)
	WRITE(3) (NaN_P(I),I=1,NP0)
	
	CLOSE (3)
	END
!************************************************************************************************************************************
!                                                      END PIV-NP
!************************************************************************************************************************************