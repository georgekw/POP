PROGRAM Structure
  USE FEMModule
  IMPLICIT NONE

  REAL, ALLOCATABLE, DIMENSION (:,:) :: sk !The stiffness matrix of the system.
  REAL, ALLOCATABLE, DIMENSION(:) :: loads !A vector of applied loads.
  
  CALL READFILE()
  CALL MAIN()

  CONTAINS
    
    SUBROUTINE CONSTRUCTSYSTEMSTIFFNESSMATRIX (element)
      IMPLICIT NONE
      REAL, DIMENSION(7) :: element
      REAL, DIMENSION(4,2) :: indices
      REAL, DIMENSION(6,6) :: k
      INTEGER :: i, j, u, v
      
      k = CONSTRUCTELEMENTSTIFFNESSMATRIX (element(1), element(2), element(3), element(4), element(5))
      u = 3*(element(6)-1) !Index increase
      v = 3*(element(7)-2) !Index increase

      DO i=1,6
         DO j=1,6
            IF(i <= 3 .AND. j <= 3) THEN
               sk(i+u,j+u) = sk(i+u,j+u) + k(i,j)
            ELSEIF(i <= 3 .AND. j >= 4) THEN
               sk(i+u,j+v) = sk(i+u,j+v) + k(i,j)
            ELSEIF(i >= 4 .AND. j <= 3) THEN
               sk(i+v, j+u) = sk(i+v, j+u) + k(i,j)
            ELSEIF(i >=4 .AND. j >= 4) THEN
               sk(i+v,j+v) = sk(i+v,j+v) + k(i,j)
            END iF
         END DO
      END DO
    END SUBROUTINE CONSTRUCTSYSTEMSTIFFNESSMATRIX

    SUBROUTINE READFILE()
      IMPLICIT NONE
      INTEGER :: i, elements, joints
      REAL, ALLOCATABLE, DIMENSION(:,:) :: elems !A matrix of the elements.
      REAL, ALLOCATABLE, DIMENSION(:) :: constraints !A vector of fixed constraints.

      OPEN (5, file='input1.dat')
      READ(5,*)
      READ(5,*) elements
      READ(5,*)
      READ(5,*) joints
      READ(5,*)
      
      ALLOCATE(elems(elements,7))
      ALLOCATE(sk(3*joints,3*joints)) !3 degress of freedom for each nodes.
      ALLOCATE(loads(3*joints))
      ALLOCATE(constraints(3*joints))
      sk = 0.0

      DO i=1, elements
         READ (5,*) elems(i,1:7)
         CALL CONSTRUCTSYSTEMSTIFFNESSMATRIX(elems(i,1:7))
      END DO
      
      READ (5,*)
      READ (5,*) loads
      READ (5,*)
      READ (5,*) constraints
      CLOSE(5)
      CALL CONSTRAINTSATISFACTION(sk, constraints, loads)
    END SUBROUTINE READFILE

    SUBROUTINE WRITEFILE(result)
      IMPLICIT NONE
      INTEGER :: elements, i
      REAL, DIMENSION(:) :: result
      REAL, DIMENSION(6) :: nodalForces !The member-end forces of an element.
      REAL, ALLOCATABLE, DIMENSION(:,:) :: elems !A matrix of the elements. 

      OPEN (5, file='input1.dat')
      READ(5,*)
      READ(5,*) elements
      READ(5,*)
      READ(5,*)
      READ(5,*)
      ALLOCATE(elems(elements,7))

      OPEN (6, file='result.dat')
      WRITE(6, *) 'RESULT FILE:'
      WRITE(6, *) 'System displacements:'
      WRITE(6, *) result
      WRITE(6, *) 'Nodal forces of each element:'

      DO i=1, elements
         READ (5,*) elems(i,1:7)
         nodalForces = CONSTRUCTMEMBERENDFORCES(elems(i,1:7), loads) !loads is now representing system deflections.
         WRITE(6, *) nodalForces
      END DO
      
      CLOSE(5)
      CLOSE(6)
    END SUBROUTINE WRITEFILE

    SUBROUTINE MAIN()
      IMPLICIT NONE
      INTEGER :: skSize
      REAL, ALLOCATABLE, DIMENSION(:,:) :: augM !An augmented matrix of the stiffness and load matrices.
      ALLOCATE(augM(SIZE(sk),SIZE(sk)+1))
      augM = 0.0
      augM = CONSTRUCTAUGMENTEDMATRIX(sk, loads)
      loads = GAUSSELIMINATION(augM) !Overwrite deflection values to load vector.
      CALL WRITEFILE(loads)
    END SUBROUTINE MAIN

END PROGRAM Structure
