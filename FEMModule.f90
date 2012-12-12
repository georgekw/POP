MODULE FEMModule

  !FEM utility procedures.

  CONTAINS

    FUNCTION CONSTRUCTELEMENTSTIFFNESSMATRIX (E, I, A, L, R) RESULT(k)
      IMPLICIT NONE
      REAL :: E,I,A,L,R
      REAL, DIMENSION(6,6) :: k !The stiffness matrix of an element.
      REAL, DIMENSION(6,6) :: C !The transformation matrix.
      REAL, DIMENSION(6,6) :: tC !The transposed transformation matrix.
      REAL, DIMENSION(6,6) :: kl !The leveled stiffness matrix.

      C = CONSTRUCTELEMENTTRANSFORMATIONSTIFFNESSMATRIX(R)
      tC = TRANSPOSE(C)
      kl = CONSTRUCTELEMENTLEVELEDSTIFFNESSMATRIX(E,I,A,L)
      k = MATMUL(MATMUL(C,kl),tC) !C*kl*tC
    END FUNCTION CONSTRUCTELEMENTSTIFFNESSMATRIX

    FUNCTION CONSTRUCTELEMENTLEVELEDSTIFFNESSMATRIX (E, I, A, L) RESULT(ek)
      IMPLICIT NONE
      REAL :: E,I,A,L
      REAL, DIMENSION (6,6) :: ek !The stiffness matrix of an element.
      REAL :: EAL !Stretching/Compression stiffness.
      REAL :: dEIL !Beam-end displacement stiffness.
      REAL :: rEIL !Beam-end rotation stiffness.
      REAL :: bEIL !Beam-end bending stiffness.

      EAL = E*A/L
      dEIL = 12*E*I/L**3
      rEIL = 6*E*I/L**2
      bEIL = E*I/L
      ek = 0.0
      ek(1,1) = EAL
      ek(4,1) = -EAL
      ek(2,2) = dEIL
      ek(3,2) = -rEIL
      ek(5,2) = -dEIL
      ek(6,2) = -rEIL
      ek(2,3) = -rEIL
      ek(3,3) = 4*bEIL
      ek(5,3) = rEIL
      ek(6,3) = 2*bEIL
      ek(1,4) = -EAL
      ek(4,4) = EAL
      ek(2,5) = -dEIL
      ek(3,5) = rEIL
      ek(5,5) = dEIL
      ek(6,5) = rEIL
      ek(2,6) = -rEIL
      ek(3,6) = 2*bEIL
      ek(5,6) = rEIL
      ek(6,6) = 4*bEIL
    END FUNCTION CONSTRUCTELEMENTLEVELEDSTIFFNESSMATRIX

    FUNCTION CONSTRUCTELEMENTTRANSFORMATIONSTIFFNESSMATRIX(theta) RESULT(C)
      IMPLICIT NONE
      REAL :: theta
      REAL, DIMENSION(6,6) :: C !The transformation matrix

      C = 0.0
      C(1,1) = cos(theta)
      C(1,2) = -sin(theta)
      C(2,1) = sin(theta)
      C(2,2) = cos(theta)
      C(3,3) = 1
      C(4,4) = cos(theta)
      C(4,5) = -sin(theta)
      C(5,4) = sin(theta)
      C(5,5) = cos(theta)
      C(6,6) = 1
    END FUNCTION CONSTRUCTELEMENTTRANSFORMATIONSTIFFNESSMATRIX

    FUNCTION CONSTRUCTAUGMENTEDMATRIX(matrixA, vectorB) RESULT(AugM)
      IMPLICIT NONE
      REAL, DIMENSION(:,:) :: matrixA
      REAL, DIMENSION(:) :: vectorB
      REAL, ALLOCATABLE, DIMENSION(:,:) :: AugM
      INTEGER :: aSize

      aSize = SIZE(matrixA,1)
      ALLOCATE(AugM(aSize,aSize+1))
      
      AugM(1:aSize,1:aSize) = matrixA
      AugM(1:aSize,aSize+1) = vectorB
      
    END FUNCTION CONSTRUCTAUGMENTEDMATRIX

    FUNCTION GAUSSElIMINATION(Aug) RESULT(X)
      IMPLICIT NONE
      REAL, DIMENSION (:,:) :: Aug
      REAL, DIMENSION (SIZE(Aug,1)):: X
      REAL :: temp, m, sums
      INTEGER :: n, i, j, k, p, max
      LOGICAL :: isSolution

      isSolution = .FALSE.
      n = SIZE(X)
      X = 0.0
      DO k = 1, n-1
         DO j = 1, n
            IF(Aug(j, 1) /= 0) isSolution = .TRUE.
         END DO
         IF(isSolution .EQV. .FALSE.) THEN
            PRINT *, 'No unique solution exists: A is singular'
            EXIT
         ELSE
            max = k
            DO j = k, n
            IF(ABS(Aug(j, k)) > ABS(Aug(max,k))) max = j
            END DO                       
            DO j = 1, n + 1
               temp = Aug(k, j)
               Aug(k,j) = Aug(max,j)
               Aug(max,j)=temp
            END DO
         ENDIF
         DO j = k+1, n
            m = Aug(j,k)/Aug(k,k)
            DO p = k+1, n+1
               Aug(j,p) = Aug(j,p) - m*Aug(k,p)
            END DO
         END DO
         IF(Aug(n,n) == 0) THEN
            PRINT *, 'No unique solution exists: Infinitely many solutions'
            EXIT
         ELSE
            X(n) = Aug(n,n+1)/Aug(n,n)
            DO i = n-1, 1, -1
               sums = 0
               DO j = i+1, n
                  sums = sums + Aug(i,j)*X(j)
               END DO
               X(i) = (Aug(i,n+1) - sums)/Aug(i,i)
            END DO
         ENDIF
      END DO
    END FUNCTION GAUSSELIMINATION

    SUBROUTINE CONSTRAINTSATISFACTION(k, constraintVector, loads)
      IMPLICIT NONE
      REAL, DIMENSION(:,:) :: k !A stiffness matrix.
      REAL, DIMENSION(:) :: constraintVector !A vector of the constraints.
      REAL, DIMENSION(:) :: loads !A vector of the loads.
      INTEGER :: vectorSize, r

      vectorSize = SIZE(constraintVector)
      DO r=1,vectorSize
         IF(constraintVector(r) == 0) THEN
            k(r,:) = 0
            k(:,r) = 0
            k(r,r) = 1
            loads(r) = 0
         END IF
      END DO
    END SUBROUTINE CONSTRAINTSATISFACTION

    FUNCTION CONSTRUCTELEMENTDISPLACEMENTVECTOR(nodeOne, nodeTwo, sysDisVec) RESULT(r)
      IMPLICIT NONE
      REAL, DIMENSION(6) :: r !The displacement vector
      REAL, DIMENSION(:) :: sysDisVec !The displacement vector of the system.
      INTEGER :: i, u, v
      REAL :: nodeOne, nodeTwo
      
      u = 3*(nodeOne-1) !Index increase
      v = 3*(nodeTwo-2) !Index increase
      DO i=1,6
         IF(i <= 3) THEN
            r(i) = sysDisVec(i+u)
         ELSE
            r(i) = sysDisVec(i+v)
         ENDIF
      END DO
    END FUNCTION CONSTRUCTELEMENTDISPLACEMENTVECTOR

    FUNCTION CONSTRUCTMEMBERENDFORCES(elems, sysDisVec) RESULT(S)
      IMPLICIT NONE
      REAL, DIMENSION(:) :: sysDisVec
      REAL, DIMENSION(7) :: elems 
      REAL, DIMENSION(6) :: v, S !S represents member-end forces.
      
      v = CONSTRUCTELEMENTDISPLACEMENTVECTOR(elems(6), elems(7), sysDisVec)
      S = MATMUL(CONSTRUCTELEMENTSTIFFNESSMATRIX (elems(1), elems(2), elems(3), elems(4), elems(5)),v)
    END FUNCTION CONSTRUCTMEMBERENDFORCES

    SUBROUTINE PRINTMATRIX(M)
      REAL M(:,:)
      INTEGER rows, cols

      rows = SIZE(M,1)
      cols = SIZE(M,2)
      DO i=1,rows
         WRITE(*,*), (M(i,j), j=1, cols) 
      END DO
    END SUBROUTINE PRINTMATRIX


END MODULE FEMModule
            
