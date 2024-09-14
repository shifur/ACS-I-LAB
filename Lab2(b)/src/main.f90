PROGRAM lab
  ! Define matrix and system parameters
  INTEGER, PARAMETER :: N = 100         ! Size of the matrix (NxN)
  INTEGER, PARAMETER :: KL = 2          ! Number of subdiagonals
  INTEGER, PARAMETER :: KU = 2          ! Number of superdiagonals
  INTEGER, PARAMETER :: LDAB = 2*KL+KU+1  ! Leading dimension for banded storage
  INTEGER, PARAMETER :: NRHS1 = 5       ! Number of systems to solve in total

  ! Declare variables
  INTEGER :: i, j, k, INFO              ! Loop indices and LAPACK return code
  REAL*8 :: A(LDAB, N), B(N), Bk(N, NRHS1), ERROR(N), l2_norm  ! Matrices and vectors
  REAL*8 :: ExactSolution(N)            ! Vector to store the exact solution
  INTEGER :: IPIV(N)                    ! Pivot indices for LU factorization

  ! External LAPACK subroutines
  EXTERNAL :: DGBTRF, DGBTRS            ! LAPACK routines for LU factorization and solving

  ! Initialize matrix A to zero
  A = 0.0d0

  ! Define matrix A in banded format
  DO i = 1, LDAB
    DO j = 1, N 

      IF (i == 3 .AND. j >= 3) THEN
        A(i,j) = 2.0d0  ! 2nd superdiagonal
      ELSE IF (i == 4 .AND. j >= 2) THEN
        A(i,j) = -1.0d0  ! 1st superdiagonal
      ELSE IF (i == 5) THEN
        A(i,j) = 8.0d0  ! Main diagonal
      ELSE IF (i == 6 .AND. j < N) THEN
        A(i,j) = 1.0d0  ! 1st subdiagonal
      ELSE IF (i == 7 .AND. j < N-1) THEN
        A(i,j) = 3.0d0  ! 2nd subdiagonal
      END IF
    END DO
  END DO

  PRINT *, "Matrix A initialized."

  ! Define the right-hand side vector B
  B = 13.0d0
  B(1) = 9.0d0
  B(2) = 10.0d0
  B(N-1) = 11.0d0
  B(N) = 12.0d0

  PRINT *, "Right-hand side B initialized."

  ! Perform LU factorization of matrix A
  CALL DGBTRF(N, N, KL, KU, A, LDAB, IPIV, INFO)

  IF (INFO /= 0) THEN
    PRINT *, "Matrix factorization failed with INFO =", INFO
    STOP
  END IF

  PRINT *, "Matrix factorization successful."

  ! Initialize Bk(:,1) as the first right-hand side
  Bk(:,1) = B

  ! Set the exact solution vector to 1's
  ExactSolution = 1.0d0

  ! Loop to solve multiple systems iteratively
  DO k = 1, NRHS1
    ! Solve the system A * X = Bk(:,k)
    CALL DGBTRS('N', N, KL, KU, 1, A, LDAB, IPIV, Bk(:,k), N, INFO)
    
    IF (INFO /= 0) THEN
      PRINT *, "System solve failed for system", k, "with INFO =", INFO
      STOP
    END IF
    
    ! Compute L2 norm of the error for this system
    ERROR = Bk(:,k) - ExactSolution  ! Error = Computed solution - Exact solution
    l2_norm = 0.0d0
    DO i = 1, N
      l2_norm = l2_norm + ERROR(i)**2
    END DO
    l2_norm = SQRT(l2_norm)

    ! Conditionally format the output based on the value of l2_norm
    IF (l2_norm < 1.0d-3) THEN
      PRINT '(A,I1,A,ES25.16)', "L2 norm of the error for system ", k, ": ", l2_norm
    ELSE
      PRINT '(A,I1,A,F25.16)', "L2 norm of the error for system ", k, ": ", l2_norm
    END IF

    IF (k < NRHS1) THEN
      ! For the next system, set Bk(:, k+1) = B + current solution Bk(:,k)
      Bk(:, k+1) = B + Bk(:,k)
    END IF
  END DO

END PROGRAM lab




