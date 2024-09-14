PROGRAM lab
  ! Define parameters for the matrix and system dimensions
  INTEGER, PARAMETER :: N = 100         ! Size of the matrix (NxN)
  INTEGER, PARAMETER :: KL = 2          ! Number of subdiagonals
  INTEGER, PARAMETER :: KU = 2          ! Number of superdiagonals
  INTEGER, PARAMETER :: LDAB = 2*KL+KU+1  ! Leading dimension for banded storage
  INTEGER, PARAMETER :: NRHS = 1        ! Number of right-hand sides (b vector)

  ! Declare variables
  INTEGER :: i, j, INFO                 ! Loop variables and info for LAPACK
  REAL*8 :: A(LDAB,N), B(N,NRHS), ERROR(N), l2_NORM
  INTEGER :: IPIV(N)                    ! Pivot indices for factorization

  ! External LAPACK subroutines
  EXTERNAL :: DGBTRF, DGBTRS            ! LAPACK routines for LU factorization and solving

  ! Initialize matrix A to zero
  A = 0.0d0

  ! Define matrix A in banded format
  DO j = 1, N
    IF (j >= 3) A(3,j) = 2.0d0             ! 2nd superdiagonal (2 places above the main diagonal)
    IF (j >= 2) A(4,j) = -1.0d0            ! 1st superdiagonal (1 place above the main diagonal)
    A(5,j) = 8.0d0                         ! Main diagonal
    IF (j < N) A(6,j) = 1.0d0              ! 1st subdiagonal (1 place below the main diagonal)
    IF (j < N-1) A(7,j) = 3.0d0            ! 2nd subdiagonal (2 places below the main diagonal)
  END DO

  ! Print the matrix A for verification
  PRINT *, "Matrix A:"
  DO i = 1, LDAB
    WRITE(*,'(1X,*(F8.3))') (A(i,j), j=1, N)
  END DO

  ! Define the right-hand side vector B
  B = 13.0d0                              ! Set all elements of B to 13.0
  B(1,1) = 9.0d0                          ! Special case for B(1,1)
  B(2,1) = 10.0d0                         ! Special case for B(2,1)
  B(N-1,1) = 11.0d0                       ! Special case for B(N-1,1)
  B(N,1) = 12.0d0                         ! Special case for B(N,1)

  ! Print the vector B for verification
  PRINT *, "Matrix B:"
  DO i = 1, N
    WRITE(*,'(1X,F8.3)') B(i,1)
  END DO

  ! Compute the LU factorization of the banded matrix A
  CALL DGBTRF(N, N, KL, KU, A, LDAB, IPIV, INFO)

  IF (INFO == 0) THEN
    PRINT *, "Matrix factorization successful"
  ELSE
    PRINT *, "Matrix factorization failed with INFO = ", INFO
    STOP
  END IF

  ! Solve the system A * X = B
  CALL DGBTRS('N', N, KL, KU, NRHS, A, LDAB, IPIV, B, N, INFO)

  IF (INFO == 0) THEN
    PRINT *, "Matrix solve successful"
  ELSE
    PRINT *, "Matrix solve failed with INFO = ", INFO
    STOP
  END IF

  ! Print the solution vector X (stored in B after solve)
  PRINT *, "Solution X:"
  DO i = 1, N
    WRITE(*,'(1X,F8.3)') B(i,1)
  END DO

  ! Compute the error vector (computed solution - exact solution)
  ERROR = B(1:N,1) - 1.0d0                ! Exact solution is 1.0 for all elements

  ! Print the error vector for verification
  PRINT *, "Error vector:"
  DO i = 1, N
    WRITE(*,'(1X,F8.3)') ERROR(i)
  END DO

  ! Compute the l2 norm of the error
  l2_NORM = 0.0d0
  DO i = 1, N
    l2_NORM = l2_NORM + ERROR(i)**2
  END DO
  l2_NORM = SQRT(l2_NORM)

  PRINT *, "The l2 norm of the errors:", l2_NORM

END PROGRAM lab
