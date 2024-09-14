!     .. Parameters ..
      INTEGER          N, NRHS
      PARAMETER        ( N = 100, NRHS = 1 ) ! Define matrix size N and number of right-hand sides NRHS
      INTEGER          LDA, LDB
      PARAMETER        ( LDA = N, LDB = N )  ! Leading dimensions of A and B arrays

!     .. Local Scalars ..
      INTEGER          INFO  ! Variable to store info from DGESV call (error status)

!     .. Local Arrays ..
      INTEGER          IPIV( N )  ! Array to store pivot indices
      DOUBLE PRECISION A( LDA, N ), B( LDB, NRHS ), X( N, NRHS ) ! Arrays for matrix A, vectors B and X

!     External Subroutines
      EXTERNAL         DGESV  ! External LAPACK routine to solve linear systems
      EXTERNAL         PRINT_MATRIX, PRINT_INT_VECTOR  ! External subroutines for printing

!     Initialize matrix A and vector X
      DO i = 1, N
            DO j = 1, N
                  IF ((i == j) .OR. (j == N)) THEN
                        A(i, j) = 1.0D0   ! Set diagonal elements and last column to 1
                  ELSE IF (j > i) THEN
                        A(i, j) = -10.0D0  ! Set elements above the diagonal to -10
                  ELSE
                        A(i, j) = 0.0D0   ! Set remaining elements to 0
                  END IF
            END DO
            X(i, NRHS) = DBLE(i) / N  ! Initialize X as i/n
      END DO

!     Compute B by multiplying A with X
      B = MATMUL(A, X)

!     .. Executable Statements ..
      WRITE(*,*) 'DGESV Example Program Results'  ! Display header for results

!     Solve the equations A*X = B
      CALL DGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)  ! Call LAPACK DGESV to solve the system

!     Check for the exact singularity
      IF (INFO .GT. 0) THEN
         WRITE(*,*) 'The diagonal element of the triangular factor of A,'   ! Error message if matrix is singular
         WRITE(*,*) 'U(', INFO, ',', INFO, ') is zero, so that'
         WRITE(*,*) 'A is singular; the solution could not be computed.'
         STOP
      END IF

!     Print the solution vector B
      CALL PRINT_MATRIX('Solution', N, 1, B, LDB)

!     Optionally print details of LU factorization (commented out)
      !CALL PRINT_MATRIX('LU Factorization', N, N, A, LDA)

!     Optionally print pivot indices (commented out)
      !CALL PRINT_INT_VECTOR('Pivot indices', N, IPIV)

      STOP
      END

!     End of DGESV Example

!  =============================================================================

!     Auxiliary routine: printing a matrix
      SUBROUTINE PRINT_MATRIX(DESC, M, N, A, LDA)
      CHARACTER*(*)    DESC  ! Description of the matrix
      INTEGER          M, N, LDA  ! Matrix dimensions
      DOUBLE PRECISION A(LDA, *)  ! Matrix to be printed

      INTEGER          I, J  ! Loop indices

      WRITE(*,*)  ! Blank line
      WRITE(*,*) DESC  ! Print description
      DO I = 1, M
          WRITE(*, 9998) (A(I, J), J = 1, N)  ! Print matrix row-wise
      END DO

!   Print matrix elements in scientific notation
 9998 FORMAT(11(1X, E12.5))  ! Formatting for scientific notation (11 elements per row, width 12, 5 decimal places)
      RETURN
      END

!     Auxiliary routine: printing a vector of integers
      SUBROUTINE PRINT_INT_VECTOR(DESC, N, A)
      CHARACTER*(*)    DESC  ! Description of the vector
      INTEGER          N  ! Length of the vector
      INTEGER          A(N)  ! Vector to be printed

      INTEGER          I  ! Loop index

      WRITE(*,*)  ! Blank line
      WRITE(*,*) DESC  ! Print description
      WRITE(*, 9999) (A(I), I = 1, N)  ! Print vector elements

!   Print integer elements (unchanged)
 9999 FORMAT(11(1X, I6))  ! Integer vector formatting (up to 11 integers per row)
      RETURN
      END
