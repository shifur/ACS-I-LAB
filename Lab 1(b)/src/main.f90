      PROGRAM MAIN
      INTEGER          N, NRHS
      PARAMETER        ( N = 2, NRHS = 1 )  ! N = number of equations, NRHS = number of right-hand sides
      
      INTEGER          LDA, LDB
      PARAMETER        ( LDA = N, LDB = N )  ! LDA, LDB define leading dimensions of matrices A and B
      
      ! Local Scalars
      INTEGER          INFO  ! INFO will store the status of the solution from DGESV
      
      ! Local Arrays
      INTEGER          IPIV( N )  ! Pivot indices for the LU factorization
      REAL(8) A( LDA, N ), B( LDB, NRHS ), C( LDB, NRHS )  ! Matrices A, B, and C

      ! External Subroutines
      EXTERNAL         DGESV  ! LAPACK routine to solve the system
      EXTERNAL         PRINT_MATRIX, PRINT_INT_VECTOR  ! Subroutines to print matrices and vectors

      ! Initialize the matrices A, B, and C with values
      A = RESHAPE([4.5, 1.6, 3.1, 1.1],[N,N])  ! Matrix A initialized as a 2x2 matrix
      B = RESHAPE([19.249, 6.843],[N,NRHS])    ! Matrix B initialized as a 2x1 column matrix
      C = RESHAPE([19.25, 6.84],[N,NRHS])      ! Matrix C initialized similarly to B

      ! Print the initial matrices A, B, and C
      CALL PRINT_MATRIX( "Matrix A", N, N, A, LDA )
      CALL PRINT_MATRIX( "Matrix B", N, NRHS, B, LDB )
      CALL PRINT_MATRIX( "Matrix C", N, NRHS, C, LDB )

      ! Solve the equations A*X = B using DGESV
      CALL DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      CALL PRINT_MATRIX( 'Solution B', N, 1, B, LDB )  ! Print the solution stored in B
      
      IF(INFO == 0) THEN
          WRITE(*, *) "Task implemented successfully for Ax = B!" 
      ELSE
          PRINT *, "INFO=", INFO
      END IF

      ! Reset matrix A and solve A*X = C
      A = RESHAPE([4.5, 1.6, 3.1, 1.1],[N,N])
      CALL DGESV( N, NRHS, A, LDA, IPIV, C, LDB, INFO )
      CALL PRINT_MATRIX( 'Solution C', N, 1, C, LDB )  ! Print the solution stored in C
      
      IF(INFO == 0) THEN
          WRITE(*, *) "Task implemented successfully for Ax = C!" 
      ELSE
          PRINT *, "INFO=", INFO
      END IF

      END PROGRAM MAIN

     ! Subroutine to print a matrix
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      REAL(8) A( LDA, * )
      INTEGER          I, J

      WRITE(*,*) ! Newline for better output formatting
      WRITE(*,*) DESC  ! Print the matrix description

      ! Loop over rows and print each row of the matrix
      DO I = 1, M
          WRITE(*,*) ( A( I, J ), J = 1, N )  ! Print each element in the row
      END DO

      END SUBROUTINE PRINT_MATRIX

     ! Subroutine to print an integer vector (used for pivot indices)
      SUBROUTINE PRINT_INT_VECTOR( DESC, N, A )
      CHARACTER*(*)    DESC
      INTEGER          N
      INTEGER          A( N )
      INTEGER          I

      WRITE(*,*)  ! Newline for better output formatting
      WRITE(*,*) DESC  ! Print the vector description
      
      ! Print the vector elements
      WRITE(*,9999) ( A( I ), I = 1, N )
      
 9999 FORMAT( 11(:,1X,I6) )  ! Format for printing integers
      RETURN
      END SUBROUTINE PRINT_INT_VECTOR

