rm ( list = ls() ) 

## THIS MODULE CONTAINS TWO FUNCTIONS THAT CARRY OUT MATRIX INVERSION.
##
## THE makeCacheMatrix() FUNCTION CREATES ROUTINES FOR STORING AND
## RETRIEVING MATRICES AND THEIR INVERSES.
##
## THE cacheSolve() ROUTINE INVERTS MATRICES USING THE BLOCKWISE
## INVERSION FORMULA (EQUATION 1) PROVIDED ON THE WIKIPEDIA 
## "INVERTIBLE MATRIX" PAGE. 
##
## AT THE BOTTOM OF THE FILE ARE SEVERAL "TEST" CALLS TO cacheSolve()
## DEMONSTRATING THAT THE SOFTWARE WORKS. 

makeCacheMatrix <- function ( m = matrix() ) {
    
    # The internal stored inverse is initialized to NULL.
    i <- NULL
    
    setmat <- function ( mat ) {
        # Set the internal stored matrix to the argument.
        # Reinitialize the outdated stored inverse to NULL.
        
        n <- dim ( mat )
        
        m <<- matrix ( NA_real_, nrow=n, ncol=n )
        m <<- mat 
    }
    
    getmat <- function () {
        # Retrieve the value of the stored matrix.
        
        m
    }
    
    setinv <- function ( inv ) {
        # Set the internal stored inverse to the argument.
        
        n <- dim ( inv )
        i <- matrix ( NA_real_, nrow=n, ncol=n )
        i <<- inv
    }
    
    getinv <- function () {
        # Retrieve the value of the stored inverse.
        
        i
    }
    
    list ( setmat = setmat, 
           getmat = getmat,
           setinv = setinv,
           getinv = getinv )
    
}

## THIS FUNCTION cacheSolve() EMPLOYS A RECURSION-LIKE FORMAT 
## TO CONSTRUCT THE INVERSE OF A NXN MATRIX X. IT WORKS AS FOLLOWS:
##
## (1) GENERATE THE INVERSE OF X2, THE TOP LEFT 2x2 SUBMATRIX OF X
## (2) USING THE BLOCK INVERSION FORMULA, INVERT X3, THE TOP LEFT
##     3x3 MATRIX, FROM THE OUTPUT OF STEP 1.
## (3) USING THE BLOCK INVERSION FORMULA, INVERT X4, THE TOP LEFT
##     4x4 MATRIX, FROM THE OUTPUT OF STEP 2.
## 
## THIS PROCESS CONTINUES UNTIL THE INVERSE OF X HAS BEEN CREATED.
## THE INVERSE FROM STEP (N-1) IS STORED, AND RETRIEVED WHEN NEEDED
## IN STEP (N).
## 
## TEST CALLS ARE INCLUDED AT THE BOTTOM OF THE FILE.

cacheSolve <- function ( x, ... ) {
    
    CacheManager <- makeCacheMatrix ()
    
    # FIRST CALL TO SETMAT()
    CacheManager$setmat ( x[1:2,1:2] )
    
    INDEX <- 1
    
    DONE <- FALSE
    
    while ( DONE == FALSE )  { 
        
        if ( nrow ( CacheManager$getmat() ) == 2 ) {
            
            CacheManager$setinv ( inv2 ( x[1:2,1:2] ) )
            
            if ( nrow ( x ) == 2 ) {
                 DONE <- TRUE
            } else  {
                 CacheManager$setmat ( x[1:3,1:3] )
            }
            
        }   else {
            
            A1  <- CacheManager$getinv()
            n   <- nrow ( A1 )
            
            TESTMAT3 <- CacheManager$getmat()
            
            B   <- matrix ( x[1:n,n+1], nrow=n, ncol=1 )
            C   <- matrix ( x[n+1,1:n], nrow=1, ncol=n )
            D   <- matrix ( x[n+1,n+1], nrow=1, ncol=1 )
            
            inv <- matrix ( rep(0,(n+1)*(n+1)), nrow=n+1, ncol=n+1 )
            
            E1  <- 1 / ( D - C%*%A1%*%B )
            
            inv[1:n,1:n] <-  A1 + A1 %*% B %*% E1 %*% C %*% A1
            inv[1:n,n+1] <- -A1 %*% B %*% E1
            inv[n+1,1:n] <- -E1 %*% C %*% A1
            inv[n+1,n+1] <-  E1
            
            CacheManager$setinv ( inv )

            TEST1 <- nrow ( CacheManager$getinv() )
            TEST2 <- nrow ( x )
            
            if ( TEST1 == TEST2 ) {
                 DONE <- TRUE
            } else  {
                 TEST2 <- nrow ( CacheManager$getinv() )
                 CacheManager$setmat ( x[1:TEST1+1,1:TEST1+1] )                
            }
            
            INDEX <- INDEX + 1
            
        }      
        
    } # END WHILE
    
    inv <- CacheManager$getinv()
}

inv2 <- function ( x ) {
    
        a     <- x[1,1]
        b     <- x[1,2]
        c     <- x[2,1]
        d     <- x[2,2]
        
        coeff <- 1/(a*d-b*c)
        
        inv <- coeff * matrix ( c(d,-c,-b,a), nrow=2, ncol=2 )
}


# INVERT A 2x2 MATRIX
A2 <- matrix ( c(1,0,-1,3), nrow=2, ncol=2 )
T2 <- cacheSolve ( A2 )
T2 %*% A2

# INVERT A 3x3 MATRIX
A3 <- matrix ( c(0,-1,4,3,2,5,-2,4,1), nrow=3, ncol=3 )
T3 <- cacheSolve ( A3 )
T3 %*% A3

# INVERT A 4x4 MATRIX
A4 <- matrix ( c(0,-2,1,4,3,3,6,-7,2,4,1,0,10,-3,-6,4), nrow=4, ncol=4 )
T4 <- cacheSolve ( A4 )
T4 %*% A4

# INVERT A 5x5 MATRIX 
A5 <- matrix ( c(1,2,3,4,5,4,7,-1,0,-5,6,3,1,-2,4,4,4,1,-2,-2,3,0,6,2,0), nrow=5, ncol=5 )
T5 <- cacheSolve ( A5 )
T5 %*% A5

# INVERT A 10x10 MATRIX
A10 <- matrix ( c(1,2,3,4,5,6,7,8,9,9,
                  6,7,4,2,1,0,5,7,9,6,
                  0,0,6,1,7,7,9,2,4,5,
                  9,5,5,7,2,7,8,8,6,4,
                  1,6,5,4,8,9,7,6,5,0,
                  3,2,7,5,4,9,7,8,5,4,
                  0,6,5,4,7,9,2,5,6,7,
                  1,3,0,9,6,8,4,5,2,7,
                  8,6,9,2,5,3,7,4,2,3,
                  1,7,4,8,2,9,6,8,7,4), nrow=10, ncol=10 )
T10 <- cacheSolve ( A10 )
T10 %*% A10

