
measure_pucker <-
function( nu0, nu1, nu2, nu3, nu4 ){
    #Ribose pucker
    pu_vec <- c( nu2, nu3, nu4, nu0, nu1)
    #Compute pucker
        sumA <- 0
        sumB <- 0
        for( rt in 1:5 ) {
              sumA <- sumA + ( pu_vec[ rt ] * cos( (4/5) * pi * ( rt-1 ) ))
              sumB <- sumB + ( pu_vec[ rt ] * sin( (4/5) * pi * ( rt-1 ) ))
        }
        A <- (2/5) * sumA
        B <- -(2/5) * sumB
        pu_amp <- round( sqrt( (A^2) + (B^2) ) ,3)
        pu_phase <- round( atan2( B, A ) * (180/pi), 3)

    #Shift 360
    pu_amp <- shift360( pu_amp )
    pu_phase <- shift360( pu_phase )
    return( list( pu_phase = pu_phase, pu_amp = pu_amp ))
}
