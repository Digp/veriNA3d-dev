
#Function to shift 360 degrees torsion angles
shift360 <- 
function(tor) {
  if( !is.null( tor ) && !is.na( tor ) && length( tor ) > 0){

    if( tor < 0 ) {
      tor_shifted <- tor + 360
    } else {
      tor_shifted <- tor
    }

    return( tor_shifted )

  } else {
    return( NA )
  }
}


