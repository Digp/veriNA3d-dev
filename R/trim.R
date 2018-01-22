#Function directly copied from bio3d. Removes initial and final empty spaces 
#from strings

.trim <- function( s, leading = TRUE, trailing = TRUE ) {
    if ( leading )
        if ( length( grep( "^ +", s ) ) != 0 )
        s <- sub( "^ +", "", s )
    if ( trailing )
        if ( length( grep( " +$", s ) ) != 0 )
        s <- sub( " +$", "", s )
    s[ (s == "") ] <- ""
    s
}

