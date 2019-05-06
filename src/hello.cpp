#include <Rcpp.h>
using namespace Rcpp;

/*
 * Test
 */

// [[Rcpp::export]]
List hello()
{

    char out = 'h';
    return List::create(out);
}
