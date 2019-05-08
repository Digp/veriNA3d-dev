#include <Rcpp.h>
#include <stdio.h>
using namespace Rcpp;

/*
 * Test
 */

// [[Rcpp::export]]
List hello(char str='s', std::string strings="")
{

    //char out[2] = {'h', 'f'};
    char out[4] = "asi";
    //out[0] = 'g';
    //out[2] = '\0';
    //out[1] = 'g';
    //char text[] = "dsf";
    //printf("Your string: %s\n", text);

    //char *text = *strings;
    FILE *file = fopen(strings.c_str(), "r");
    if (file == NULL)
    {
        printf("Could not open %s.\n", strings.c_str());
        // Close text
        fclose(file);
        return 1;
    }

    // Close text
    fclose(file);

    return List::create(out, str, strings);
}

// List hello(std::vector< std::string > strings)
