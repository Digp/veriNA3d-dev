#include <Rcpp.h>
#include <stdio.h>
using namespace Rcpp;

/*
 * Test
 */

// [[Rcpp::export]]
List hello(char str='s')
{

    //char out[2] = {'h', 'f'};
    char out[4] = "asi";
    //out[0] = 'g';
    //out[2] = '\0';
    //out[1] = 'g';
    //char text[] = "dsf";
    //printf("Your string: %s\n", text);

    //char *text = *strings;
    //FILE *file = fopen(text, "r");
    //if (file == NULL)
    //{
    //    printf("Could not open %s.\n", text);
    //    //unload();
    //    return 1;
    //}

    //// Close text
    //fclose(file);

    return List::create(out, str);
}

// List hello(std::vector< std::string > strings)
