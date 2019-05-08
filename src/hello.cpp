#include <Rcpp.h>
#include <stdio.h>
using namespace Rcpp;

/*
 * Test
 */

#define maxchar 10
    
// [[Rcpp::export]]
List hello(std::string strings="")
{
    // Open file
    FILE *file = fopen(strings.c_str(), "r");
    // Check if pointer to file is correct, or stop function
    if (file == NULL)
    {
        printf("Could not open %s.\n", strings.c_str());
        // Close text
        //fclose(file);
        return 1;
    }

    // Define object line and index to save the contents of the file
    char line[maxchar];
    int index = 0;

    // Iterate over the characters of the file and save them
    for (int c = fgetc(file); c != EOF; c = fgetc(file))
    {
        printf("%c\n", c);
        if (index < maxchar)
        {
            line[index] = c;
            index++;
        }
    }

    // Terminate array
    line[index] = '\0';


    // Close text
    fclose(file);

    // Return output
    return List::create(line, strings);
}

// List hello(std::vector< std::string > strings)
    //char out[2] = {'h', 'f'};
//    char out[4] = "asi";
    //out[0] = 'g';
    //out[2] = '\0';
    //out[1] = 'g';
    //char text[] = "dsf";
    //printf("Your string: %s\n", text);
