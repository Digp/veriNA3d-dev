#include <Rcpp.h>
#include <stdio.h>
using namespace Rcpp;

/*
 * Function can now print the whole file in screen and save the first
 * 100 characters, returning them to R.
 * gzip files are incorrectly read!
 */

#define maxchar 100
    
// [[Rcpp::export]]
List cifParserC(std::string strings="")
{
    // Open file
    FILE *file = fopen(strings.c_str(), "r");
    // Check if pointer to file is correct, or stop function
    if (file == NULL)
    {
        printf("Could not open %s\n", strings.c_str());
        return 1;
    }

    // Define object line and index to save the contents of the file
    char line[maxchar];
    int index = 0;
    int newsection = 0;

    // Iterate over the characters of the file and save them
    for (int c = fgetc(file); c != EOF; c = fgetc(file))
    {
        // Detect new mmCIF section based on lines like '# \n'
        newsection = 0;
        if (c == '#') {
            c = fgetc(file);
            if (c == ' ') {
                c = fgetc(file);
                if (c == '\n') {
                    newsection = 1;
                }
            }
        }

        if (newsection) 
        {
            c = fgetc(file);
            //printf("%c", c);
            if (index < maxchar)
            {
                line[index] = ' ';
                index++;
            }
        }
    }

    // Terminate array
    line[index] = '\0';

    // Close text
    fclose(file);

    // Return output
    return List::create(line);
}
