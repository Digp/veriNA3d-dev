#include <Rcpp.h>
#include <stdio.h>
using namespace Rcpp;

/*
 * Function can now detect new sections and return to R the first character
 * of each section.
 * gzip files are incorrectly read!
 */

// Temporal macro definition
#define maxchar 100

// Helper function to detect new mmCIF section based on lines like '# \n'
int newsec(FILE *file, int c)
{
    if (c == '#') {
        c = fgetc(file);
        if (c == ' ') {
            c = fgetc(file);
            if (c == '\n') {
                return 1;
            }   
        }   
    }

    return 0;
}
    
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

    // Iterate over the characters of the file
    for (int c = fgetc(file); c != EOF; c = fgetc(file))
    {
        newsection = newsec(file, c);

        // Parse section
        if (newsection) 
        {
            c = fgetc(file);
            //printf("%c", c);
            if (index < maxchar)
            {
                line[index] = 'q';
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
