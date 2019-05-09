#include <Rcpp.h>
#include <string.h>
#include <stdio.h>
using namespace Rcpp;

/*
 * Function can now detect new sections and return to R the first character
 * of each section.
 * gzip files are incorrectly read!
 */

// Temporal macro definition
#define maxchar 10000

// Helper function to detect new mmCIF section based on lines like '# \n'
// Designed to be executed after a \n character or beggining of file
int newsec(FILE *file)
{
    // Read three characters
    char line[4];
    line[0] = fgetc(file);
    line[1] = fgetc(file);
    line[2] = fgetc(file);
    // Terminate array
    line[3] = '\0';

    // If it does not detect array '# \n', move back file pointer 3 places
    if (strcmp(line, "# \n\0")) {
        fseek(file, -3, SEEK_CUR);
        // Return 0 to be used as bool false
        return 0;
    } else {
        // Return 1 to be used as bool true
        return 1;
    }
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
    // char* pl = NULL;
    char line[maxchar];
    char line2[maxchar] = "asfasdgfsadfasdf";
    // pl = &line2;
    //*pl = 'b';
    //free(pl);

    int index = 0;
    int newsection = newsec(file);

    // Iterate over the characters of the file
    for (int c = fgetc(file); c != EOF; c = fgetc(file))
    {

        // Parse section
        if (newsection) 
        {
            line[index] = c;
            index++;
            while ((c = fgetc(file)) != EOF && c != '\n')
            {
                //c = fgetc(file);
                //printf("%c", c);
                if (index < maxchar)
                {
                    line[index] = c;
                    index++;
                }
            }
        }

        // After newline, check if a new section starts
        newsection = newsec(file);
    }

    // Terminate array
    line[index] = '\0';

    // Close text
    fclose(file);

    // Return output
    return List::create(_["line"] = line, line2);
}
