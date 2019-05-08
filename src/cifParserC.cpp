#include <Rcpp.h>
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
int newsec(FILE *file, int c)
{
    // Read first character
    c = fgetc(file);
    // Detect # character or move back file pointer 1 place
    if (c == '#') {
        // Read second character
        c = fgetc(file);
        // Detect space character or move back file pointer 2 places
        if (c == ' ') {
            // Read third character
            c = fgetc(file);
            // Detect new line of move back file pointer 3 places
            if (c == '\n') {
                // Return 1 to be used as bool true
                return 1;
            } else {
                fseek(file, -3, SEEK_CUR);
            }
        } else {
            fseek(file, -2, SEEK_CUR);
        }
    } else {
        fseek(file, -1, SEEK_CUR);
    }

    // Return 0 to be used as bool false
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
    int newsection = 1;

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
        newsection = newsec(file, c);
    }

    // Terminate array
    line[index] = '\0';

    // Close text
    fclose(file);

    // Return output
    return List::create(line);
}
