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

// Helper function to find and return the entry section of the mmCIF
Rcpp::StringVector entry(FILE *file) {
    // Read three characters
    char line[8];

    for (int i = 0; i < 7; i++)
    {
        line[i] = fgetc(file);
    }
    //line[2] = fgetc(file);
    // Terminate array
    line[8] = '\0';

    if (strcmp(line, "_entry.\0") == 0) {
        printf("%s\n", line);
        Rcpp::StringVector myvector(1);
        myvector[0] = line;
        return myvector;
    } else {
        fseek(file, -7, SEEK_CUR);
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

    //int index = 0;

    fseek(file, 9, SEEK_CUR);
    line2[0] = fgetc(file);
    line2[1] = '\0';
    //printf("%i\n", c);
    int newsection = newsec(file);
    printf("%i\n", newsection);
    Rcpp::StringVector line3 = entry(file);

    for (int i = 0; i < 7; i++)
    {
        line[i] = fgetc(file);
    }
    line[7] = '\0';
    

    // Iterate over the characters of the file
    //for (int c = fgetc(file); c != EOF; c = fgetc(file))
    //do {

        // Parse section
//        if (newsection) 
//        {
//            line[index] = c;
//            index++;
//            while ((c = fgetc(file)) != EOF && c != '\n')
//            {
//                //c = fgetc(file);
//                //printf("%c", c);
//                if (index < maxchar)
//                {
//                    line[index] = c;
//                    index++;
//                }
//            }
//        }
//
//        // After newline, check if a new section starts
//        newsection = newsec(file);
    //}
    //    c = fgetc(file);
    //} while (c != EOF);

    // Terminate array
//    line[index] = '\0';

    // Close text
    fclose(file);

    // Return output
    return List::create(_["line"] = line, line2, line3);
}
