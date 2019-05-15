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
int newsec(FILE *file, int c)
{
    // Read three characters in array 'line'
    char line[4];
    line[0] = c;
    for (int i = 1; i < 3; i++)
    {
        c = fgetc(file);
        line[i] = c;
    }
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
Rcpp::StringVector entry(FILE *file, int c)
{
    // Read 7 characters in array 'line' to recognize section
    char line[8];
    line[0] = c;
    for (int i = 1; i < 7; i++)
    {
        line[i] = fgetc(file);
    }
    // Terminate array
    line[8] = '\0';

    // Check if section is the one desired
    if (strcmp(line, "_entry.\0") == 0) 
    { //yes
        // Skip unnecesary chars
        fseek(file, 5, SEEK_CUR);

        // Read 4 characters with pdb ID
        for (int i = 0; i < 4; i++)
        {
            line[i] = fgetc(file);
        }
        // Terminate array
        line[4] = '\0';

        // Create Rcpp string vector
        Rcpp::StringVector myvector(1);
        // Assign resulting char string
        myvector[0] = line;
        // Assign names attribute
        myvector.attr("names") = "id";

        // Skip file pointer to next section/line
        int c;
        while ((c = fgetc(file)) != '\n');

        // Return Rcpp string vector
        return myvector;

    } else { //no: 
        // Move file pointer back
        fseek(file, -7, SEEK_CUR);

        // Return empty string
        Rcpp::StringVector myvector(1);
        return myvector;
        //return 1;
    }
}

// Helper function to detect the end of the mmCIF file
int is_end(FILE *file, int *c)
{
    // Move file pointer two positions ahead
    fseek(file, 2, SEEK_CUR);
    int tmp = *c;

    // Read byte and assign to variable 'c' using its pointer
    *c = fgetc(file);

    // Check if the byte was end of file
    if (*c == EOF) 
    { // yes
        return 1;
    } else { // no
        // Move file pointer back
        fseek(file, -2, SEEK_CUR);
        // Change back 'c' variable
        *c = tmp;
        return 0;
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
    Rcpp::StringVector line3;
    Rcpp::StringVector sec1;
    //char line3[maxchar] = "asfasdgfsadfasdf";
    line[0] = line2[0];

//    fseek(file, 10, SEEK_CUR);
//    int c = fgetc(file);
//    int newsection = newsec(file, c);
//    c = fgetc(file);
//    Rcpp::StringVector line3 = entry(file, c);
//    line[0] = fgetc(file);
//    line[1] = '\0';
//
    
    int c = fgetc(file);
    while ((c = fgetc(file)) != EOF && c != '\n');
    int newsection;
    //int newsection = newsec(file);
    // Iterate over the characters of the file
    //for (int c = fgetc(file); c != EOF; c = fgetc(file))
    //c = fgetc(file);
    do {
        c = fgetc(file);
        newsection = newsec(file, c);
        line2[0] = c;
        line2[1] = '\0';
        printf("%s", line2);
        // Parse section
        if (newsection) 
        {
            c = fgetc(file);
            line3 = entry(file, c);
            if (line3[0] != "") 
            {
                sec1 = line3;
            }
            Rcpp::Rcout << line3[0] << '\n';

            line[0] = c;
            for (int i = 0; i < 2; i++)
            {
                line[i] = fgetc(file);
            }
            line[2] = '\n';
            line[3] = '\0';
            printf("%s", line);

            //c = fgetc(file);
            //line2[0] = c;
            //line2[1] = '\0';
            //printf("%s", line2);
            //fseek(file, -1, SEEK_CUR);
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
        //} else {
            //printf();
        }
        while ((c = fgetc(file)) != EOF && c != '\n');
        //c = fgetc(file);
        //if (c == EOF)
        //{
        //    line2[0] = c;
        //    line2[1] = '\0';
        //    printf("%s", line2);
//    } while ((c = fgetc(file)) != EOF);

        // Check for EOF. Second argument is the direciton of 'c' in memory
        is_end(file, &c);
    } while (c != EOF);

    // Terminate array
//    line[index] = '\0';

    // Close text
    fclose(file);

    // Return output
    //return List::create(_["entry"] = line3, line2, line);
    return List::create(_["entry"] = sec1, line2, line);
}
