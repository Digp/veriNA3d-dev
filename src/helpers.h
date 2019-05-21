#ifndef HELPERS_H
#define HELPERS_H

#include <Rcpp.h>
#include <string.h>
#include <stdio.h>
using namespace Rcpp;


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

    // If line is '# \n'
    if (strcmp(line, "# \n\0") == 0) {
        // Return 1 to be used as bool true
        return 1;
    }

    // Terminate array one before
    line[2] = '\0';

    // If line is '#\n'
    if (strcmp(line, "#\n\0") == 0) {
        // Move file pointer one back
        fseek(file, -1, SEEK_CUR);
        // Return 1 to be used as bool true
        return 1;
    }

    // If it does not detect array '# \n', move back file pointer 3 places
    fseek(file, -3, SEEK_CUR);
    // Return 0 to be used as bool false
    return 0;
}

// Helper function to detect the end of the mmCIF file
int is_end(FILE *file, int *c)
{
    // Define counter
    int tmp = 0;

    do {
        // Check if the byte is the end of file
        if (*c == EOF)
        { // yes
            return 1;
        } //else { // no
        tmp++;

        // Read byte and assign to variable 'c' using its pointer
        *c = fgetc(file);
    } while (tmp < 4);

    // Move file pointer back
    fseek(file, -tmp, SEEK_CUR);
    // Change back 'c' variable2
    //*c = tmp; // Not necessary
    return 0;
}

// Helper function to find and return the entry section of the mmCIF
Rcpp::StringVector entry(FILE *file, int c)
{
    // Read 7 characters in array 'line' to recognize section
    char line[maxchar];
    line[0] = c;
    for (int i = 1; i < 7; i++)
    {
        line[i] = fgetc(file);
    }
    // Terminate array
    line[7] = '\0';

    // Check if section is the one desired
    if (strcmp(line, "_entry.\0") == 0)
    { //yes
        // Skip unnecesary chars
        fseek(file, 2, SEEK_CUR);

        // Read pdbID or string
        int i = 0;
        while ((c = fgetc(file)) != '\n') {
            if (c != ' ')
            {
                line[i] = c;
                i++;
            }
        }
        // Terminate array
        line[i] = '\0';
        //printf("%s", line);

        // Create Rcpp string vector
        Rcpp::StringVector myvector(1);
        // Assign resulting char string
        myvector[0] = line;
        // Assign names attribute
        myvector.attr("names") = "id";

        // Move file pointer one back to stay in same line
        // Calling function needs it this way
        fseek(file, -1, SEEK_CUR);

        // Return Rcpp string vector
        return myvector;

    } else { //no: 
        // Move file pointer back
        fseek(file, -7, SEEK_CUR);

        // Return empty string
        //Rcpp::StringVector myvector(1);
        //return myvector;
        return 1;
    }
}

// Helper function to find and return the entry section of the mmCIF
Rcpp::StringVector audit_conform(FILE *file, int c)
{
    // Read 15 characters in array 'line2' to recognize section
    char line2[maxchar];
    line2[0] = c;
    for (int i = 1; i < 15; i++)
    {
        line2[i] = fgetc(file);
    }
    // Terminate array
    line2[15] = '\0';
    //printf("%s", line2);

    // Check if section is the one desired
    if (strcmp(line2, "_audit_conform.\0") == 0)
    { //yes
        // Create Rcpp string vector
        Rcpp::StringVector myvector2(3);

        // Move file pointer back
        fseek(file, -15, SEEK_CUR);

        int i = 0;
        while (i < 3) {
            // Skip unnecesary chars
            fseek(file, 28, SEEK_CUR);

            // Read pdbID or string
            int j = 0;
            while ((c = fgetc(file)) != '\n') {
                if (c != ' ')
                {
                    line2[j] = c;
                    j++;
                }
            }
            // Terminate array
            line2[j] = '\0';

            // Assign resulting char string
            myvector2[i] = line2;

            i++;
        }

        // Assign names attribute
        myvector2.attr("names") = CharacterVector::create("dict_name", "dict_version", "dict_location");

        // Move file pointer one back to stay in same line
        // Calling function needs it this way
        fseek(file, -1, SEEK_CUR);

        // Return Rcpp string vector
        return myvector2;

    } else { //no: 
        // Move file pointer back
        fseek(file, -15, SEEK_CUR);

        // Return empty string
        return 1;
    }
}

#endif
