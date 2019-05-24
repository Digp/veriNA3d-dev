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

// Helper function to find and return the entry section of the mmCIF
Rcpp::DataFrame database_2(FILE *file, int c)
{
    // Read 15 characters in array 'line2' to recognize section
    char line2[maxchar];
    line2[0] = c;
    for (int i = 1; i < 18; i++)
    {
        line2[i] = fgetc(file);
    }
    // Terminate array
    line2[18] = '\0';
    //printf("%s", line2);

    // Check if section is the one desired
    if (strcmp(line2, "loop_\n_database_2.\0") == 0)
    { //yes
        // Move file pointer back
        fseek(file, -11, SEEK_CUR);

        // Create Rcpp string vector with unknown length
        Rcpp::StringVector myvec;

        // As long as the first character of the line is '_' do:
        do {
            // Skip 12 characters
            fseek(file, 11, SEEK_CUR);

            // Read characters until a newline is found
            int i = 0;
            while ((c = fgetc(file)) != '\n') {
                line2[i] = c;
                i++;
            }
            line2[i-1] = '\0';

            // Resize Rcpp vector to be able to add a new string
            // Add new string into the Rcpp vector of strings
            myvec.push_back(line2);

        // The loop will finish when the first character of the line is not '_'
        } while ((c = fgetc(file)) == '_');

        // Move file pointer one back - NOT NECESSARY SINCE IT'S DONE LATER

        // Create list of vectors of unknown length or somewhat similar
        //Rcpp::List tmp(myvec.size());
        // Define index that will be used in next loop

        // As long as the first character of the line is NOT "#" do:
            // Set to 0 an index that will count the entries in the line
            // Move file pointer one back
            // As long as index < length of Rcpp list
                // Read first character

                // Check if it is a ' or a ", or a ;, or something different.

                // Define the character to detect end of entry:
                // If first character was ', end = '
                // Else if first character was ", end = "
                // Else if first character was ;, end = ';'
                // Else, Save first character of line as first character of vector and define end = " "

                // Keep reading and saving characters until a character matches the variable end, and ignoring newlines

                // Add string to Rcpp list, vector selected by variable index 
                // Add +1 to index

            // If current character is not '\n', keep reading until finding a newline
            // Read next character to check if it's '#'

        // Creating vector v
        NumericVector v = {1,2};
        // Creating DataFrame df
        DataFrame df = DataFrame::create(Named("V1") = v,         // simple assign
                                         Named("V2") = clone(v)); // using clone()
        //df.attr("names") = CharacterVector::create("dict_name", "dict_version");
        df.attr("names") = myvec;
        // Changing vector v
        v = v * 2;
        return df;

    } else { //no: 
        // Move file pointer back
        fseek(file, -18, SEEK_CUR);

        DataFrame df = DataFrame::create(Named("V1") = "");         // simple assign

        // Return empty string
        return df;
    }
}

#endif
