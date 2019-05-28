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

// Helper function to read the sections starting with _loop
Rcpp::StringVector core_nonloop(FILE *file, int skip)
{
    // Create Rcpp string vector with unknown length
    Rcpp::StringVector mynames;
    Rcpp::StringVector mycontent;

    // Create necessary variables to contain data
    char line2[maxchar];
    char c = fgetc(file);
    skip = skip - 1;

    // As long as the first character of the line is '_' do:
    do {
        // Skip desired characters
        fseek(file, skip, SEEK_CUR);

        // Parse name
        // Read characters until a space is found
        int i = 0;
        while ((c = fgetc(file)) != ' ') {
            line2[i] = c;
            i++;
        }
        // Terminate array
        if (line2[i-1] == ' ')
        {
            line2[i-1] = '\0';
        } else {
            line2[i] = '\0';
        }

        // Resize Rcpp vector to be able to add a new string
        // Add new string into the Rcpp vector of strings
        mynames.push_back(line2);

        // Skip space characters
        // If current character is ' ', keep reading until finding something else
        while ((c = fgetc(file)) == ' ');

        // Parse content
        // Read characters until the end of line is found
        i = 0;
        char end;

        // If it's a new line, skip character
        if (c == '\n')
        {
            c = fgetc(file);
        }

        // Check if it is a ' or a ", or a ;, or something different.
        // Then, define the character to detect end character
        if (c == '"')
        { // If first character was ", end = "
            // Read next character
            c = fgetc(file);
            end = '"';
        } else if (c == ';') { // Else if first character was ;, end = ';'
            // Read next character
            c = fgetc(file);
            end = ';';
        } else if (c == '\'') { // Else if first character was ', end = '
            // Read next character
            c = fgetc(file);
            end = '\'';
        } else { // Else, just define end = " "
            end = ' ';
        }

        // Keep reading and saving characters until a character matches the variable end
        do {
            // Ignore newlines
            if (c != '\n')
            {
                line2[i] = c;
                i++;
            }
            c = fgetc(file);
        } while (c != end);
        // Terminate array
        if (line2[i-1] == ' ')
        {
            line2[i-1] = '\0';
        } else {
            line2[i] = '\0';
        }

        // Resize Rcpp vector to be able to add a new string
        // Add new string into the Rcpp vector of strings
        mycontent.push_back(line2);
        // If current character is not '\n', keep reading until finding a newline
        while ((c = fgetc(file)) != EOF && c != '\n');

    // The loop will finish when the first character of the line is not '_'
    } while ((c = fgetc(file)) == '_');

    // Move file pointer one back to stay in same line
    // Calling function needs it this way
    fseek(file, -1, SEEK_CUR);

    // Assign attribute names to the data.frame
    mycontent.attr("names") = mynames;

    return mycontent;
}

// Helper function to parse non_loop sections of the mmCIF
Rcpp::StringVector parse_nonloop(FILE *file, int c, char title[maxchar])
{
    // Measure title length
    int len = strlen(title);

    // Read 15 characters in array 'line2' to recognize section
    char line2[maxchar];
    line2[0] = c;
    for (int i = 1; i < len; i++)
    {
        line2[i] = fgetc(file);
    }
    // Terminate array
    line2[len] = '\0';

    // Check if section is the one desired
    if (strcmp(line2, title) == 0)
    { //yes
        //printf("%s\n", line2);
        // Move file pointer back
        fseek(file, -len, SEEK_CUR);

        // Parse _loop section
        Rcpp::StringVector myvector = core_nonloop(file, len);

        // Return data frame
        return myvector;

    } else { //no: 
        // Move file pointer back
        fseek(file, -len, SEEK_CUR);

        // Return empty string
        return 1;
    }
}

// Helper function to read the sections starting with _loop
Rcpp::DataFrame core_loop(FILE *file, int skip)
{
    // Create Rcpp string vector with unknown length
    Rcpp::StringVector myvec;

    // Create necessary variables to contain data
    char line2[maxchar];
    char c;

    // As long as the first character of the line is '_' do:
    do {
        // Skip desired characters
        fseek(file, skip, SEEK_CUR);

        // Read characters until a newline is found
        int i = 0;
        while ((c = fgetc(file)) != '\n') {
            line2[i] = c;
            i++;
        }
        // Terminate array
        if (line2[i-1] == ' ')
        {
            line2[i-1] = '\0';
        } else {
            line2[i] = '\0';
        }

        // Resize Rcpp vector to be able to add a new string
        // Add new string into the Rcpp vector of strings
        myvec.push_back(line2);

    // The loop will finish when the first character of the line is not '_'
    } while ((c = fgetc(file)) == '_');
    //int kk = myvec.size();
    //printf("%i", kk);

    // Count lines
    int maxline = 0;
    double maxchars = 0;
    do {
        do {
            maxchars++;
        } while ((c = fgetc(file)) != '\n');
        maxchars++;
        maxline++;
    } while ((c = fgetc(file)) != '#');
    // Move file pointer back
    fseek(file, -maxchars, SEEK_CUR);

    // Create list of vectors with the length equal to the number of lines
    Rcpp::List tmp(myvec.size());
    for (int w = 0; w < myvec.size(); w++)
    {
        tmp[w] = Rcpp::StringVector(maxline);
    }

    // Define index that will be used in next loop
    int i;
    int j;
    char line[maxchar];
    char end;
    int lineind = 0;

    // As long as the first character of the line is NOT "#" do:
    do {
        // Set to 0 an index that will count the entries in the line
        i = 0;
        // Move file pointer one back
        fseek(file, -1, SEEK_CUR);
        // As long as index < length of Rcpp list
        while (i < myvec.size()) {
            // Read first character
            c = fgetc(file);

            // If it's a new line, skip character
            if (c == '\n')
            {
                c = fgetc(file);
            }

            // Set index to 0
            j = 0;

            // Check if it is a ' or a ", or a ;, or something different.
            // Then, define the character to detect end character
            if (c == '"')
            { // If first character was ", end = "
                // Read next character
                c = fgetc(file);
                end = '"';
            } else if (c == ';') { // Else if first character was ;, end = ';'
                // Read next character
                c = fgetc(file);
                end = ';';
            } else if (c == '\'') { // Else if first character was ', end = '
                // Read next character
                c = fgetc(file);
                end = '\'';
            } else { // Else, just define end = " "
                end = ' ';
            }

            // Keep reading and saving characters until a character matches the variable end
            do {
                // Ignore newlines
                if (c != '\n')
                {
                    line[j] = c;
                    j++;
                }
                c = fgetc(file);
            } while (c != end);
            // Terminate array
            line[j] = '\0';
            //printf("%i", j);
            //printf("%s", line);

            // Add string to Rcpp list, vector selected by variable index 
            as<StringVector>(tmp[i])[lineind] = line;

            // Add +1 to index
            i++;
            // If current character is ' ', keep reading until finding something else
            while ((c = fgetc(file)) == ' ');
            // Move file pointer one back
            fseek(file, -1, SEEK_CUR);
        }
        // If current character is not '\n', keep reading until finding a newline
        while ((c = fgetc(file)) != EOF && c != '\n');

    lineind++;
    // Read next character to check if it's '#'
    } while ((c = fgetc(file)) != '#');

    // Move file pointer one back to stay in same line
    // Calling function needs it this way
    fseek(file, -1, SEEK_CUR);

    // Create a dataframe from the list of vectors
    Rcpp::DataFrame df(tmp);

    // Assign attribute names to the data.frame
    df.attr("names") = myvec;
    return df;
}

// Helper function to parse _loop sections of the mmCIF
Rcpp::DataFrame parse_loop(FILE *file, int c, char title[maxchar])
{
    // Measure title length
    int len = strlen(title);

    // Read 15 characters in array 'line2' to recognize section
    char line2[maxchar];
    line2[0] = c;
    for (int i = 1; i < len; i++)
    {
        line2[i] = fgetc(file);
    }
    // Terminate array
    line2[len] = '\0';

    // Check if section is the one desired
    if (strcmp(line2, title) == 0)
    { //yes
        //printf("%s\n", line2);
        // Move file pointer back
        len = len - 7;
        fseek(file, -len, SEEK_CUR);

        // Parse _loop section
        Rcpp::DataFrame df = core_loop(file, len);

        // Return data frame
        return df;

    } else { //no: 
        // Move file pointer back
        fseek(file, -len, SEEK_CUR);

        // Create empty df
        DataFrame df = DataFrame::create(Named("V1") = "");

        // Return empty string
        return df;
    }
}
#endif
