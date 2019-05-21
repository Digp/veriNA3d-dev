#include <Rcpp.h>
#include <string.h>
#include <stdio.h>
#include "helpers.h"

/*
 * Function can now parse the tow first sections (entry and audit_conform)
 * gzip files are incorrectly read!
 */

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

    // Define containers for sections 
    Rcpp::StringVector presec1;
    Rcpp::StringVector presec2;
    Rcpp::StringVector sec1;
    Rcpp::StringVector sec2;

    int c;
    //int c = fgetc(file);
    int newsection;
    //int newsection = newsec(file);
    do {
        // After newline, check if a new section starts
        c = fgetc(file);
        newsection = newsec(file, c);
        // Parse section
        if (newsection) 
        {
            c = fgetc(file);
            presec1 = entry(file, c);
            if (presec1[0] != "") 
            {
                sec1 = presec1;
            }
            //Rcpp::Rcout << presec1[0];// << '\n';
            c = fgetc(file);
            presec2 = audit_conform(file, c);
            if (presec2[0] != "") 
            {
                sec2 = presec2;
            }
            //Rcpp::Rcout << presec2[0];// << '\n';
        }

        // Iterate over the characters of the file until a newline is found
        while ((c = fgetc(file)) != EOF && c != '\n');

        // Check for EOF. Second argument is the direciton of 'c' in memory
        is_end(file, &c);
    } while (c != EOF);

    // Close text
    fclose(file);

    // Return output
    return List::create(_["entry"] = sec1, _["audit_conform"] = sec2);
}
