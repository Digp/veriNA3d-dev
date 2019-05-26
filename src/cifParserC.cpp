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
    Rcpp::StringVector tmpsec;
    Rcpp::DataFrame tmpsec_df;
    Rcpp::StringVector sec1;
    Rcpp::StringVector sec2;
    Rcpp::DataFrame sec3;
    Rcpp::DataFrame sec5;

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
            // Check if it's "_entry" section and parse it
            c = fgetc(file);
            tmpsec = entry(file, c);
            if (tmpsec[0] != "") 
            {
                sec1 = tmpsec;
            }
            //Rcpp::Rcout << tmpsec[0];// << '\n';

            // Check if it's "_audit_conform" section and parse it
            c = fgetc(file);
            tmpsec = audit_conform(file, c);
            if (tmpsec[0] != "") 
            {
                sec2 = tmpsec;
            }

            // Check if it's "_database_2" section and parse it
            c = fgetc(file);
            tmpsec_df = database_2(file, c);
            if (tmpsec_df.size() > 1)
            {
                sec3 = tmpsec_df;
            }

            // Check if it's "_pdbx_database_status" section and parse it
            // Check if it's "_audit_author" section and parse it
            c = fgetc(file);
            tmpsec_df = audit_author(file, c);
            if (tmpsec_df.size() > 1)
            {
                sec5 = tmpsec_df;
            }

            // Check if it's "_entity" section and parse it
            // Check if it's "_chem_comp" section and parse it
            // Check if it's "_exptl" section and parse it
            // Check if it's "_struct" section and parse it
            // Check if it's "_struct_keywords" section and parse it
            // Check if it's "_struct_asym" section and parse it
            // Check if it's "_atom_sites" section and parse it
            // Check if it's "_atom_type" section and parse it
            // Check if it's "_atom_site" section and parse it
        }

        // Iterate over the characters of the file until a newline is found
        while ((c = fgetc(file)) != EOF && c != '\n');

        // Check for EOF. Second argument is the direciton of 'c' in memory
        is_end(file, &c);
    } while (c != EOF);

    // Close text
    fclose(file);

    // Return output
    return List::create(_["entry"] = sec1, _["audit_conform"] = sec2, _["database_2"] = sec3,
                        _["audit_author"] = sec5);
}
