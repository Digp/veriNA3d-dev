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
    Rcpp::DataFrame sec7;
    Rcpp::DataFrame sec11;
    Rcpp::DataFrame sec13;
    Rcpp::DataFrame sec14;

    int c;
    //int c = fgetc(file);
    int newsection;
    //int newsection = newsec(file);
    do {
        // After newline, check if a new section starts
        parse_new:
        c = fgetc(file);
        newsection = newsec(file, c);

        // Parse section
        if (newsection) 
        {
            // Check if it's "_entry" section and parse it
            c = fgetc(file);
            char title1[] = "_entry.\0";
            //tmpsec = entry(file, c);
            tmpsec = parse_nonloop(file, c, title1);
            if (tmpsec[0] != "") 
            {
                sec1 = tmpsec;
                goto parse_new;
            }
            //Rcpp::Rcout << tmpsec[0];// << '\n';

            // Check if it's "_audit_conform" section and parse it
            c = fgetc(file);
            tmpsec = audit_conform(file, c);
            if (tmpsec[0] != "") 
            {
                sec2 = tmpsec;
                goto parse_new;
            }

            // Check if it's "_database_2" section and parse it
            c = fgetc(file);
            char title3[] = "loop_\n_database_2.\0";
            tmpsec_df = parse_loop(file, c, 18, title3);
            //tmpsec_df = database_2(file, c);
            if (tmpsec_df.size() > 1 || tmpsec_df.nrows() > 1)
            {
                sec3 = tmpsec_df;
                goto parse_new;
            }

            // Check if it's "_pdbx_database_status" section and parse it
            // Check if it's "_audit_author" section and parse it
            c = fgetc(file);
            char title5[] = "loop_\n_audit_author.\0";
            tmpsec_df = parse_loop(file, c, 20, title5);
            if (tmpsec_df.size() > 1 || tmpsec_df.nrows() > 1)
            {
                sec5 = tmpsec_df;
                goto parse_new;
            }

            // Check if it's "_entity" section and parse it
            // Check if it's "_chem_comp" section and parse it
            c = fgetc(file);
            char title7[] = "loop_\n_chem_comp.\0";
            tmpsec_df = parse_loop(file, c, 17, title7);
            if (tmpsec_df.size() > 1 || tmpsec_df.nrows() > 1)
            {
                sec7 = tmpsec_df;
                goto parse_new;
            }

            // Check if it's "_exptl" section and parse it
            // Check if it's "_struct" section and parse it
            // Check if it's "_struct_keywords" section and parse it
            // Check if it's "_struct_asym" section and parse it
            c = fgetc(file);
            char title11[] = "loop_\n_struct_asym.\0";
            tmpsec_df = parse_loop(file, c, 19, title11);
            if (tmpsec_df.size() > 1 || tmpsec_df.nrows() > 1)
            {
                sec11 = tmpsec_df;
                goto parse_new;
            }

            // Check if it's "_atom_sites" section and parse it
            // Check if it's "_atom_type" section and parse it
            c = fgetc(file);
            char title13[] = "loop_\n_atom_type.\0";
            tmpsec_df = parse_loop(file, c, 17, title13);
            if (tmpsec_df.size() > 1 || tmpsec_df.nrows() > 1)
            {
                sec13 = tmpsec_df;
                goto parse_new;
            }

            // Check if it's "_atom_site" section and parse it
            c = fgetc(file);
            char title14[] = "loop_\n_atom_site.\0";
            tmpsec_df = parse_loop(file, c, 17, title14);
            if (tmpsec_df.size() > 1 || tmpsec_df.nrows() > 1)
            {
                sec14 = tmpsec_df;
                goto parse_new;
            }

        }

        // Iterate over the characters of the file until a newline is found
        while ((c = fgetc(file)) != '\n');

        // Check for EOF. Second argument is the direciton of 'c' in memory
        is_end(file, &c);
    } while (c != EOF);

    // Close text
    fclose(file);

    // Return output
    return List::create(_["entry"] = sec1, _["audit_conform"] = sec2, _["database_2"] = sec3,
                        _["audit_author"] = sec5, _["chem_comp"] = sec7,
                        _["struct_asym"] = sec11, _["atom_type"] = sec13, _["atom_site"] = sec14);
}
