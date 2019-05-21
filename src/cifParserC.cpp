#include <Rcpp.h>
#include <string.h>
#include <stdio.h>
#include "helpers.h"

/*
 * Function can now detect new sections and return to R the first character
 * of each section.
 * gzip files are incorrectly read!
 */

// [[Rcpp::export]]
List cifParserC(std::string strings="")
{
    // Open file
    FILE *file = fopen(strings.c_str(), "r");
    fseek(file, 0, SEEK_SET);
    // Check if pointer to file is correct, or stop function
    if (file == NULL)
    {
        printf("Could not open %s\n", strings.c_str());
        return 1;
    }

    // Define object line and index to save the contents of the file
    // char* pl = NULL;
    //char line[maxchar];
    //char line2[maxchar] = "asfasdgfsadfasdf";
    Rcpp::StringVector presec1;
    Rcpp::StringVector presec2;
    Rcpp::StringVector sec1;
    Rcpp::StringVector sec2;
    //line[0] = line2[0];

    int c;
    //int c = fgetc(file);
    //while ((c = fgetc(file)) != EOF && c != '\n');
    int newsection;
    //int newsection = newsec(file);
    // Iterate over the characters of the file
    //for (int c = fgetc(file); c != EOF; c = fgetc(file))
    //c = fgetc(file);
    //        line[0] = c;
    //        line[1] = '\n';
    //        line[2] = '\0';
    //        printf("%s", line);
    //        fseek(file, -1, SEEK_CUR);
    do {
        c = fgetc(file);
        newsection = newsec(file, c);
        //line2[0] = c;
        //line2[1] = '\n';
        //line2[1] = '\0';
        //printf("%s", line2);
        //printf("%i", newsection);
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
            //presec1 = entry(file, c);
            //c = fgetc(file);
            presec2 = audit_conform(file, c);
            if (presec2[0] != "") 
            {
                sec2 = presec2;
            }
            //Rcpp::Rcout << presec2[0];// << '\n';

            //for (int i = 0; i < 2; i++)
            //{
            //    line[i] = fgetc(file);
            //}

            //c = fgetc(file);
            //line2[0] = c;
            //line2[1] = '\0';
            //printf("%s", line2);
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
//        c = fgetc(file);
//        line2[0] = c;
//        line2[1] = '\0';
//        printf("%s", line2);
//        fseek(file, -3, SEEK_CUR);

        // Check for EOF. Second argument is the direciton of 'c' in memory
        is_end(file, &c);
    } while (c != EOF);

    // Terminate array
//    line[index] = '\0';

    fseek(file, 0, SEEK_END);
    // Close text
    fclose(file);

    // Return output
    return List::create(_["entry"] = sec1, _["audit_conform"] = sec2);
    //return List::create(_["entry"] = sec1, sec2);
    //return List::create(_["entry"] = sec1);
}
